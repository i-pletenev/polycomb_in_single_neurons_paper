import numpy as np
import pandas as pd
import cooler
import cooltools
from cooltools.api.snipping import ObsExpSnipper
from coolpuppy.lib.numutils import get_enrichment
import bioframe as bf
from scipy.stats import chi2_contingency

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib import ticker
from matplotlib.ticker import EngFormatter
import matplotlib.gridspec as gridspec
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from chromosome_plotting.chromosome_plotting import Chromosome


def bedpe_to_bed(df):
    """
    Extract unique bed regions from bedpe dataframe. Index is ignored.
    """
    df_reg1 = df[['chrom1', 'start1', 'end1']]
    df_reg2 = df[['chrom2', 'start2', 'end2']]

    df_reg1.columns = df_reg1.columns.str.rstrip('1')
    df_reg2.columns = df_reg2.columns.str.rstrip('2')

    df_bed = pd.concat([df_reg1, df_reg2], ignore_index=True)\
               .drop_duplicates(ignore_index=True)
    
    return df_bed   


def regs_to_bins(regs_df, binsize, mode='middle', return_index=True):
    """
    Assign bins to each genomic region
    
    Parameters
    ----------
    regs_df : pd.DataFrame
        Dataframe with genomic regions. 
        Required columns are: 'chrom', 'start', 'end'
    binsize : int
    mode : 'middle' or 'all'
        If 'middle': for each region, get a bin that overlap the middle of the region.
        If 'all': for each region, get all bins that overlap this region and write to separate rows 
    return_index : bool
        If True and mode='all', return index of input dataframe.
        
    Returns
    -------
    bins_df : pd.DataFrame
        Dataframe with bin coordinates written in
        'chrom', 'start', 'end' columns
    """
    
    if mode == 'middle':
        bins_df = regs_df.copy(deep=True)
        bins_df.loc[:, 'start'] = (bins_df['start'] + bins_df['end']) / 2 \
        - ((bins_df['start'] + bins_df['end']) / 2) % binsize
        bins_df.loc[:, 'end'] = bins_df['start'] + binsize
        bins_df.loc[:, ['start', 'end']] = bins_df.loc[:, ['start', 'end']].astype(int)
        
    elif mode == 'all':
        regs_cp = regs_df.copy(deep=True)
        regs_cp.loc[:, 'first_bin_start'] = regs_cp['start'] - (regs_cp['start'] % binsize)
        regs_cp.loc[:, 'last_bin_end'] = regs_cp['end'] - 1 - (regs_cp['end'] - 1) % binsize + binsize
        regs_cp.loc[:, 'nbins'] = (regs_cp['last_bin_end'] - regs_cp['first_bin_start']) // binsize
        regs_cp.loc[:, 'idx'] = regs_cp.index
        
        all_bins = np.empty((regs_cp['nbins'].sum(), 3), dtype=int)
        i = 0
        for reg_id in range(regs_cp.shape[0]):
            n, start, end, idx = regs_cp.loc[
                reg_id, ['nbins', 'first_bin_start', 'last_bin_end', 'idx']]
            all_bins[i:i+n, 0] = np.arange(start, end, binsize)
            all_bins[i:i+n, 1] = np.arange(start+binsize, end+binsize, binsize)
            all_bins[i:i+n, 2] = np.tile(idx, n)
            i+=n
            
        bins_df = pd.DataFrame(data=all_bins, columns=['start', 'end', 'idx'])
        regs_cp_drop_cols = ['start', 'end', 'first_bin_start', 'last_bin_end', 'nbins']
        bins_df = bins_df.merge(regs_cp.drop(columns=regs_cp_drop_cols), on='idx', how='left')
        bins_df.insert(0, 'chrom', bins_df.pop('chrom'))
        
        if not return_index:
            bins_df = bins_df.drop(columns='idx')
        else:
            bins_df.insert(bins_df.shape[1]-1, 'idx', bins_df.pop('idx'))
        
    else:
        raise ValueError()
        
    return bins_df


def reg_obs_exp(clr, reg, exp):
    """
    Take all pixels that correspond to interactions within a region 
    and return DataFrame with balanced counts normalized by expected. 
    For example, if region is "chr1:10000-40000" and bin size is 10_000, 
    then the following triangle from a Hi-C matrix is considered:
    
          10000 20000 30000 40000
    10000   +     +     +     +
    20000         +     +     +
    30000               +     +
    40000                     +

    Parameters
    ----------
    clr : cooler
        Cooler object to fetch data from
    reg : str
        UCSC string - "chrom:start-end"
    exp : pd.DataFrame
        output of cooltools.expected_cis()

    Returns
    -------
    obs_exp_df : pd.DataFrame
        DataFrame with observed normalized by expected
    """
    chrom = reg.split(':')[0]

    # Get pixels from cooler
    pix = clr.pixels().fetch(reg)
    bins = clr.bins().fetch(reg)
    # pix_ann = cooler.annotate(pix, bins[['weight']], replace=False)
    pix_ann = cooler.annotate(pix, bins, replace=False)
    max_bin = pix_ann['bin1_id'].max()
    pix_ann = pix_ann.loc[pix_ann['bin2_id'] <= max_bin]
    pix_ann['bal_count'] = pix_ann['count'] * pix_ann['weight1'] * pix_ann['weight2']
    pix_ann['dist'] = pix_ann['bin2_id'] - pix_ann['bin1_id']

    # Merge pixels with expected
    exp_chr_df = exp.loc[exp['region1'] == chrom, ['balanced.avg', 'dist']]
    exp_chr_ser = pd.Series(exp_chr_df['balanced.avg'].values, index=exp_chr_df['dist'])
    pix_ann['exp'] = pix_ann['dist'].map(exp_chr_ser)
    # pix_ann.dropna(inplace=True)
    pix_ann['bal_obs_exp'] = pix_ann['bal_count'] / pix_ann['exp']

    return pix_ann.reset_index(drop=True)


def obs_for_bin_pairs(bin_pairs, clr, obs_colname='obs', balance=True):
    """
    Get observed Hi-C values for pairs of bins
    
    Parameters
    ----------
    bin_pairs : pd.DataFrame
        Dataframe with pairs of genomic regions. 
        Required columns are: 
        'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'.
    clr : cooler.Cooler
        cooler to extract Hi-C counts
    obs_colname : str
        Name of column to store Hi-C values
    balance : bool
        Passed to 'clr.matrix()'
        
    Returns
    -------
    bin_pairs_obs : pd.DataFrame
        Input dataframe with new column 'obs'
    """
    bin_pairs_obs = bin_pairs
   
    mtx = clr.matrix(balance=balance)
    bin_pairs_obs[obs_colname] = bin_pairs_obs.apply(
        lambda crds: mtx.fetch([crds['chrom1'], crds['start1'], crds['end1']],
                               [crds['chrom2'], crds['start2'], crds['end2']])[0][0],
        axis=1)
    
    return bin_pairs_obs

############################################################################
# Functions related to plotting individual PcG contacts with the same anchor
############################################################################

def plot_gene_dots(gene, hand_anch, stat_gb, savefig=False):
    """
    For a given annotation of PcG contacts, plot snippets of contact
    matrix for a particular anchor ("gene") with all other anchors on the same chromosome.

    Parameters
    ----------
    gene : str
        gene in "genes" column of stat_gb that defines an anchor of
        PcG contacts
    hand_anch : pd.DataFrame
        List of PcG contact anchors with corresponding genes
    stat_gb : pd.DataFrame
        Dataframe with values to plot
    savefig : str or False
        Name of the file to save the figure
    """
    hg38_chromsizes = bf.fetch_chromsizes('hg38')
    hg38_cens = bf.fetch_centromeres('hg38')\
        .set_index('chrom')\
        .rename(columns={'mid': 'cent'})\
        .loc[:, 'cent']

    plot_df = stat_gb.set_index('genes').loc[gene]
    
    if (plot_df.shape[0] > 1) and (isinstance(plot_df, pd.DataFrame)):
        print(f"Gene {gene} has more than one dot, selecting the first dot")
        plot_df = plot_df.iloc[0]
        
    vals_idx = [idx for idx in plot_df.index if idx.startswith('vals_')]
    ndots = len(plot_df[vals_idx[0]])
    target_genes = hand_anch.loc[plot_df['idx2'], 'genes'].str.split(',').str[0].values
    
    wr = [1] * ndots
    wr.append(0.2)
    fig, axs = plt.subplots(len(vals_idx)+1, ndots+1, dpi=200, figsize=[1 + 0.65*ndots, 0.9 * len(vals_idx)],
                            gridspec_kw={'width_ratios': wr})
    for ax in axs[:, ndots]:
        ax.remove()
    
    # Plot chromosome
    for ax in axs[0, :ndots]:
        ax.remove()
    gs = axs[0, 0].get_gridspec()
    ax_chr = fig.add_subplot(gs[0, :-1])
    loci = hand_anch.loc[plot_df['idx2'], 'start'].sort_values()
    chrom = hand_anch.loc[plot_df['idx2'], 'chrom'].iloc[0]
    chrom_plot = Chromosome(length=hg38_chromsizes[chrom], 
                            name=chrom, 
                            centromere=hg38_cens[chrom],
                            loci=loci)
    chrom_plot.plot_chromosome(ax=ax_chr, height=3, linewidth=0.4)
    ax_chr.text(0.95, 0.95, gene, fontsize=8, transform = ax.transAxes)
    
    # Plot snips
    for i, idx in enumerate(vals_idx, start=1):
        arr = np.array(plot_df[idx])

        axs[i, 0].set_ylabel(idx[5:], rotation='horizontal', fontsize=6, ha='right', va='center')
            
        for j in range(ndots):
            ax = axs[i, j]
            pup = ax.imshow(
                arr[j, :, :], 
                cmap='coolwarm',
                norm=LogNorm(vmax=5, vmin=1/5)
            )
            ax.set_aspect(1)
            ax.set_xticks([])
            ax.set_yticks([])
            
            if i == len(vals_idx):
                ax.set_xlabel(target_genes[j], rotation=15, fontsize=6)
    
    ax_cb = fig.add_subplot(gs[1:, ndots])          
    fig.colorbar(pup, cax=ax_cb, shrink=0.25)

    if savefig:
        plt.savefig(savefig)

    return
    
#####################################################################
# Functions related to plotting dots with top and bottom PCA loadings
#####################################################################

def parse_pca_df(pca_loads_path, hand_anch):
    pca_load = pd.read_csv(pca_loads_path, index_col=0)
    pca_load[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']] = \
        pca_load.index.astype(str).to_series()\
        .str.replace(r"\)|\(", "", regex=True)\
        .str.split(', ', expand=True)
    pca_load[['start1', 'end1', 'start2', 'end2']] = \
        pca_load[['start1', 'end1', 'start2', 'end2']].astype(int)
    pca_load['chrom1'] = pca_load['chrom1'].str.replace("'", "")
    pca_load['chrom2'] = pca_load['chrom2'].str.replace("'", "")
    del_chroms = ['chrX', 'chrY', 'chrM']
    pca_load = pca_load.loc[~pca_load['chrom1'].isin(del_chroms) & 
                            ~pca_load['chrom2'].isin(del_chroms)]

    if 'genes' not in pca_load.columns:
        pca_load = pca_load.merge(hand_anch, how='left', left_on=['chrom1', 'start1', 'end1'], 
                                  right_on=['chrom', 'start', 'end'])\
                           .drop(columns=['chrom', 'start', 'end'])\
                           .rename(columns={'genes': 'gene1'})\
                           .merge(hand_anch, how='left', left_on=['chrom2', 'start2', 'end2'], 
                                  right_on=['chrom', 'start', 'end'])\
                           .drop(columns=['chrom', 'start', 'end'])\
                           .rename(columns={'genes': 'gene2'})
        pca_load['genes'] = pca_load['gene1'] + ' <-> ' + pca_load['gene2']
        
    return pca_load


def ucsc2tup(reg):
    """
    Convert UCSC-like coordinates to tuple
    Ex. chr1:10-20 -> ('chr1', 10, 20)
    """
    chrom = reg.split(':')[0]
    start, end = reg.split(':')[1].split('-')
    start, end = int(start), int(end)
    return (chrom, start, end)


def subt_snips(regs_df, clrs, cvd, view_df, nbin_flank=6, print_ct=True):
    """
    Get obs/exp snippets from Hi-C map for each cell types

    Parameters
    ----------
    regs_df : pd.DataFrame
        Coordinates of loops
    clrs : dict
        Dictionary with cell types as keys and coolers as values
    cvd : dict
        Dictionary with cell types as keys and expected values 
        obtained from cooltools.expected_cis()
    view_df : pd.DataFrame
        Dataframe used for cooltools and coolpuppy 
    nbin_flank : int
        Number of bins to flank from each side

    Returns
    -------
    snips : dict
        Dictionary with cell types as keys and lists with numpy 
        matrices as values
    """
    cts = list(clrs.keys())
    snips = {ct: [] for ct in cts}
    chroms = regs_df['chrom1'].unique()

    for i, ct in enumerate(cts):
        if print_ct:
            print(ct)
        oes = ObsExpSnipper(clrs[ct], cvd[ct], view_df=view_df)
        res = clrs[ct].binsize
        for chrom in chroms:
            mtx = oes.select(chrom, chrom)
            regs_chrom = regs_df.loc[regs_df['chrom1'] == chrom]
            for idx, reg in regs_chrom.iterrows():
                flank = nbin_flank * res
                # get middle bin
                mid1 = (reg.start1 + reg.end1) / 2
                mid2 = (reg.start2 + reg.end2) / 2
                start1 = mid1 - (mid1 % res) - flank
                end1 = mid1 - (mid1 % res) + res + flank
                start2 = mid2 - (mid2 % res) - flank
                end2 = mid2 - (mid2 % res) + res + flank
                snip = oes.snip(
                    mtx, chrom, chrom, 
                    tup=(int(start1), int(end1), int(start2), int(end2))
                )
                snips[ct].append(snip) 

    return snips


def plot_loads(regs_df, pl_grps, snips, figsize=None, vmax=20):
    if figsize == None:
        figsize = [0.2*regs_df.shape[0], 0.2*len(pl_grps)]
        
    fig, axs = plt.subplots(len(pl_grps), regs_df.shape[0], dpi=200, 
                            figsize=figsize)
    
    for i, grp in enumerate(pl_grps):
        for j in range(regs_df.shape[0]):
            snip_idx = np.argsort(regs_df.index)[j]
        
            ax = axs[i, j]
            ax.imshow(
                snips[grp][snip_idx], 
                cmap='coolwarm',
                norm=LogNorm(vmax=vmax, vmin=1/vmax)
            )
            ax.set_xticks([])
            ax.set_yticks([])
    
            if j == 0:
                ax.set_ylabel(grp, fontsize=6, rotation=0, ha='right', va='center')
    
            if i == len(pl_grps) - 1:
                genes = regs_df.loc[j, "genes"]
                ax.set_xlabel(genes, fontsize=6, rotation=90)
    
    plt.subplots_adjust(hspace=0, wspace=0)


def format_ticks_zoom_out(ax, base=2e7, x=True, y=True, rotate=True):
    bp_formatter = ticker.EngFormatter('b')
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base))
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base))
        ax.xaxis.tick_bottom()


def plot_hic_reg_zoom_out(clr, reg, ax, vmin, vmax, snips=None, 
                          title=None, cbar=True, cmap='default'):
    """
    Plot region of a Hi-C heatmap. 
    Optionally, contour parts of a heatmap (snips) with rectangles.
    
    Parameters
    ----------
    clr : cooler
        Cooler object to fetch data from
    reg : list or tuple
        List or tuple with chromosome name, region start 
        and region end of the form (chrom, start, end)
    ax : matplotlib Axes
        Axes object to draw the plot onto
    vmin, vmax, cmap
        Passed on to `ax.matshow()`
    snips : list
        List of snips; each snip is a list itself
        conaining two region coordinates of the form 
        (chrom, start, end)
    title : str
        Title for the plot
    cbar : bool
        Whether to draw a colorbar
    cmap : str
        cmap to draw
    """
    # Create custom colormap
    if cmap == 'default':
        cmap_colors = np.array([
            [255, 255, 255, 255],
            [245, 166, 35, 255],
            [208, 2, 27, 255],
            [0, 0, 0, 255]
        ]) / 255
        cmap_nodes = [0, 0.35, 0.6, 1]
        cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(cmap_nodes, cmap_colors)))
    
    # Plot heatmap
    chrom, start, end = reg
    im = ax.matshow(
        clr.matrix(balance=True).fetch(reg),
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=cmap,
        extent=(start, end, end, start), 
        interpolation=None,
    );
    
    # Set title
    if title:
        ax.set_title(title)
    
    format_ticks_zoom_out(ax)
    
    # Draw rectangles
    if snips:
        for j, (reg1, reg2) in enumerate(snips):
            rect = patches.Rectangle(
                (reg2[1], reg1[1]), reg2[2]-reg2[1], reg1[2]-reg1[1], 
                edgecolor='k', facecolor='none', linewidth=1
            )
            ax.add_patch(rect)
    
    # Draw colorbar
    if cbar:
        cb = plt.colorbar(im, ax=ax, shrink=0.5, location='right')
        cb.ax.set_ylabel('Contact probability', rotation=270, fontsize=8)
    
    return im


def plot_av_dot_2d(pup: dict, vmax: float, **subplot_kwargs):
    bp_formatter = EngFormatter('b')

    # get cts from pup
    grps1 = list(pup.keys())
    grps2 = list(pup[grps1[0]].keys())

    nrows, ncols = len(grps1), len(grps2)
    width_ratios = [1] * ncols
    width_ratios.append(0.05)

    subplot_kwargs.setdefault('dpi', 400)
    subplot_kwargs.setdefault('figsize', (ncols / 2, nrows / 2))
    
    fig, axs = plt.subplots(nrows, ncols+1, gridspec_kw={'width_ratios': width_ratios}, **subplot_kwargs)
    
    # vmax = {'short': 6.5, 'long': 2, 'trans': 8}
    vmin = 1/vmax
    # vmin = {k: 1/val for k, val in vmax.items()}
    
    for i, grp1 in enumerate(grps1):
        for j, grp2 in enumerate(grps2):
            ax = axs[i, j] if nrows > 1 else axs[j]
                
            arr = np.nanmean(np.array(pup[grp1][grp2]['data'].to_list()), axis=0)
            im = ax.imshow(
                arr, 
                cmap='coolwarm',
                norm=LogNorm(vmax=vmax, vmin=vmin),
                extent=(-150000, 150000, 150000, -150000)
            )
            
            # Enrichment
            en = get_enrichment(arr, 0)
            ax.text(0.02, 0.98, round(en, 2), va='top', 
                    transform=ax.transAxes, fontsize=6)
            
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            
            if i==0:
                ax.set_title(grp2, fontsize=6)
        
        # Colorbar
        cbar_ax = axs[i, ncols] if nrows > 1 else axs[ncols]
        plt.colorbar(im, cax=cbar_ax, ticks=[vmin, 1, vmax], 
                     format=ticker.FuncFormatter(lambda x, pos: f"{x:.2g}"))
        cbar_ax.set_box_aspect(20)
        cbar_ax.minorticks_off()

    if nrows == 1:
        axs[0].set_ylabel(grps1[0], fontsize=8)
    else:
        for i, grp1 in enumerate(grps1): 
            axs[i, 0].set_ylabel(grp1, fontsize=8)

    return


class Pb_de:
    def __init__(self, pb_df, hand_pc):
        self.pb_df = pb_df
        self.results_df = {}
        self.cont_tab = {}
        self.cont_chisq = {}
        self.hand_pc = hand_pc
        

    def age_deseq2(self, age):
        # Get samples
        age2samp = {
            '2T': ["ga22", "ga24"],
            'fetal': ["ga22", "ga24", "ga34"],
            'infant': ["118d", '179d'],
            'adult': ['20yr', '25yr']
        }
        samps = age2samp[age]
        pb_cols = [samp + '_' + ct for samp in samps for ct in ['PN', 'IN']] 
        pb_cts = np.array([ct for samp in samps for ct in ['PN', 'IN']])
        counts_df = self.pb_df[pb_cols].T
        metadata = pd.DataFrame(data=pb_cts.T, index=pb_cols, columns=['ct'])

        # Run pyDESEQ2
        inference = DefaultInference(n_cpus=5)
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design="~ct",  # compare samples based on the "condition"
            # column ("B" vs "A")
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()
        ds = DeseqStats(
            dds,
            contrast=['ct', 'PN', 'IN'],
            alpha=0.05,
            cooks_filter=True,
            independent_filter=False,
        )
        ds.summary()

        ds.results_df.loc[:, 'hand'] = ds.results_df.index.isin(self.hand_pc)
        ds.results_df.loc[:, 'min_log10_padj'] = -np.log10(ds.results_df['padj'])
        # ds.results_df.loc[:, 'lfc_sign'] = ds.results_df['log2FoldChange'] > 0
        ds.results_df.loc[:, 'lfc_sign'] = ds.results_df['log2FoldChange']\
            .apply(lambda lfc: 1 if lfc > 0 else -1)
        
        self.results_df[age] = ds.results_df


    def cont_table(self, age):
        self.cont_tab[age] = self.results_df[age]\
            .loc[self.results_df[age]['padj'] < 0.05, ['lfc_sign', 'hand']]\
            .value_counts()\
            .to_frame()\
            .reset_index()\
            .pivot(index='lfc_sign', columns='hand')\
            .to_numpy()
        
        self.cont_chisq[age] = chi2_contingency(self.cont_tab[age]).pvalue


    def plot_cont(self, age):
        plt.figure(figsize=[1.5, 1.5], dpi=150)
        sns.heatmap(self.cont_tab[age], square=True, annot=True, fmt='g', cbar=False, cmap='Reds')
        plt.xticks([0.5, 1.5], ['Not at Polycomb dot', 'At Polycomb dot'], rotation=30)
        plt.yticks([0.5, 1.5], ['More expressed in IN', 'More expressed in EN'], rotation=30)
        plt.title(f'DEGs at {age}, p={round(self.cont_chisq[age], 2)}')


    def plot_volc(self, age):
        plt.figure(figsize=[2, 2], dpi=150)
        plot_df = self.results_df[age].loc[self.results_df[age]['padj'] < 0.05].sort_values('hand')
        sns.scatterplot(data=plot_df, 
                        x='log2FoldChange', y='min_log10_padj', hue='hand', 
                        s=10, palette={True: 'red', False: 'grey'})
        plt.title(age)