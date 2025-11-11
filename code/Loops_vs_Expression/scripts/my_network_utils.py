import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import networkx as nx


my_pal = {"EN": "#3924B1",
          "IN": "#C1118C",
          "NN": "#ffc875"}

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
       'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chrX'] + ['trans']
chroms_colors = {chroms[i] : plt.get_cmap('tab20b')((i+1) / len(chroms)) for i in range(len(chroms[:-1])) }
chroms_colors['trans'] = (0,0,0)


def print_network_characteristics(G):
    print('-- General network characteristics --')

    clusters = list(nx.connected_components(G))
    clustering = list(nx.clustering(G).values())
    degrees=[d for n, d in G.degree()]
    
    print('Nodes N:', G.number_of_nodes())
    print('Edges E:', G.number_of_edges())
    print('Density rho:', nx.density(G))
    print('Av. degree <k>:', np.mean(degrees))
    print('N conected components NCC:', nx.number_connected_components(G))
    
    diameter_in_components = []
    for C in (G.subgraph(c).copy() for c in clusters):
        diameter_in_components.append(nx.diameter(C))
    
    print('Max diameter d:',np.max(diameter_in_components))
    
    number_of_nodes = []
    for C in (G.subgraph(c).copy() for c in clusters):
        number_of_nodes.append(C.number_of_edges())
    print('Max number of nodes :',np.max(number_of_nodes))

def plot_heatmap_for_clusters(df_for_heatmap, dpi=120, mins_maxs=None, totals=True, add_chroms=True, rw=0.6):
    nrows, ncols = df_for_heatmap.shape
    fig = plt.figure(figsize=(14, rw*(nrows-1)), dpi=dpi)
    ax = plt.gca()
    

    if mins_maxs is None:
        val_min, val_max = 1e-10, 1e-8
        val = df_for_heatmap[[c for c  in df_for_heatmap.columns if c != 'Total']]
        val_min, val_max = 1e-16 + val.min().min(), val.max().max()

        if totals: 
            tot_min, tot_max = 1e-7, .001
            tot = df_for_heatmap.T.Total
            tot_min, tot_max = 1e-16 + tot.min(), tot.max()
            
    
        # print(f'{val_min=}, {val_max=}, {tot_min=}, {tot_max=}')
        # print(f'{val_min}, {val_max}, {tot_min}, {tot_max}')
    else: 
        if totals:
            val_min, val_max, tot_min, tot_max = mins_maxs
        else:
            val_min, val_max = mins_maxs
        
    val_norm = mpl.colors.Normalize(vmin=val_min, vmax=val_max)
    
    col_chr_color  = df_for_heatmap.columns.to_series().apply(lambda x: x.get_chrom()).map(chroms_colors)
    
    trues = np.ones(df_for_heatmap.shape, dtype=bool)
    mask = ~trues
    if totals:
        mask[3, :] = True
    h1 = sns.heatmap(df_for_heatmap, 
                cmap="RdYlGn_r",
                cbar = False, 
                norm = val_norm,
                mask=mask, ax=ax)

    if totals:
        tot_norm = mpl.colors.LogNorm(vmin=tot_min, vmax=tot_max)
        
        mask = trues
        mask[3, :] = False
        h2 = sns.heatmap(df_for_heatmap, 
                    cmap="YlGnBu",
                    cbar = False,
                    norm = tot_norm,
                    mask=mask, ax=ax)

    h = h2 if totals else h1
    
    transform = ax.transData

    if add_chroms:
        for i, color in enumerate(col_chr_color):
            h.add_patch(plt.Rectangle(xy=(i, nrows), height=.5, width=1, color=color, lw=0,
                                       transform=transform, clip_on=False))#ax.get_xaxis_transform()
        handles = []
        for chr, color in chroms_colors.items():
            handles.append(mpatches.Patch(color=color, label=chr))
            
        ax.legend(handles=handles, ncols=7, frameon=False,loc ='upper left', 
                  bbox_to_anchor=(ncols*.3, nrows+1), bbox_transform=transform)
    
    ax.tick_params(left=False, bottom=False)
    plt.yticks(rotation=0, ha='right');
    ax.set_xticks([])
    ax.set_xlabel('', labelpad=10)#Loop
    ax.xaxis.set_label_position('top') 
    


    l = .12
    b = -.2
    w = .08
    h = .1
    cax1 = fig.add_axes([l, b, w, h], transform=ax.transAxes)
    cb1 = mpl.colorbar.ColorbarBase(cax1, orientation='horizontal', cmap='RdYlGn_r', norm=val_norm, ticks=[val_min, (val_max-val_min)/2, val_max])
    if totals:
        cax2 = fig.add_axes([l + .08 + .02, b, w, h], transform=transform)
        cb2 = mpl.colorbar.ColorbarBase(cax2, orientation='horizontal', cmap='YlGnBu', norm=tot_norm)#, ticks=[1, 100, t_max])
    return ax