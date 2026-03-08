# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research paper repository: "Ultra-long-range Polycomb-coupled interactions underlie subtype identity of human cortical neurons" (Pletenev, Vaulin et al., bioRxiv 2025). The code analyzes single-nucleus multiome (snm3C-seq) Hi-C and scRNA-seq data to study Polycomb-mediated chromatin loops in human neuronal subtypes.

## Conda Environments

Each notebook header specifies which environment to use. Restore an environment with:
```bash
conda env create -f conda_envs/<env>.yml
conda activate <env_name>
```

| File | Environment name | Primary use |
|------|-----------------|-------------|
| `conda_envs/cool.yml` | `cool` | Main Hi-C analysis (cooler, cooltools, coolpuppy, bioframe, pydeseq2, scanpy) |
| `conda_envs/coolpup.yml` | `coolpup` | Older coolpuppy-based analysis |
| `conda_envs/scrnaseq.yml` | `scrnaseq_env` | scRNA-seq analysis |
| `conda_envs/hictk.yml` | `hictk` | hictk-based Hi-C processing |
| `conda_envs/hic_nikita.yml` | `hic-env` | Hi-C preprocessing |

Run Jupyter notebooks with:
```bash
conda activate <env>
jupyter notebook
```

## Code Architecture

### Shared utilities
- `code/utils.py` — shared Python functions for Hi-C analysis imported by multiple notebooks:
  - `bedpe_to_bed()`, `regs_to_bins()` — genomic interval helpers
  - `reg_obs_exp()`, `obs_for_bin_pairs()` — fetch observed/expected Hi-C values from cooler files
  - `subt_snips()` — extract obs/exp Hi-C snippets per cell type using `cooltools.ObsExpSnipper`
  - `plot_hic_reg_zoom_out()`, `plot_av_dot_2d()`, `plot_gene_dots()`, `plot_loads()` — Hi-C visualization
  - `Pb_de` class — pseudobulk differential expression using pyDESEQ2 (PN vs IN cell types)

- `code/Loops_vs_Expression/scripts/my_loops_module.py` — data structures for loop analysis:
  - `Loop` class — hashable, orderable representation of genomic loop (bedpe format); always normalizes so anchor1 < anchor2
  - `Chromosome` class — orderable chromosome representation with chr1–chrX + trans ordering

- `code/Loops_vs_Expression/scripts/my_network_utils.py` — network/graph utilities for loop analysis

### Analysis pipeline (notebook sequence)

**Data preprocessing:**
1. `tian2023.make_cool_from_pairs.ipynb` — build .cool files from Tian et al. snm3C-seq pairs
2. `make_clr_heffel.ipynb` — build .cool files from Heffel et al. data

**Loop annotation:**
3. `Loop annotation.ipynb` — call Polycomb loops
4. `combine_dot_annotation.ipynb` — merge loop calls

**Adult brain analysis** (Tian et al. Hi-C + Siletti et al. scRNA-seq):
5. `Loops_vs_Expression/` notebooks 1–7 — call PcG dots, extract intensities per cell type, get expression, heatmaps, PCA
6. `tian2023_analysis.v2.ipynb` — main adult brain analysis

**Developing brain analysis** (Heffel et al. Hi-C + Herring et al. scRNA-seq):
7. `3D_genome_in_development.Heffel2024.ipynb`
8. `brain_scRNA-seq.herring2022.ipynb`
9. `Loops_vs_Expression/` notebooks 8–9

**ChIP-seq analysis:**
10. `cgi_motifs.ipynb`, `ctcf_and_polycomb_loops_compare.ipynb`, `interval_clusterisation_by_chipseq_signal.ipynb`, `polycomb_dots_amount_with_chipseq_signal.ipynb`

### Key data formats
- Hi-C data: `.cool`/`.mcool` (cooler format), accessed via `cooler.Cooler`
- Loops: bedpe format (chrom1, start1, end1, chrom2, start2, end2)
- Genome: hg38; chromosomes chr1–chr21, chrX (chrY/chrM excluded from most analyses)
- Key output: 263 Polycomb anchors in `data/polycomb_dot_anchors.5kb.16_06_25.csv`
- Cell types distinguished: EN (excitatory/pyramidal neurons, PN) vs IN (inhibitory interneurons)
