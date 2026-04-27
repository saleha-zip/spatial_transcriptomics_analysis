# Spatial Transcriptomics Analysis — 10x Genomics Visium & Xenium


## Project Overview

This repository documents a hands-on exploration of spatial transcriptomics using 10x Genomics platforms. Spatial transcriptomics extends conventional RNA sequencing by preserving the physical location of each measurement within the tissue, allowing gene expression to be interpreted in the context of tissue architecture and cellular organization. Four analysis notebooks are included, each following an official tutorial and extended with biological interpretation and additional visualizations.

---

## Platforms

### Visium
Visium is a sequencing-based spatial transcriptomics platform. A tissue section is placed on a slide containing thousands of spatially barcoded spots (55 µm in diameter), each capturing RNA from the tissue above it. The resulting data pairs a gene expression matrix with the tissue image — either H&E stained or fluorescence-labeled — enabling gene expression to be mapped back onto tissue morphology. Because each spot may contain multiple cells, Visium operates at near-single-cell rather than true single-cell resolution.

### Xenium
Xenium is an imaging-based in situ platform that achieves true single-cell resolution. Rather than extracting RNA for sequencing, Xenium detects transcripts directly within the tissue using fluorescent probes, assigning each detected transcript to an individually segmented cell. This allows precise spatial mapping of gene expression at the level of individual cells, at the cost of a smaller, targeted gene panel compared to sequencing-based approaches.

---

## Tools

### Scanpy
[Scanpy](https://scanpy.readthedocs.io/) is a Python library for single-cell gene expression analysis built around the `AnnData` data structure. It provides functions for quality control, normalization, dimensionality reduction (PCA, UMAP), and clustering. In this project it is used for the foundational Visium analysis pipeline in Notebook 1.

### Squidpy
[Squidpy](https://squidpy.readthedocs.io/) extends Scanpy with tools designed specifically for spatial omics data. It builds spatial neighborhood graphs from spot or cell coordinates and provides methods for spatial statistics including neighborhood enrichment, co-occurrence analysis, and image feature extraction. Squidpy is used across Notebooks 2, 3, and 4.

---

## Repository Structure

```
spatial_transcriptomics_analysis/
│
├── 01_scanpy_basic/          # Visium — basic QC, clustering, spatial visualization
├── 02_squidpy_visium_fluo/   # Visium — fluorescence image segmentation and feature extraction
├── 03_squidpy_visium_hne/    # Visium — spatial graph statistics on H&E data
└── 04_squidpy_xenium/        # Xenium — single-cell spatial analysis
```

Each subdirectory contains the analysis notebook, a detailed README walkthrough, and saved output figures.

---

## Notebooks

| # | Notebook | Dataset | Key Methods |
|---|---|---|---|
| 1 | Basic Scanpy Spatial | Human Lymph Node (Visium) | QC, normalization, Leiden clustering, spatial visualization, marker genes |
| 2 | Visium Fluorescence | Mouse Brain (Visium + fluorescence) | Nucleus segmentation, image feature extraction |
| 3 | Visium H&E | Mouse Brain (Visium + H&E) | Neighborhood enrichment, co-occurrence analysis |
| 4 | Xenium | Human Lung Cancer (Xenium) | Single-cell spatial clustering, centrality scores, Moran's I |

---

## Dependencies

```bash
# Notebooks 1–3
pip install scanpy squidpy igraph leidenalg scikit-image dask

# Notebook 4
pip install spatialdata spatialdata-io spatialdata-plot
```

All notebooks were developed and tested on Google Colab.

---

## References

- Piñeiro AJ, Houser AE, Ji AL. Research Techniques Made Simple: Spatial Transcriptomics. J Invest Dermatol. 2022 Apr;142(4):993-1001.e1. doi: 10.1016/j.jid.2021.12.014. PMID: 35331388; PMCID: PMC8969263.
- Wolf et al. (2018) SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. https://doi.org/10.1186/s13059-017-1382-0
- Palla et al. (2022) Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*. https://doi.org/10.1038/s41592-021-01358-2
- [Scanpy spatial tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)
- [Squidpy Visium fluorescence tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html)
- [Squidpy Visium H&E tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)
- [Squidpy Xenium tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)


