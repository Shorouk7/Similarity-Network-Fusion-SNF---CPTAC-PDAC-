# Multi-Omics Tumor Subtyping using MOFA + Similarity Network Fusion (SNF)

### CPTAC Pancreatic Ductal Adenocarcinoma (CPTAC-PDAC)

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge\&logo=r\&logoColor=white)
![MOFA](https://img.shields.io/badge/MOFA-Multi--Omics-blueviolet?style=for-the-badge)
![SNF](https://img.shields.io/badge/SNF-Network%20Fusion-orange?style=for-the-badge)
![ClusterProfiler](https://img.shields.io/badge/clusterProfiler-Enrichment-green?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-Completed-brightgreen?style=for-the-badge)

---

## Author

**Assistant Lecturer:** Shorouf Aldiarabi
Faculty of Science, Port Said University, Egypt

---

## Overview

This project implements a **multi-omics integration framework** to identify **tumor molecular subtypes** in **CPTAC Pancreatic Ductal Adenocarcinoma (PDAC)**.

The pipeline combines:

* MOFA latent factor modeling
* Multi-omics feature extraction
* Similarity Network Fusion (SNF)
* Spectral clustering
* Functional enrichment analysis

We demonstrate that **latent tumor programs (Factor1)** capture **core oncogenic processes** and enable **stage-independent tumor stratification**.

---

## Pipeline Workflow

<p align="center">
  <img src="Figures/diagram.png" width="750">
</p>

---

## Objectives

* Extract latent tumor programs using MOFA
* Identify key omics features driving tumor biology
* Integrate multi-omics data using SNF
* Detect molecular tumor subtypes
* Interpret biological mechanisms per cluster

---

## Methods

### Multi-Omics Integration

8 omics layers were integrated:

* Proteomics
* Phosphoproteomics (Gene / Site / Multi-site)
* RNA-seq
* DNA methylation
* Glycoproteomics (Gene / Site)

---

### Feature Selection

* Top **300 positive Factor1 features** selected per omics layer

---

### Network Construction

* Z-score normalization
* Euclidean distance matrices
* Affinity matrices:

  * `K = 20`
  * `sigma = 0.5`

---

### Similarity Network Fusion (SNF)

* Iterative network fusion (`t = 20`)
* Integration of multi-layer similarity networks

---

### Clustering

* Spectral clustering (K = 2–6 tested)
* Final selection: **K = 3 tumor subtypes**

---

## Results Summary

### Identified Tumor Subtypes

| Cluster   | Biological Theme         | Interpretation                |
| --------- | ------------------------ | ----------------------------- |
| Cluster 1 | Immune / Inflammatory    | Immune-rich tumor             |
| Cluster 2 | Glycosylation remodeling | Invasive glyco-driven subtype |
| Cluster 3 | Signaling & Metabolism   | Highly active tumor program   |

---

## Key Findings

* **3 robust molecular subtypes (n = 117)**
* No association with tumor stage *(Chi-square p = 0.137)*
* Factor1 captures **core tumor program activity**

Tumor program activation:

* **Cluster 3 → Highest**
* **Cluster 1 → Intermediate**
* **Cluster 2 → Lowest**

---

## Biological Insights

* **Cluster 1:** Immune activation & stromal interaction
* **Cluster 2:** Glycosylation-driven invasion
* **Cluster 3:** Metabolic & signaling activation

---

## Figures

### Pipeline

![Pipeline](Figures/diagram.png)

### SNF Heatmap

![Heatmap](Figures/Heatmap-SNF.png)

### Tumor Program Score

![Score](Figures/Tumor-Prog-Score.png)

### GO Enrichment

![GO](Figures/GO-ENRICH.png)

---

## Project Structure

```
SNF-ANALYSIS/
│
├── Figures/
├── Script/
│   ├── Data/
│   ├── TOP-FEATURES/
│   └── SNF-ANALYSIS.R
```

---

## How to Run

```r
source("Script/SNF-ANALYSIS.R")
```

---

## Requirements

```r
MOFA2
SNFtool
tidyverse
ggplot2
pheatmap
clusterProfiler
enrichplot
pathview
```

---

## Notes

* Large raw datasets are excluded via `.gitignore`
* Pipeline is fully reproducible
* Designed for publication-grade analysis

---

## License

Academic and research use only.
