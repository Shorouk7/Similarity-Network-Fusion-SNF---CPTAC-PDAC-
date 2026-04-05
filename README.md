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

This project implements a **multi-omics integration pipeline** to identify tumor molecular subtypes in CPTAC PDAC.

The pipeline integrates:

* MOFA latent factor analysis
* Multi-omics feature selection
* Similarity Network Fusion (SNF)
* Spectral clustering
* Functional enrichment analysis

---

## Pipeline Workflow

<p align="center">
  <img src="Figures/diagram.png" width="750">
</p>

---

## Objectives

* Extract tumor programs using MOFA
* Identify key multi-omics features
* Integrate omics layers using SNF
* Detect molecular tumor subtypes
* Characterize biological mechanisms

---

## Results Summary

### Tumor Subtypes

| Cluster   | Biological Theme       | Interpretation       |
| --------- | ---------------------- | -------------------- |
| Cluster 1 | Immune / Inflammatory  | Immune-rich subtype  |
| Cluster 2 | Glycosylation          | Invasive subtype     |
| Cluster 3 | Signaling & Metabolism | Active tumor program |

---

## Figures

---

### Cluster Composition

![Cluster](Figures/Cluster-composition.png)

---

### Tumor Stage Distribution

![Stage](Figures/Tumor-Stage-SNF.png)

---

### Tumor Stage vs Clusters (Alternative View)

![Stage2](Figures/TumorStages-Clusters.png)

---

### Tumor Program Score

![Score](Figures/Tumor-Prog-Score.png)

---

### Factor1 vs Tumor Stage

![Factor](Figures/Factor1-TumorStage.png)

---

### Omics Layer Activity

![Layer](Figures/Omic-Layer-Activity.png)

---

### SNF Heatmap

![Heatmap](Figures/Heatmap-SNF.png)

---

### Feature Heatmap (Key Drivers)

![Heatmap2](Figures/Heatmap.png)

---

### Top Features Driving Clusters

![TopFeatures](Figures/TopFeature-SNF.png)

---

### GO Enrichment (Dotplot)

![GO](Figures/GO-ENRICH.png)

---

### GO Enrichment (Barplot)

![GOBar](Figures/GO-ENRICH-BAR.png)

---

## Project Structure

```
SNF-ANALYSIS/
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
clusterProfiler
enrichplot
pathview
pheatmap
```

---

## License

Academic use only.

