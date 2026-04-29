# MTM1R2therDNM2_TA_SPS26

# Proteomics Data Analysis Pipeline (Perseus + R)

##  Overview

This repository contains a complete workflow for the analysis of proteomics data obtained from mass spectrometry (MS). The pipeline integrates:

* **Raw protein intensity data processing** using *Perseus*
* **Statistical analysis and visualization** using *R*

The goal of this project is to identify significantly regulated proteins across different experimental conditions and visualize them using volcano plots.

---

##  Project Structure

```
├── data/
│   ├── 260218_disease_control.csv
│   ├── 260218_treated_disease.csv
│   └── 260218_treated_control.csv
│
├── scripts/
│   └── 260218_analysis.r
│
├── perseus/
│   └── Persus_session.sps
│   └──persus_session.csv #input raw file for the persus session
│
└── README.md
```

---

## 🧪 Workflow Summary

### 1. Raw Data Processing (Perseus)

* Input: Protein intensity data obtained from mass spectrometry (MS)
* Steps performed in Perseus:

  * Data filtering and cleaning
  * Log2 transformation 
  * Normalization 
  * Imputing missing values
  * Statistical testing (e.g., t-test)
* Output:

  * Processed datasets with:

    * log2 fold change
    * p-values
    * -log10(p-value)
    * analysed the level of rescue by the treatment

The full session is provided in:

```
perseus/Persus_session.sps
```

---

### 2. Downstream Analysis (R)

The R script performs:

* Data import
* Threshold-based classification of proteins
* Volcano plot generation
* Highlighting top significantly regulated proteins
* Analysed the rescue by the treatment using the formula followed by further downstream analysis of the rescued proteins
* Go enrichment (clusteprofiler, Reactome, GSEA)

Script location:

```
scripts/260218_analysis.r
```

---

## 📊 Key Features of the R Analysis

### Volcano Plot Generation

* Thresholds used:

  * **log2 Fold Change cutoff:** 0.5
  * **Significance cutoff (-log10 p-value):** 1.3
    *   breaks = c(-Inf, 0, 30, 80, 120, Inf),
  labels = c("Worsened", "Not rescued", "Partially rescued", "Rescued", "over-rescued")


* Protein classification:

  * **Upregulated**
  * **Downregulated**
  * **Worsened**
  * **Not rescued**
  * **Partially rescued**
  * **Rescued**
  * **over-rescued**

Author
Supriya Priyadarshani SWAIN
