# JUMPsem <img src="https://img.shields.io/badge/R-Package-blue.svg" alt="R badge" align="right"/>

## Tool to Calculate Enzyme Activity

**JUMPsem** is an R package developed by the **Wang Lab (UTHSC)** to infer enzyme activity from quantitative omics data (e.g., phosphoproteomics, ubiquitinomics, or acetylomics).  
It uses *structural equation modeling (SEM)* to estimate latent enzyme activity based on the coordinated regulation of its substrates, allowing researchers to translate complex site-level data into interpretable enzyme-level insights.

---

## 🚀 Key Features

- 🔬 **SEM-based enzyme activity inference** — integrates enzyme–substrate relationships into a latent variable model.
- 🧩 **Compatible with multiple data types** — supports phospho-, ubiquityl-, or acetyl-proteomics datasets.
- 📊 **Flexible model fitting** — supports covariance- or data-based SEM fitting with KMO filtering and residual correction.
- 🧠 **Single- or multi-enzyme modeling** — enables both individual enzyme activity inference (`singleEnzymeSEM`) and enzyme–enzyme interaction inference (`JUMPeei` integration).
- 🔁 **Bootstrap validation** — estimates the robustness and reproducibility of enzyme activity across replicates.
- 🎨 **Publication-ready outputs** — produces clean plots (activity curves, heatmaps, and ROC curves) and summarized reports.
- 🧱 **Integrates easily with tidyverse** — uses data frames as input/output for smooth integration with `dplyr`, `ggplot2`, and downstream analytics.

---

## 🧭 Installation

Install the latest released version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Wanglab-UTHSC/JUMPsem")
