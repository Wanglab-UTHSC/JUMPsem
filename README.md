<p align="center">
  <img src="doc/LOGO.png" alt="JUMPsem logo" width="360">
</p>

Toolbox for **inference of enzyme activity** using structural equation modeling (SEM) on phospho-/ubiquityl-/acetyl-proteomics data.

> Simple to install. Minimal inputs. Reproducible outputs with clear model‑fit diagnostics.

---

## Key Features

* **Unified API** for kinase, ligase, and HAT activity inference
* **SEM under the hood** (via `lavaan`) with optional residual‑covariance refinement
* **Automatic pre‑processing**: log2 transforms, normalization, PCA‑based substrate selection (optional)
* **Model‑fit reports**: CFI/TLI/RMSEA/SRMR and modification indices
* **Reproducible outputs**: timestamped file names and structured result lists
* **Friendly examples** bundled with the package

---

# Quick Installation of JUMPsem
**First, install devtools (for installing GitHub packages) if it isn't already installed:**
``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

**Then, install all dependencies:**
``` r
packages <- c("lavaan", "dplyr", "tidyr", "devtools", "psych", "MASS", "Matrix", "tidyverse", “EFAtools”, "data.table")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
}

```

**Lastly, install JUMPsem:**
``` r
devtools::install_github("Wanglab-UTHSC/JUMPsem")
```

## Quick Start

The package ships with small example datasets so you can run a complete workflow in seconds.

```r
library(JUMPsem)

# example inputs bundled with the package
data(input_psp_example)        # phospho substrate matrix
data(motif_example)            # motif annotations (optional)
data(wholeProteome_example)    # whole-proteome table (optional)

# run
res <- JUMPsem(
  input  = input_psp_example,
  datatype = "psp",                         # "psp" | "ubi" | "ace"
  organism = "mouse",                       # sample organism
  enzyme.organism = c("human","mouse","rat"),
  input.log2.norm = TRUE,                   # log2 transform substrate matrix
  whole.log2.trans = TRUE,                  # log2 transform whole-proteome (if provided)
  motif = motif_example,                    # optional motif table
  whole.proteome = wholeProteome_example,   # optional whole-proteome
  output.folder = "JUMPsem_results"         # output folder
)

# explore results
names(res)
#> e.g. c("Activity", "Affinity", "Evaluations")

head(res$Activity)
res$Evaluations
```

What you should see:

* An activity table per enzyme
* An affinity table per enzyme with z‑scores / p‑values (if applicable)
* A model‑fit summary (CFI, TLI, RMSEA, SRMR)
* A folder ` e.g. JUMPsem_results/` containing tab‑delimited outputs like `Activity_<timestamp>.txt`

---

## Example Usage

### 1) Minimal run (phospho‑proteomics)

```r
res <- JUMPsem(
  input = input_psp_example,
  datatype = "psp",
  organism = "mouse",
  enzyme.organism = c("human","mouse","rat")
)
```

### 2) With motif and whole‑proteome controls

```r
res <- JUMPsem(
  input = input_psp_example,
  datatype = "psp",
  organism  = "mouse",
  enzyme.organism  = c("human","mouse","rat"),
  motif = motif_example,
  whole.proteome = wholeProteome_example,
  input.log2.norm = TRUE,
  whole.log2.trans = TRUE
)
```

### 4) Single‑enzyme SEM (advanced)

If you want to fit or inspect a single enzyme model explicitly (e.g., PRKCD):

```r

se_res <- JUMPsem(
  input = input_psp_example,
  datatype = "psp",
  organism  = "mouse",
  enzyme.organism  = c("human","mouse","rat"),
  enzyList = "PRKCD",
  motif = motif_example,
  whole.proteome = wholeProteome_example,
  input.log2.norm = TRUE,
  whole.log2.trans = TRUE
)

```

---

## Inputs

* **`input`**: data frame / matrix (rows = substrates/sites, columns = samples); numeric values.
* **`datatype`**: one of `"psp"` (phosphorylation), `"ubi"` (ubiquitination), `"ace"` (acetylation).
* **`organism`**: substrate species for the experiment (e.g., `"human"`, `"mouse"`, `"rat"`).
* **`enzyme.organism`**: vector of enzyme species to consider when mapping (e.g., `c("human", "mouse", "rat")`).
* **`database`** *(optional; default `NULL`)*: custom enzyme–substrate reference; if `NULL`, uses the internal database.
* **`kmo.off`** *(optional; default `0`)*: KMO cutoff in `[0, 1]` to filter weak variables before SEM.
* **`mdsite`** *(optional; default `TRUE`)*: map using precise modification sites (e.g., S/T/Y positions); if `FALSE`, use gene-level mapping.
* **`enzyList`** *(optional; default `NULL`)*: restrict calculation to these enzyme IDs; if `NULL`, compute all available enzymes.
* **`input.log2.norm`** *(optional; default `FALSE`)*: apply log2 transform and normalization to the PTM matrix.
* **`relative.norm.p`** *(optional; default `TRUE`)*: apply relative normalization to the PTM input.
* **`whole.log2.trans`** *(optional; default `FALSE`)*: log2 transform the whole-proteome (ignored if `whole.proteome` is `NULL`).
* **`whole.proteome`** *(optional; default `NULL`)*: whole-proteome matrix used for normalization/adjustment.
* **`relative.norm.w`** *(optional; default `TRUE`)*: apply relative normalization to the whole-proteome input.
* **`motif`** *(optional; default `NULL`)*: motif table with additional kinase–substrate pairs for filtering/enrichment.
* **`output.folder`** *(optional; default `"."`)*: directory to save outputs; files are timestamped to avoid overwrites.

---


## Outputs

**Returned object** (list; field names may evolve):

* `Activity`: main per‑enzyme activity estimates
* `Affinity`: main per‑enzyme affinity estimates with z‑scores / p‑values (if applicable)
* `Evaluations`: a model‑fit summary (CFI, TLI, RMSEA, SRMR)

**Files written** into `output.folder`:

* `Activity_<timestamp>.txt`
* `Affinity_<timestamp>.txt`
* `Eval_table_<timestamp>.txt`

---

## Tips & Notes

* **Normalization**: set `relative.norm.p = TRUE` and `input.log2.norm = TRUE`to perform log2 + relative normalization of the input matrix
* **PCA selection**: optional PCA can winnow substrates to the most informative set (controls available via function arguments)

---

## RShiny website

Explore an interactive preview:

* JUMP Shiny: [https://jumpshiny.genenetwork.org/](https://jumpshiny.genenetwork.org/)

> Use the same input conventions as the examples above.

---

## Reproducibility

* Output files are **timestamped** to avoid accidental overwrites
* Randomness is minimal; set a seed if you use resampling‑based options

```r
set.seed(1)
```

---

## How to Cite

If you use JUMPsem in your work, please cite:

> Kong, Dehui et al. “A computational tool to infer enzyme activity using post-translational modification profiling data.” Communications biology vol. 8,1 103. 21 Jan. 2025, doi:10.1038/s42003-025-07548-4
---

## Getting Help

* **Issues & bugs**: please open a GitHub Issue with a minimal reproducible example
* **Questions / feature requests**: open a Discussion or Issue

---

## Contributing

Pull requests are welcome! Please:

1. Fork the repo and create a feature branch
2. Add tests or an example when relevant
3. Run `devtools::check()` before submitting

---

## License

MIT License (see `LICENSE`).
