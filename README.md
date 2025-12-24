# Single-Cell RNA-seq Analysis Pipeline

## 📋 Overview

This is an **optimized, production-ready** single-cell RNA-seq analysis pipeline for comparing Wild-Type (WT) vs Knockout (KO) conditions across multiple tissues (Spleen and Bone Marrow).

**Key Features:**
- ✅ Complete end-to-end analysis (QC → DEG → Enrichment → Comparison)
- ✅ Centralized parameter management
- ✅ Publication-ready visualizations
- ✅ Fully reproducible with documented parameters

---

## 🚀 Quick Start

### Prerequisites

```r
# Install required packages
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "reshape2", 
                   "ggrepel", "RColorBrewer", "ggalluvial", "compositions"))

BiocManager::install(c("scDblFinder", "SingleR", "celldex", 
                       "SingleCellExperiment", "clusterProfiler", 
                       "org.Mm.eg.db", "fgsea", "msigdbr"))

# CACOA (if not installed)
devtools::install_github("kharchenkolab/cacoa")
```

### Running the Analysis

1. **Prepare your data:**
   ```
   DATA/
   ├── WT_SPsample_filtered_feature_bc_matrix.h5
   ├── KO_SP_sample_filtered_feature_bc_matrix.h5
   ├── WT_BM_sample_filtered_feature_bc_matrix.h5
   └── KO_BM_sample_filtered_feature_bc_matrix.h5
   ```

2. **Open RStudio:**
   - Open `scRNAseq_OPTIMIZED.Rmd`
   - Click "Knit"
   - Wait ~45 minutes
   - View results in `results/`

3. **Customize parameters (optional):**
   - Edit the `PARAMS` list in Section 2
   - Change QC thresholds, clustering resolution, etc.
   - Re-knit to see effects

---

## 📊 What Does It Do?

### For Each Tissue (Spleen & Bone Marrow):

1. **Quality Control**
   - Remove low-quality cells
   - Detect and remove doublets
   - Filter by MT%, ribosomal%, complexity

2. **Preprocessing**
   - Normalization
   - Cell cycle scoring & regression
   - Integration of WT and KO samples

3. **Dimensionality Reduction & Clustering**
   - PCA
   - UMAP
   - Leiden clustering

4. **Cell Type Annotation**
   - Automated annotation using SingleR (ImmGen reference)
   - Manual validation with canonical markers

5. **Differential Expression**
   - DEG analysis per cell type (KO vs WT)
   - Statistical testing (Wilcoxon rank-sum)
   - FDR correction

6. **Functional Enrichment**
   - GO enrichment (up/downregulated genes)
   - GSEA (pathway analysis)

7. **Compositional Analysis**
   - Cell type proportion changes
   - CACOA analysis (expression shifts)

8. **Cross-Tissue Comparison**
   - Correlation of DEGs across tissues
   - Scatter plots with statistics
   - Shared vs tissue-specific changes

---

## 📁 Output Structure

```
project/
├── scRNAseq.R                      # Main analysis pipeline (R)
├── scRNAseq.Rmd                    # Main analysis pipeline (RMarkdown)
├── scRNAseq.html                   # Rendered HTML report
├── results/                        # All analysis outputs (see below)
│   ├── spleen/                     # Spleen-specific results
│   ├── bonemarrow/                 # Bone marrow-specific results
│   ├── tissue_comparison/          # Cross-tissue comparison
│   └── cacoa/                      # CACOA analysis results
├── data/                           # Input data (not included in repo)
│   ├── spleen_wt.h5                # Spleen WT raw data
│   ├── spleen_ko.h5                # Spleen KO raw data
│   ├── bonemarrow_wt.h5            # Bone marrow WT raw data
│   └── bonemarrow_ko.h5            # Bone marrow KO raw data
└── README.md                       # This file
```

---

## ⚙️ Key Parameters

All parameters are defined in Section 2 (`PARAMS` list):

### QC Thresholds
```r
qc = list(
  min_features = 500,           # Min genes per cell
  max_mt_spleen = 10,           # Max mitochondrial % (spleen)
  max_mt_bonemarrow = 12,       # Max mitochondrial % (bone marrow)
  min_complexity = 0.8          # Min log10(genes)/log10(UMI)
)
```

### Clustering
```r
clustering = list(
  resolution = 0.5,             # Higher = more clusters
  dims = 1:30                   # PCs to use
)
```

### Differential Expression
```r
deg = list(
  test_use = "wilcox",          # Statistical test
  min_pct = 0.1,                # Min % cells expressing
  logfc_threshold = 0.25,       # Min log2FC
  padj_cutoff = 0.05            # Significance threshold
)
```

### CACOA
```r
cacoa = list(
  n_pseudo_reps = 3,            # Pseudo-replicates
  ref_level = "WT",             # Reference condition
  target_level = "KO"           # Test condition
)
```

**To change parameters:**
1. Edit the `PARAMS` list in Section 2
2. Save
3. Re-knit

---

## 🎨 Visualizations

### Plots Generated:

**Per Tissue:**
- UMAP (condition, clusters, cell types)
- QC violin plots (genes, UMI, MT%, ribosomal%)
- Doublet detection plots
- DEG volcano plots
- GO enrichment dot plots
- GSEA bar plots
- Compositional volcano plots
- Alluvial plots (cell type flow WT→KO)
- CACOA plots (abundance changes, expression shifts)

**Cross-Tissue:**
- Correlation scatter plots (per cell type)
- Combined panel plot (all cell types)
- R² summary bar plot
- Correlation heatmap

**All plots are:**
- Publication-ready (high resolution, 300 DPI)
- Saved as PDF and PNG
- Consistently styled
- Properly annotated

---

## 🔬 Scientific Workflow

```
Raw Data (h5 files)
    ↓
QC & Filtering (remove low-quality cells, doublets)
    ↓
Normalization & Integration (combine WT + KO)
    ↓
Dimensionality Reduction (PCA, UMAP)
    ↓
Clustering (Leiden algorithm)
    ↓
Cell Type Annotation (SingleR + ImmGen)
    ↓
Differential Expression (per cell type, KO vs WT)
    ↓
Functional Enrichment (GO, GSEA)
    ↓
Compositional Analysis (cell type proportion changes)
    ↓
Cross-Tissue Comparison (shared vs specific changes)
    ↓
Publication-Ready Results!
```

---

## 💡 Advanced Features

### 1. Master Wrapper Function

Instead of repeating 500 lines of code for each tissue, we use ONE function:

```r
results <- analyze_tissue_complete(
  wt_seurat = wt_data,
  ko_seurat = ko_data,
  tissue_name = "YourTissue",
  params = PARAMS
)

# Returns everything:
# - results$seurat (annotated Seurat object)
# - results$deg (differential expression)
# - results$go (GO enrichment)
# - results$gsea (GSEA results)
# - results$compositional (cell type changes)
# - results$cacoa (CACOA object)
```

**Benefits:**
- Add new tissue in 35 lines (vs 500)
- Guaranteed consistency
- Easy to maintain

### 2. Parameter Management

All settings in ONE place:

```r
# Try different resolutions
PARAMS$clustering$resolution <- 0.3  # Fewer clusters
# Re-knit to see effect

# Or run multiple resolutions
for (res in c(0.3, 0.5, 0.8)) {
  PARAMS$clustering$resolution <- res
  # Run analysis...
}
```

### 3. Reproducibility

Everything is saved:
```r
# Load exact parameters used
params <- readRDS("results/analysis_parameters.rds")

# Load complete results
spleen_results <- readRDS("results/spleen/spleen_complete_results.rds")

# Reload Seurat object
spl_integrated <- readRDS("results/spleen/spleen_integrated_annotated.rds")
```

---

## 🐛 Troubleshooting

### Common Issues:

**Error: "object 'deg_results_spl_df' not found"**
- Solution: Make sure Section 7 (Spleen Analysis) ran successfully before Section 10 (Tissue Comparison)

**Error: "cannot open file 'DATA/...'"**
- Solution: Check that your h5 files are in the DATA/ directory
- Or update the file paths in Section 3

**Memory Error**
- Solution: Reduce n_variable_features in PARAMS (2000 → 1000)
- Or increase system RAM
- Or reduce number of cells analyzed

**Different results than original**
- Check: Are PARAMS identical to your original settings?
- Verify: Same random seed? (default: 1234)
- Note: Small differences are expected due to random initialization

---

## 📖 Code Structure

### Main Sections:

1. **Packages** - Load all required libraries
2. **Global Parameters** - Define all analysis settings
3. **Data Import** - Read h5 files with error checking
4. **QC Functions** - Quality control utilities
5. **QC Metrics & Filtering** - Apply QC to all samples
6. **Analysis Functions** - All analysis tools + wrapper
7. **Cell Cycle Genes** - Prepare S/G2M genes
8. **Spleen Analysis** - Complete analysis in 35 lines!
9. **Bone Marrow Analysis** - Complete analysis in 35 lines!
10. **Cross-Tissue Comparison** - Correlation analysis
11. **Session Info** - R version, package versions

### Key Functions:

- `apply_qc_filtering()` - QC filtering with auto-detection
- `annotate_immune_cells()` - SingleR annotation
- `find_degs_per_celltype()` - DEG analysis
- `run_go_enrichment()` - GO enrichment
- `run_gsea_per_celltype()` - GSEA
- `run_cacoa_analysis()` - CACOA compositional
- `run_compositional_analysis()` - Fisher's exact test
- `analyze_tissue_complete()` - **MASTER WRAPPER** ⭐
- `plot_tissue_comparison_improved()` - Cross-tissue correlation

---

## 🤝 Contributing

To extend this pipeline:

### Add a New Tissue:

```r
# Just call the wrapper!
results_liver <- analyze_tissue_complete(
  wt_seurat = liver_WT_seurat,
  ko_seurat = liver_KO_seurat,
  tissue_name = "Liver",
  params = PARAMS,
  s_genes = s.genes,
  g2m_genes = g2m.genes
)
```

### Add a New Analysis:

1. Create function in Section 6
2. Call it in wrapper function
3. Add to exports
4. Done!

### Modify Parameters:

1. Edit `PARAMS` list in Section 2
2. Re-knit
3. Compare results

---

## 📞 Contact

**Author:** Sergi Roig-Soucase
**Email:** serroisoupv@gmail.com
**Date:** December 24, 2025

---

## 🙏 Acknowledgments

- Seurat team (Hao et al., 2021)
- SingleR developers (Aran et al., 2019)
- CACOA authors (Barkas et al., 2019)
- All package maintainers

---

## 📊 References

- Hao et al. (2021) Integrated analysis of multimodal single-cell data. Cell.
- Aran et al. (2019) Reference-based analysis of lung single-cell sequencing. Nat Immunol.
- Barkas et al. (2019) Joint analysis of heterogeneous single-cell RNA-seq datasets. Nat Methods.
- Squair et al., (2021) Confronting false discoveries in single-cell differential expression. Nat Commun
---

**Happy Analyzing! 🧬✨**
