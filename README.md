# Protocols

This repository contains pipelines for biological data analysis.

## 🚀 Workflow

### 0. Environment Setup
```bash
# 1) clone repository
git clone https://github.com/sungjune222/protocols.git
cd protocols

# 2) sync pixi dependencies
pixi env sync

# 3) activate environment
pixi shell

# 4) install R dependencies (You must use --vanilla to avoid errors)
Rscript --vanilla install.R
exit

# 5) update R dependencies
pixi shell
```
## scRNA-seq Analysis

### 1. Data Preparation
First, open and run the **`notebooks/data_prep.ipynb`** notebook.
* This step loads the raw data files (e.g., `.rds`, `.csv`, `.mtx`, ...) and converts them into the **`.h5ad`** format for efficient processing.

### 2. Main Analysis
Once the data preparation is complete, run the **`pipeline/main.py`** script to execute the full analysis pipeline.

### 3. Cell Type Annotation
Open and run the **`notebooks/celltype_annotation.ipynb`** to perform cell type annotation based on marker genes.
