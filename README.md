# scRNA-seq Analysis Pipeline

This repository contains a pipeline for Single-cell RNA sequencing data analysis.

## 🚀 Workflow

The analysis process consists of two main steps:

### 1. Data Preparation
First, open and run the **`notebooks/data_prep.ipynb`** notebook.
* This step loads the raw data files (e.g., `.rds`, `.csv`, `.mtx`, ...) and converts them into the **`.h5ad`** format for efficient processing.

### 2. Main Analysis
Once the data preparation is complete, run the **`pipeline/main.py`** script to execute the full analysis pipeline.

```bash
pixi shell
python main.py
```

### 3. Cell Type Annotation
Open and run the **`notebooks/celltype_annotation.ipynb`** to perform cell type annotation based on marker genes.

