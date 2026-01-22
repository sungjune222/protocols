# scRNA-seq Analysis Pipeline

This repository contains a pipeline for Single-cell RNA sequencing data analysis.

## 🚀 Workflow

The analysis process consists of two main steps:

### 1. Data Preparation
First, open and run the **`data_prep.ipynb`** notebook.
* This step loads the raw data files (e.g., `.rds`, `.csv`, `.mtx`, ...) and converts them into the **`.h5ad`** format for efficient processing.

### 2. Main Analysis
Once the data preparation is complete, run the **`main.py`** script to execute the full analysis pipeline.

```bash
pixi shell
python main.py