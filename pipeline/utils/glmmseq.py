import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
from typing import Union, Iterable
from anndata import AnnData
from rpy2.robjects import pandas2ri
from rpy2.robjects import conversion, default_converter
from pipeline.utils.pseudobulk import pseudobulk
from pipeline.config.constants import CPU_CORE_COUNT

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

def pseudobulk_glmmseq_comp(
    adata: AnnData, 
    condition_col: str, group_test: str, group_control: str, 
    random_effect_col: str, # column name for random effect (e.g., patient ID)
    covariates: Union[str, Iterable[str]] = [],
    filter_dict: dict = {},
    min_cells: int = 30
) -> pd.DataFrame:
    group_keys = ["sample", condition_col]
    covariates = [covariates] if isinstance(covariates, str) else list(covariates)
    
    cond_mask = np.array(adata.obs[condition_col].isin([group_test, group_control]), dtype=bool)
    for col, vals in filter_dict.items():
        cond_mask &= adata.obs[col].isin(vals)
        
    counts_df, pb_meta = pseudobulk(adata[cond_mask], group_keys=group_keys, min_cells=min_cells)
    
    min_count_per_sample = 2
    min_samples = int(np.ceil(counts_df.shape[0] * 0.10))
    keep_genes = (counts_df >= min_count_per_sample).sum(axis=0) >= min_samples
    pb_counts = counts_df.loc[:, keep_genes]

    keep_genes_var = pb_counts.var(axis=0) > 0
    pb_counts = pb_counts.loc[:, keep_genes_var]
    
    assert isinstance(adata.obs, pd.DataFrame)
    obs_meta = adata.obs.drop_duplicates(subset=group_keys).copy()
    obs_meta["pseudobulk_id"] = obs_meta[group_keys].astype(str).agg("__".join, axis=1)
    meta = obs_meta.set_index("pseudobulk_id").loc[pb_counts.index]
    meta["n_cells"] = pb_meta["n_cells"]
    meta["library_size"] = pb_meta["library_size"]
    
    design_cols = list(set(covariates + [condition_col, random_effect_col]))
    ok = meta[design_cols].notna().all(axis=1)
    meta = meta.loc[ok]

    for col in design_cols:
        if meta[col].dtype == 'object' or meta[col].dtype.name == 'category':
            meta[col] = meta[col].astype(str).str.replace(r'[^a-zA-Z0-9_.]', '_', regex=True).astype('category')
    pb_counts = pb_counts.loc[meta.index]
 
    formula_str_fixed = "~ " + " + ".join(covariates + [condition_col]) if covariates else f"~ {condition_col}"
    formula_str_glmm = formula_str_fixed + f" + (1 | {random_effect_col})"
    formula_str_reduced = "~ " + " + ".join(covariates) if covariates else "~ 1"
    formula_str_reduced = formula_str_reduced + f" + (1 | {random_effect_col})"

    # In R, the count matrix should be in the format [genes x samples], so we apply Transpose (.T)
    count_data_r = pb_counts.T 

    r_script = """
        suppressPackageStartupMessages({
            library(glmmSeq)
            library(DESeq2)
            library(ashr)
            library(parallel)
        })
        
        meta_df[[condition_col]] <- relevel(as.factor(meta_df[[condition_col]]), ref=group_control)
        meta_df[[random_effect_col]] <- as.factor(meta_df[[random_effect_col]])

        dds <- DESeqDataSetFromMatrix(as.matrix(countdata), meta_df, as.formula(formula_str_fixed))
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds, quiet=TRUE)

        robust_disp <- setNames(dispersions(dds), rownames(countdata))
        keep_disp <- is.finite(robust_disp) & !is.na(robust_disp) & (robust_disp > 0)
        if (!any(keep_disp)) stop("No genes with valid dispersion estimates.")
        countdata <- countdata[keep_disp, , drop=FALSE]
        robust_disp <- robust_disp[keep_disp]

        size_factor <- sizeFactors(dds)
        
        ctrl_samps <- rownames(meta_df)[meta_df[[condition_col]] == group_control]
        test_samps <- rownames(meta_df)[meta_df[[condition_col]] == group_test]
        
        cpm_mat <- sweep(as.matrix(countdata), 2, as.numeric(meta_df$library_size), "/") * 1e6
        mean_ctrl <- rowMeans(cpm_mat[, ctrl_samps, drop=FALSE])
        mean_test <- rowMeans(cpm_mat[, test_samps, drop=FALSE])

        genes <- rownames(countdata)
        num_chunks <- min(n_cores, length(genes))
        gene_chunks <- split(genes, cut(seq_along(genes), num_chunks, labels=FALSE))
        
        target_col <- paste0(condition_col, group_test)
        clust <- makeCluster(num_chunks, type="PSOCK")
        clusterExport(clust, c("countdata", "meta_df", "formula_str_glmm", "formula_str_reduced", "random_effect_col", 
                               "robust_disp", "target_col", "size_factor"), envir=environment())
        clusterEvalQ(clust, suppressPackageStartupMessages(library(glmmSeq)))
        
        res_list <- parLapply(clust, gene_chunks, function(g_chunk) {
            res <- data.frame(gene=g_chunk, coefs=NA_real_, se=NA_real_, pval=NA_real_, stringsAsFactors=FALSE)
            rownames(res) <- g_chunk
            
            fit <- tryCatch({
                glmmSeq(as.formula(formula_str_glmm), as.matrix(countdata[g_chunk, , drop=FALSE]),
                        metadata = meta_df, id=random_effect_col, dispersion=robust_disp[g_chunk], 
                        sizeFactors = size_factor, reduced = as.formula(formula_str_reduced), 
                        returnList = TRUE, cores=1, progress=FALSE)
            }, error = function(e) NULL)
            
            if (is.null(fit)) return(res)
            
            for (i in seq_along(fit)) {
                one_fit <- fit[[i]] 
                one_gene <- g_chunk[i]
                
                if (length(one_fit$coef) == 1 && is.na(one_fit$coef)) next
                if (length(one_fit$stdErr) == 1 && is.na(one_fit$stdErr)) next
                
                if (target_col %in% names(one_fit$coef)) { 
                    res[one_gene, "coefs"] <- as.numeric(one_fit$coef[target_col])
                    res[one_gene, "se"] <- as.numeric(one_fit$stdErr[target_col])
                } 
                
                if (!is.null(one_fit$chisq) && !is.null(one_fit$df) && 
                    !is.na(one_fit$chisq[1]) && !is.na(one_fit$df[1])) {
                    res[one_gene, "pval"] <- pchisq(one_fit$chisq[1], df=one_fit$df[1], lower.tail=FALSE) 
                } 
            } 
            return(res)
        })
        stopCluster(clust)
        
        comb <- do.call(rbind, res_list)
        padj_res <- as.numeric(p.adjust(comb$pval, method="BH"))
        
        lfc <- as.numeric(comb$coefs / log(2))
        lfc_se <- as.numeric(comb$se / log(2))
        
        ash_lfc <- rep(NA_real_, nrow(comb))
        ash_se <- rep(NA_real_, nrow(comb))
        valid <- !is.na(lfc) & !is.na(lfc_se)
        
        if (any(valid)) {
            ash_res <- ash(lfc[valid], lfc_se[valid], mixcompdist="normal")
            ash_lfc[valid] <- as.numeric(ash_res$result$PosteriorMean)
            ash_se[valid] <- as.numeric(ash_res$result$PosteriorSD)
        }
        
        res_df <- data.frame(
            gene = comb$gene,
            baseMean_control = as.numeric(mean_ctrl[comb$gene]),
            baseMean_test = as.numeric(mean_test[comb$gene]),
            log2FoldChange = lfc, 
            lfcSE = lfc_se,
            log2FoldChange_shrunk = ash_lfc,
            lfcSE_shrunk = ash_se, 
            pvalue = as.numeric(comb$pval),
            padj = padj_res,
            contrast = paste0(group_test, "_vs_", group_control),
            stringsAsFactors = FALSE
        )
    """

    with conversion.localconverter(default_converter + pandas2ri.converter):
        ro.globalenv['countdata'] = count_data_r
        ro.globalenv['meta_df'] = meta
        ro.globalenv['formula_str_fixed'] = formula_str_fixed
        ro.globalenv['formula_str_glmm'] = formula_str_glmm
        ro.globalenv['formula_str_reduced'] = formula_str_reduced
        ro.globalenv['condition_col'] = condition_col
        ro.globalenv['group_control'] = group_control
        ro.globalenv['group_test'] = group_test
        ro.globalenv['random_effect_col'] = random_effect_col
        ro.globalenv['n_cores'] = CPU_CORE_COUNT
        
        ro.r(r_script)
        result_df = ro.globalenv['res_df']

    return result_df