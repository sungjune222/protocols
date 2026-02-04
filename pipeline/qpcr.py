import matplotlib

# Supports file saving only; GUI rendering is not available
matplotlib.use("Agg")

import numpy as np
import os
import pandas as pd
import pingouin as pg
import zipfile
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator
from pipeline.utils.env import find_env_dir
import seaborn as sns
import matplotlib.pyplot as plt

eds_name = "2026-01-27_234145_PUP_PIEZO.eds"
control_gene = "Actb"
control_group_keywords = ["control", "ctrl", "ctr", "wildtype", "wt", "normal", "nor", "healthy", "con"]

eds_qpcr_location = find_env_dir("EDS_DATA")
eds_path = os.path.join(eds_qpcr_location, eds_name)

z = zipfile.ZipFile(eds_path, 'r')
file = z.open("apldbio/sds/analysis_result.txt")
lines = [line.decode('utf-8').strip() for line in file.readlines()]

header = None
current_meta = None
data_list = []

for i, line in enumerate(lines):
    if not line: continue
    parts = line.split('\t')

    if parts[0] == "Well":
        header = parts
        continue
    
    if header and parts[0].isdigit():
        if len(header) != len(parts):
                raise ValueError(f"Header and data length mismatch at line {i}: {len(header)} vs {len(parts)}")
        current_meta = dict(zip(header, parts))
        continue
    
    if parts[0] == "Rn values" and current_meta is not None:
        rn_vals = [float(x) for x in parts[1:] if x.strip() != '']
        current_meta['Rn_Curve'] = np.array(rn_vals) #type: ignore
        data_list.append(current_meta)
        current_meta = None
if current_meta: data_list.append(current_meta)

qpcr_data = pd.DataFrame(data_list)
qpcr_data = qpcr_data[["Sample Name", "Detector", "Ct"]].copy()
qpcr_data["Ct"] = pd.to_numeric(qpcr_data["Ct"])
qpcr_data = qpcr_data[qpcr_data["Ct"] < 40]
qpcr_data.dropna(subset=["Ct"], inplace=True)

ct_values = qpcr_data.pivot(
    index='Sample Name',
    columns='Detector',
    values='Ct'
)

if control_gene not in ct_values.columns:
    raise ValueError(f"Control gene '{control_gene}' not found in the data.")

# ΔCt = Ct(target) - Ct(control)
d_ct = ct_values.subtract(ct_values[control_gene], axis=0) #type: ignore
assert isinstance(d_ct, pd.DataFrame)
d_ct = d_ct.drop(columns=[control_gene])

target_genes = d_ct.columns.tolist()
def clean_name(name):
    if "." in name:
        name = name.split(".")[-1]
    
    if "_" in name:
        name = name.rsplit("_", 1)[0]
        
    return name.strip()

d_ct["group"] = [clean_name(name) for name in d_ct.index]

results = []
plot_data_list = []
p_value_annotations = []

for gene in target_genes:
    qpcr_gene = d_ct[[gene, "group"]].dropna().rename(columns={gene: "dCt"})

    unique_groups = qpcr_gene["group"].unique().tolist()
    ctrl_name = next((g for g in unique_groups if g.lower() in control_group_keywords), None)
    if not ctrl_name:
        print(f"Skipping {gene}: Control group not found.")
        continue

    # Welch's ANOVA
    anova = pg.welch_anova(dv="dCt", between="group", data=qpcr_gene)
    anova_p_val = anova["p-unc"].values[0]

    other_groups = sorted([g for g in unique_groups if g != ctrl_name])
    control_vals = qpcr_gene[qpcr_gene["group"] == ctrl_name]["dCt"]
    
    ttest_p_values = []
    comparisons_ordered = []

    for group in other_groups:
        group_vals = qpcr_gene[qpcr_gene["group"] == group]["dCt"]

        # Control vs Others (Welch's t-test)
        ttest = pg.ttest(control_vals, group_vals, correction=True) #type: ignore
        ttest_p_val = ttest["p-val"].values[0]
        ttest_p_values.append(ttest_p_val)
        comparisons_ordered.append(group)
    
    # FDR Correction
    rejects, adj_p_values, _, _ = multipletests(ttest_p_values, alpha=0.05, method='fdr_bh')
    
    control_mean_d_ct = qpcr_gene[qpcr_gene["group"] == ctrl_name]["dCt"].mean()
    qpcr_gene["ddCt"] = qpcr_gene["dCt"] - control_mean_d_ct
    qpcr_gene["Fold Change"] = 2 ** (-qpcr_gene["ddCt"])
    qpcr_gene["Gene"] = gene

    for i, group in enumerate(comparisons_ordered):
        p_adj = adj_p_values[i]
        
        def get_star(p):
            if p > 0.05: return "ns"
            elif p <= 0.0001: return "****"
            elif p <= 0.001: return "***"
            elif p <= 0.01: return "**"
            else: return "*"
        
        mean_fc = qpcr_gene[qpcr_gene["group"] == group]["Fold Change"].mean()

        results.append({
            'Gene': gene,
            'ANOVA (Welch) P': anova_p_val,
            'Comparison': f"{group} vs Control",
            'Fold Change': mean_fc,
            'Welch T-test P': ttest_p_values[i],
            'FDR Adj P-value': p_adj,
            'Significant': get_star(p_adj)
        })

        p_value_annotations.append({
            "pairs": ((gene, ctrl_name), (gene, group)),
            "p_value": p_adj,
        })
    
    plot_data_list.append(qpcr_gene)

qpcr_results_location = find_env_dir("QPCR_RESULTS")
qpcr_results_location = os.path.join(qpcr_results_location, eds_name.replace(".eds", ""))
os.makedirs(
    qpcr_results_location,
    exist_ok=True,
)

results = pd.DataFrame(results)
results.to_csv(os.path.join(qpcr_results_location, "qpcr_analysis_results.csv"), index=False)

final_plot_df = pd.concat(plot_data_list)

all_groups = final_plot_df["group"].unique()
ctrl_name_final = next((g for g in all_groups if g.lower() in control_group_keywords), "Control")
other_groups_final = sorted([g for g in all_groups if g.lower() not in control_group_keywords])
plot_order = [ctrl_name_final] + other_groups_final

unique_genes = final_plot_df["Gene"].unique()
for gene in unique_genes:
    gene_data = final_plot_df[final_plot_df["Gene"] == gene]
    gene_annotations = [item for item in p_value_annotations if item["pairs"][0][0] == gene]
    gene_pairs = [item["pairs"] for item in gene_annotations]
    gene_p_values = [item["p_value"] for item in gene_annotations]

    plt.figure(figsize=(6, 6))
    ax = sns.barplot(x="Gene", y="Fold Change", hue="group", 
                        data=gene_data, 
                        hue_order=plot_order,
                        palette="tab10", 
                        errorbar='sd', 
                        capsize=0.1, 
                        err_kws={"color": "black", "linewidth": 1.5}, 
                        edgecolor="black", 
                        linewidth=1.0,
                        alpha=0.9)
    
    sns.stripplot(x="Gene", y="Fold Change", hue="group", 
                    data=gene_data, 
                    hue_order=plot_order,
                    dodge=True, 
                    jitter=0.2, 
                    color="black", 
                    marker="s",
                    size=6, 
                    alpha=0.7,
                    legend=False)

    annotator = Annotator(ax, gene_pairs, data=gene_data, x="Gene", y="Fold Change", hue="group", hue_order=plot_order)
    annotator.configure(test=None, text_format="star", loc="inside", verbose=False)
    annotator.set_pvalues_and_annotate(gene_p_values)

    plt.axhline(1, color="gray", linestyle="--", linewidth=1)
    plt.title("Relative Gene Expression (Fold Change)", fontsize=14, fontweight="bold")
    plt.ylabel("Fold Change ($2^{-\Delta\Delta C_t}$)") # pyright: ignore[reportInvalidStringEscapeSequence]

    handles, labels = ax.get_legend_handles_labels()
    l = len(plot_order)
    plt.legend(handles[:l], labels[:l], bbox_to_anchor=(1.05, 1), loc="upper left", title="Group")
    plt.tight_layout()
    plt.savefig(os.path.join(qpcr_results_location, f"{gene}_fold_change_plot.svg"), format="svg")
    plt.close()