import matplotlib
matplotlib.use("Agg")

import os
import pandas as pd
import pingouin as pg
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator
from pipeline.utils.env import find_env_dir
import seaborn as sns
import matplotlib.pyplot as plt

data_dict = {
        "Con":     [26.808, 27.268, 27.319, 27.115, 27.710, 28.316, 30.394, 29.354],
        "Tnfa_24": [27.611, 27.543, 27.757, 27.633, 28.194, 28.113, 28.520, 28.654],
        "IL1b_24": [28.822, 28.819, 29.478, 29.293, 28.190, 28.149, 29.159, 28.887],
        "IL27_24": [28.678, 28.625, 29.033, 28.622, 28.531, 29.468, 29.603, 28.410],
        "Tgfb_24": [29.424, 29.945, 28.342, 28.179, 27.603, 28.265, 29.210, 28.354],
        "Ifnr_24": [27.987, 27.576, 28.360, 28.144, 28.877, 27.747, 28.576, 28.427]
    }
gene_name = "Tafa5_raw" 
project = "Cytokine_24h"
ctrl_name = "Con"  

records = []
for group, values in data_dict.items():
    for val in values:
        if pd.notna(val): 
            records.append({
                "group": group,
                "Value": val,
                "Gene": gene_name
            })

df = pd.DataFrame(records)

header = None
current_meta = None
data_list = []

unique_groups = df["group"].unique().tolist()
if ctrl_name not in unique_groups:
    raise ValueError(f"Control group '{ctrl_name}' not found in the data.")

# Welch's ANOVA
anova = pg.welch_anova(dv="Value", between="group", data=df)
anova_p_val = anova["p-unc"].values[0]

# T-test
other_groups = [g for g in unique_groups if g != ctrl_name]
control_vals = df[df["group"] == ctrl_name]["Value"]

ttest_p_values = []
comparisons_ordered = []

for group in other_groups:
    group_vals = df[df["group"] == group]["Value"]
    
    # Control vs Others (Welch's t-test)
    ttest = pg.ttest(control_vals, group_vals, correction=True) #type: ignore
    ttest_p_values.append(ttest["p-val"].values[0])
    comparisons_ordered.append(group)

# FDR correction (Benjamini-Hochberg)
_, adj_p_values, _, _ = multipletests(ttest_p_values, alpha=0.05, method='fdr_bh')

control_mean_val = control_vals.mean()
df["ddCt"] = df["Value"] - control_mean_val
df["Fold Change"] = 2 ** (-df["ddCt"])

results = []
p_value_annotations = []

for i, group in enumerate(comparisons_ordered):
    p_adj = adj_p_values[i]
    
    def get_star(p):
        if p > 0.05: return "ns"
        elif p <= 0.0001: return "****"
        elif p <= 0.001: return "***"
        elif p <= 0.01: return "**"
        else: return "*"
    
    mean_fc = df[df["group"] == group]["Fold Change"].mean()

    results.append({
        'Gene': gene_name,
        'ANOVA (Welch) P': anova_p_val,
        'Comparison': f"{group} vs {ctrl_name}",
        'Fold Change': mean_fc,
        'Welch T-test P': ttest_p_values[i],
        'FDR Adj P-value': p_adj,
        'Significant': get_star(p_adj)
    })

    p_value_annotations.append({
        "pairs": ((gene_name, ctrl_name), (gene_name, group)),
        "p_value": p_adj,
    })

qpcr_results_location = find_env_dir("QPCR_RESULTS")
qpcr_results_location = os.path.join(qpcr_results_location, project)
os.makedirs(
    qpcr_results_location,
    exist_ok=True,
)
pd.DataFrame(results).to_csv(os.path.join(qpcr_results_location, f"{gene_name}.csv"), index=False)

plot_order = [ctrl_name] + other_groups
gene_pairs = [item["pairs"] for item in p_value_annotations]
gene_p_values = [item["p_value"] for item in p_value_annotations]

plt.figure(figsize=(8, 6))

# Barplot
ax = sns.barplot(x="Gene", y="Fold Change", hue="group", 
                 data=df, 
                 hue_order=plot_order,
                 palette="tab10", 
                 errorbar='sd', 
                 capsize=0.1, 
                 err_kws={"color": "black", "linewidth": 1.5}, 
                 edgecolor="black", 
                 linewidth=1.0,
                 alpha=0.9)

# Stripplot
sns.stripplot(x="Gene", y="Fold Change", hue="group", 
              data=df, 
              hue_order=plot_order,
              dodge=True, 
              jitter=0.2, 
              color="black", 
              marker="s",
              size=5, 
              alpha=0.7,
              legend=False)

if gene_pairs:
    annotator = Annotator(ax, gene_pairs, data=df, x="Gene", y="Fold Change", hue="group", hue_order=plot_order)
    annotator.configure(test=None, text_format="star", loc="inside", verbose=False)
    annotator.set_pvalues_and_annotate(gene_p_values)

plt.axhline(1, color="gray", linestyle="--", linewidth=1)
plt.title(f"Relative Gene Expression: {gene_name}", fontsize=14, fontweight="bold")
plt.ylabel("Fold Change ($2^{-\Delta\Delta C_t}$)") # pyright: ignore[reportInvalidStringEscapeSequence]

handles, labels = ax.get_legend_handles_labels()
l = len(plot_order)
plt.legend(handles[:l], labels[:l], bbox_to_anchor=(1.05, 1), loc="upper left", title="Group")

plt.tight_layout()

output_path = os.path.join(qpcr_results_location, f"{gene_name}_fold_change_plot.svg")
plt.savefig(output_path, format="svg")
plt.close()