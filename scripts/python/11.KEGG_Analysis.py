import sys
import os
import warnings
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import re

warnings.filterwarnings("ignore")


if len(sys.argv) < 3:
    print("Error: Missing file arguments!")
    sys.exit(1)

blood_csv_file = sys.argv[1]
output_dir     = sys.argv[2]

def clean_term_name(term):
   
    return re.sub(r"\s*\(GO:\d+\)", "", term).strip()

def dotplot(df_sig, title, filename, top_n=15, figsize=(9, 7)):
    
    d = df_sig.head(top_n).copy()
    d["Term_clean"] = d["Term"].apply(clean_term_name)
    d["GeneCount"] = d["Overlap"].apply(lambda x: int(x.split("/")[0]))
    d["neg_log10_p"] = -np.log10(d["Adjusted P-value"])
    d = d.iloc[::-1].reset_index(drop=True)  # most significant at top

    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(
        d["neg_log10_p"], d["Term_clean"],
        s=d["GeneCount"] * 60,                 # dot size scaled by gene count
        c=d["neg_log10_p"],              # color by significance
        cmap="viridis_r",
        edgecolors="black", linewidths=0.5, alpha=0.85, zorder=3
    )
    ax.hlines(d["Term_clean"], xmin=0, xmax=d["neg_log10_p"], color="grey",
              linewidth=0.8, alpha=0.5, zorder=1)

    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.tick_params(axis='y', labelsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=1, alpha=0.6)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("-log10(adj. P)", fontsize=9)

    # legend for dot size (gene count)
    
    for count in sorted(d["GeneCount"].unique())[:3]:
        ax.scatter([], [], s=count * 60, c="grey", edgecolors="black",
                   alpha=0.7, label=f"{count} genes")
    ax.legend(scatterpoints=1, frameon=True, labelspacing=1.2,
             title="Gene Count", loc="upper left",
             bbox_to_anchor=(1.55, 0.30), fontsize=8, title_fontsize=9)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

def run_kegg_final():
   
    try:
        df = pd.read_csv(blood_csv_file)
    except FileNotFoundError:
        print(f"Error: '{blood_csv_file}' file not found!")
        return


    # EXCLUDE protein-coding genes only

    exclude_patterns = [
        r"^MIR\d",     # microRNAs
        r"^SNOR",      # snoRNAs
        r"^SAPCD1$",   # readthrough partner alone
        r"-SAPCD1$",   # fusion/readthrough transcript
        r"^AK\d{6}$",  # GenBank clone/EST IDs
    ]

    def is_excluded(gene):
        return any(re.search(pat, gene) for pat in exclude_patterns)

    # Process genes with exclusion
    gene_list = []
    for s in df['Gene_Symbol'].dropna():
        parts = str(s).replace(' /// ', ',').replace('/', ',').split(',')
        for g in parts:
            g = g.strip()
            if g and g != 'N/A' and not is_excluded(g):
                gene_list.append(g)

    unique_genes = list(set(gene_list))
    print(f"Total Unique Genes found: {len(unique_genes)}")

    # Run KEGG Analysis
    print("Running KEGG Analysis...")
    enr = gp.enrichr(gene_list=unique_genes,
                     gene_sets=['KEGG_2021_Human'],
                     organism='human',
                     outdir=None,
                     cutoff=0.05)

    full_results = enr.results

    # Save Results as CSV 
    
    csv_path = os.path.join(output_dir, "Full_KEGG_Analysis_62_Genes.csv")
    full_results.to_csv(csv_path, index=False)
    print(f" Analysis saved to '{csv_path}'.")

    # Filter significant results
    kegg_sig = full_results[full_results["Adjusted P-value"] < 0.05].reset_index(drop=True)
    print(f"KEGG significant pathways: {len(kegg_sig)}")

   
    if len(kegg_sig) > 0:
        dotplot(
            kegg_sig,
            "KEGG Pathway Enrichment\n(Blood-Detectable Panel, 57 Protein-Coding Genes)",
            os.path.join(output_dir, "KEGG_Pathways_Top10.png")
        )
        print(f" Top 10 chart saved to '{os.path.join(output_dir, 'KEGG_Pathways_Top10.png')}'.")
    else:
        print("No significant KEGG pathways found (adjusted p < 0.05)")

if __name__ == "__main__":
    run_kegg_final()
