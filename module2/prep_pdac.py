import pandas as pd
from scipy import stats

def diff_expr(tumor_path, normal_path, out_path):
    tumor = pd.read_csv(tumor_path, sep="\t", index_col=0)
    normal = pd.read_csv(normal_path, sep="\t", index_col=0)

    # strip Ensembl version suffix (ENSG00000000003.15 -> ENSG00000000003)
    tumor.index = tumor.index.str.split(".").str[0]
    normal.index = normal.index.str.split(".").str[0]

    common = tumor.index.intersection(normal.index)
    tumor, normal = tumor.loc[common], normal.loc[common]

    log2fc, pval = {}, {}
    for g in common:
        t = tumor.loc[g].dropna()
        n = normal.loc[g].dropna()
        if len(t) < 3 or len(n) < 3:
            continue
        log2fc[g] = t.mean() - n.mean()
        pval[g] = stats.mannwhitneyu(t, n, alternative="two-sided").pvalue

    out = pd.DataFrame({
        "gene": list(log2fc),
        "log2fc": pd.Series(log2fc),
        "pval": pd.Series(pval)
    }).reset_index(drop=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"{out_path} -> {out.shape}")

print("Processing RNA (tumor vs normal)...")
diff_expr(
    "PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt",
    "data/pdac_transcriptomics.tsv"
)

print("Processing protein (tumor vs normal)...")
diff_expr(
    "PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt",
    "data/pdac_proteomics.tsv"
)

print("Processing mutation frequency (genomic layer)...")
mut = pd.read_csv("PDAC_somatic_mutation_gene_level_binary.txt", sep="\t", index_col=0)
mut.index = mut.index.str.split(".").str[0]
mut_freq = mut.mean(axis=1)
mut_freq.name = "mut_freq"
out = mut_freq.reset_index().rename(columns={"idx": "gene"})
out.to_csv("data/pdac_mutation.tsv", sep="\t", index=False)
print(f"data/pdac_mutation.tsv -> {out.shape}")

print("\nDone. All three gene-level files are in data/.")