# BIOT6900 Module 2: Multi-Omics Target Discovery in Pancreatic Ductal Adenocarcinoma (PDAC)

**Name:** Deepthi Dineshkumar
**Disease:** Pancreatic Ductal Adenocarcinoma (PDAC)
**Full write up:** `report.pdf`

## 1. Disease and data choice

### Why PDAC
PDAC was chosen because it is one of the poorest five year survival rates of any solid tumor. Its core driver genes
(KRAS, TP53, SMAD4, CDKN2A) are well characterized, but the disease still has very few
effective targeted therapies. That gap between known biology and available drugs makes it a
good disease for target discovery. The well known drivers also give a clear way to check
whether the pipeline recovers real biology.

CPTAC measures RNA, protein, and somatic mutations in the same tumors (about 140 patients).
It also provides normal tissue samples for RNA and protein. That made a real tumor vs normal
comparison possible for each gene, which gives a stronger differential signal than ranking
genes by raw abundance.

### Datasets
All data came from LinkedOmicsKB (CPTAC Pan-Cancer Data Freeze v1.2),
https://kb.linkedomics.org/download. Every file is **open access**. No data use agreement or
application was required.

| Layer | Dataset | Source link | Access |
|---|---|---|---|
| Transcriptomics (RNA), tumor | RNAseq gene level, Tumor | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | Open |
| Transcriptomics (RNA), normal | RNAseq gene level, Normal | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | Open |
| Proteomics, tumor | Proteomics gene level abundance, Tumor | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | Open |
| Proteomics, normal | Proteomics gene level abundance, Normal | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | Open |
| Genomics | Somatic mutation, gene level (binary) | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_somatic_mutation_gene_level_binary.txt | Open |

**Samples:** about 140 tumors (matched across RNA, protein, and mutation), 21 RNA normals,
and 75 protein normals. The normal groups are not matched to the tumor patients or to each
other, so the layers are integrated at the gene level rather than the sample level.

## 2. Loading and reducing each layer to gene level (`prep_pdac.py`)

| Layer | Genes in raw file | Gene level value computed | Output file |
|---|---|---|---|
| RNA | 60,579 | Tumor vs normal log2 fold change and p value (Mann Whitney U test) | `data/pdac_transcriptomics.tsv` (`gene`, `log2fc`, `pval`) |
| Protein | 11,057 | Tumor vs normal log2 fold change and p value (Mann Whitney U test) | `data/pdac_proteomics.tsv` (`gene`, `log2fc`, `pval`) |
| Mutation | 3,362 | Mutation frequency: fraction of tumors with a somatic mutation in the gene | `data/pdac_mutation.tsv` (`gene`, `mut_freq`) |

**Note on the genomic layer:** CPTAC provided tumor mutation calls, and not GWAS statistics. So I
used mutation frequency as the genomic evidence layer. It measures how often a gene is mutated
within tumors, not its association with inherited disease risk, and I did not treat it as GWAS data.

## 3. Harmonization, concordance, and multi evidence score (notebook Part 3)

1. **Harmonization:** Ensembl gene IDs had their version suffix removed (for example
   ENSG00000133703.12 became ENSG00000133703) so the same gene matches across all three files.
2. **Join:** The three tables were inner joined on gene ID. Gene symbols were then added with
   the `mygene` package, so every gene in the output carries both its Ensembl ID and its symbol.
   I joined on the stable Ensembl ID rather than the symbol, because symbols can be ambiguous
   or outdated across datasets.
3. **Surviving genes:** **1,880 genes** were present in all three layers. The mutation file
   was the limiting factor, since it only contains genes with at least one observed mutation.
4. **Concordance:** A gene is concordant when its RNA and protein log2 fold changes have the
   same sign. **1,190 of 1,880 genes (63%)** were concordant. This flag does not affect the score.
5. **Multi evidence score:** Each layer (absolute RNA log2FC, absolute protein log2FC, and
   mutation frequency) was converted to a 0 to 1 rank percentile. This puts values measured in
   different units on the same scale. The three percentiles were then combined with
   **equal weights (one third each)**. I chose equal weighting because the mutation layer
   covers far fewer genes, and weighting it more heavily would favor genes just because they
   had mutation data.

## 4. Ranked target table

The full ranked list is exported to **`targets_pdac.csv`** (1,880 genes, sorted by score from
highest to lowest).

| Column | Meaning |
|---|---|
| `gene` | Ensembl gene ID (version removed) |
| `gene_symbol` | HGNC gene symbol |
| `rna_lfc` | RNA log2 fold change, tumor vs normal |
| `prot_lfc` | Protein log2 fold change, tumor vs normal |
| `mut_freq` | Fraction of tumors with a somatic mutation |
| `concordant` | True if RNA and protein change in the same direction |
| `score` | Equal weighted multi evidence score (0 to 1) |

### Top 10 targets
| Rank | Gene | RNA log2FC | Protein log2FC | Mut. freq. | Concordant | Score |
|---|---|---|---|---|---|---|
| 1 | MUC16 | 3.696 | 1.117 | 0.057 | True | 0.979 |
| 2 | CFTR | -2.438 | -1.819 | 0.029 | True | 0.978 |
| 3 | COL11A1 | 2.417 | 1.402 | 0.029 | True | 0.973 |
| 4 | CEL | -6.348 | -2.550 | 0.019 | True | 0.962 |
| 5 | LRP1B | -1.561 | -1.060 | 0.086 | True | 0.957 |
| 6 | BRSK2 | -3.083 | -1.955 | 0.019 | True | 0.956 |
| 7 | ACADL | -3.025 | -1.983 | 0.019 | True | 0.956 |
| 8 | NOS1 | -2.205 | -0.978 | 0.029 | True | 0.955 |
| 9 | XDH | 2.788 | 0.836 | 0.029 | True | 0.952 |
| 10 | KRT19 | 2.688 | 1.433 | 0.019 | True | 0.950 |

## 5. Biological interpretation

### Checking top hits against known PDAC biology
- **MUC16 (CA-125)** is a known tumor marker and is strongly increased at both RNA and protein level.
- **KRT19** is a pancreatic ductal cell marker that is commonly expressed in PDAC.
- **COL11A1** is linked to the dense stromal and fibrotic environment typical of PDAC.
- **CEL, CFTR, and ACADL** are strongly decreased. They mostly reflect loss of normal
  pancreatic tissue as tumor replaces it, so they are better seen as markers than as drug targets.

**Non obvious hit:** CFTR ranks second because of strong downregulation. It relates to normal
duct cell function and has been reported as altered in pancreatic cancer, but it is not usually
considered a major PDAC target. I treat it as a gene for further investigation, which could be
checked in a curated resource such as Open Targets.

### Recovery of known PDAC drivers
None of the major drivers appeared in the top 15, so I checked their ranks in the full table:

| Rank | Gene | RNA log2FC | Protein log2FC | Mut. freq. | Note |
|---|---|---|---|---|---|
| 32 | CDKN2A | 1.206 | 0.601 | 0.229 | Concordant, signal in all three layers |
| 69 | KRAS | 0.685 | 0.414 | 0.962 | Mutated in about 96% of tumors, modest RNA and protein change |
| 381 | TP53 | -0.046 | 0.803 | 0.829 | Discordant (see Section 6) |
| 398 | ARID1A | -0.370 | -0.193 | 0.057 | Weak signal in all layers |
| 530 | SMAD4 | -0.727 | -0.017 | 0.190 | Weak protein signal |
| n/a | RNF43 | present | missing | present | Dropped by the join, absent from proteomics |

KRAS drives PDAC through activating mutations, not by being made in much larger amounts. Its
importance shows clearly in the mutation layer but is diluted by equal weighting. This is why
combining layers matters: no single layer captures every type of driver.

## 6. Discordance and limitations

### Discordant gene: TP53
All top 15 genes were concordant, so I looked further down the ranking. TP53 (rank 381) has
almost no RNA change (-0.046) but a clear protein increase (0.803), and it is mutated in 83% of
tumors. Many TP53 missense mutations produce a stable mutant protein that escapes normal MDM2
driven degradation and accumulates in the cell. The protein level rises without any change in
transcription. This RNA vs protein mismatch is a real biological signal, not noise.

### Limitations
- **Unmatched samples:** Tumors and normals come from different people, and the RNA and protein
  normal groups differ in size. Results describe group level differences per gene, not changes
  within the same patient.
- **Gene coverage:** Only 1,880 genes survived the three way join, and known drivers like RNF43 were lost.
- **Mutation frequency is not GWAS data:** It reflects how often a gene is mutated in tumors,
  not inherited risk.
- **Stromal content:** PDAC tumors contain a lot of stroma, so bulk RNA and protein signals are
  not pure tumor cell biology.
- **Equal weighting:** This is a simple choice. Weighting mutation more heavily would raise KRAS
  but would also favor genes just for having mutation data.

## 7. Repository contents

| File | Description |
|---|---|
| `BIOT6900_Module2_Starter.ipynb` | Main notebook. Part 3 is adapted for PDAC |
| `prep_pdac.py` | Converts the five raw downloads into the three gene level tables in `data/` |
| `targets_pdac.csv` | Exported ranked target table |
| `report.pdf` | Full written report |
| `README.md` | This file |

The raw downloads and the `data/` folder are not committed because they are large and fully
reproducible with the steps below.

## 8. How to reproduce

**1. Clone the repo**

    git clone https://github.com/Deepthi-Dineshkumar/biot6900.git
    cd biot6900/module2

**2. Set up the environment**

    conda create -n biot6900 python=3.11
    conda activate biot6900
    pip install pandas numpy scipy mygene jupyter

**3. Download the raw data.** Save all five files from the Datasets table into `module2/`, next
to `prep_pdac.py`, and keep the original file names. Right click each link and choose
"Save Link As", because opening the link may show the file as text in the browser instead.

**4. Build the gene level tables**

    python prep_pdac.py

This creates `data/pdac_transcriptomics.tsv`, `data/pdac_proteomics.tsv`, and `data/pdac_mutation.tsv`.

**5. Run the notebook.** Start Jupyter from inside `module2/` so the `data/` paths resolve:

    jupyter notebook

Open `BIOT6900_Module2_Starter.ipynb` and choose Kernel, then Restart & Run All. An internet
connection is needed for the `mygene` symbol lookup.

**Expected results:** 1,880 genes after the join, 1,190 concordant (63%), MUC16 ranked first
with a score of about 0.979, and `targets_pdac.csv` written to `module2/`.
