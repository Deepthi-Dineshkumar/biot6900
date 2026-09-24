# Module 2 Multi-Omics Target Discovery

**Name:** Deepthi Dineshkumar
**Disease:** Pancreatic Ductal Adenocarcinoma (PDAC)

## Datasets

All data was sourced from LinkedOmicsKB (CPTAC Pan-Cancer Data Freeze v1.2), 
https://kb.linkedomics.org/download. Everything here is open access, so no data use agreement was needed.

| Layer | Dataset | Source link | Access |
|---|---|---|---|
| Transcriptomics (RNA), tumor | RNAseq gene level, Tumor | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | Open |
| Transcriptomics (RNA), normal | RNAseq gene level, Normal | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | Open |
| Proteomics, tumor | Proteomics gene level abundance, Tumor | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | Open |
| Proteomics, normal | Proteomics gene level abundance, Normal | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | Open |
| Genomics | Somatic mutation, gene level (binary) | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_somatic_mutation_gene_level_binary.txt | Open |

The cohort is CPTAC PDAC, with around 140 tumor samples that have matched RNA, protein, and mutation calls, compared against separate normal cohorts for each layer. Integration happens at the gene level since the tumor vs normal comparison is unmatched, rather than using sample level correlation.

## Data processing
The raw matrices were turned into three gene level tables using `prep_pdac.py`:
- For RNA and protein, I computed the per gene tumor vs normal differential effect using a Mann Whitney U test. This produced `data/pdac_transcriptomics.tsv` and `data/pdac_proteomics.tsv`, each with the columns `gene`, `log2fc`, and `pval`.
- For genomics, I computed per gene mutation frequency, meaning the fraction of tumors that had a mutation in that gene. This produced `data/pdac_mutation.tsv` with the columns `gene` and `mut_freq`, which is used as the genomic evidence layer in place of a GWAS association score.
- The Ensembl gene IDs had their version suffix stripped before joining across layers, and the final top hit gene symbols were resolved using the `mygene` package.

## What's in this repo
- `BIOT6900_Module2_Starter.ipynb`, the notebook with Part 3 adapted for PDAC
- `prep_pdac.py`, the script that generates the three `data/pdac_*.tsv` files from the raw downloads
- `targets_pdac.csv`, the full ranked target table
- `report.pdf`, the written report

## Known limitations
The mutation file only covers around 3,362 genes, which bottlenecks the three way join down to 1,880 genes. RNF43, a known PDAC driver, was dropped from the join entirely, since it was present in the RNA and mutation data but missing from the proteomics file.
