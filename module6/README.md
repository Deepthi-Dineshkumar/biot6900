# Module 6 · Computational Lab 6 — PRS Portability on Real Data (PDAC)

**BIOT 6900 · AI-Powered Biotech Workflows · Fall 2026 · Deepthi Dineshkumar**

This folder applies a published **pancreatic cancer polygenic score (PGS000725, 22 variants)** to **real 1000 Genomes phase 3 genotypes (2,504 people, 5 ancestry groups)** and asks who the score works for, who it fails, and why. It continues the PDAC theme from Modules 2 to 4.

The full write-up is in `Deepthi_Dineshkumar_Report.pdf`.

> ⚠️ **The phenotype is simulated.** 1000 Genomes has no disease labels. The PRS and its distribution across ancestries use real genotypes and real published weights; AUC and calibration use a simulated case label and are illustrative only.

---

## Folder contents

```
module6/
├── README.md                              this file
├── Deepthi_Dineshkumar_Report.pdf   full report with figures 
├── BIOT6900_Module6_Part3_PDAC.ipynb      graded notebook (run top to bottom)
├── build_inputs.sh                        downloads the score + 1000G genotypes, then runs make_data.py
├── make_data.py                           builds the four input CSVs and simulates the phenotype
├── data/
│   ├── genotypes.csv                      sample + 22 variant columns, dosage 0/1/2 of the ALT allele
│   ├── variants.csv                       variant, ref, alt (GRCh37)
│   ├── score.csv                          variant, effect_allele, other_allele, weight (PGS000725)
│   ├── meta.csv                           sample, ancestry, age, phenotype
│   └── SIMULATION_NOTES.txt               exact simulation settings and seed
└── figures/
    ├── prs_distribution_by_ancestry.png
    ├── calibration_by_ancestry.png
    ├── cutoff_by_ancestry.csv
    └── auc_by_ancestry.csv
```

## Data sources

| Input | Source |
|---|---|
| Score | [PGS000725](https://www.pgscatalog.org/score/PGS000725/) "PRS_Pancreas", Jia et al. 2020 (*JNCI Cancer Spectrum*); 22 genome-wide significant variants, beta weights, GWAS of 21,536 European-ancestry people. Harmonized GRCh37 scoring file from the PGS Catalog FTP. |
| Genotypes and REF/ALT | 1000 Genomes phase 3 (release 20130502, GRCh37), extracted at the 22 score positions with `bcftools` (biallelic SNPs only) |
| Ancestry | 1000 Genomes super-population panel (`integrated_call_samples_v3.20130502.ALL.panel`) |
| Age, phenotype | Simulated by `make_data.py` (seed 6900) |

**Simulation design.** Age ~ N(65, 8). Case status follows a logistic model with OR 1.5 per SD of genetic liability plus 0.06 log-odds per year of age. Liability = ρ × (PRS centred within group) + noise, with ρ = EUR 1.00, AMR 0.79, SAS 0.79, EAS 0.71, AFR 0.47 (the score's assumed tagging efficiency, based on Martin et al. 2019). Each group has 30% cases, so true risk is equal across ancestries by design. Exact settings are in `data/SIMULATION_NOTES.txt`.

## How to reproduce

### 1. Rebuild the data (optional; `data/` is already included)

Requires internet access, `bcftools` (with network support), and Python with numpy and pandas. On the Explorer HPC this was run inside a conda env:

```bash
conda create -n bcftools_env -c conda-forge -c bioconda python=3.11 numpy pandas scipy matplotlib scikit-learn ipykernel bcftools
conda activate bcftools_env
bash build_inputs.sh PGS000725
```

The script downloads the harmonized score, streams only the 22 positions from the remote 1000 Genomes VCFs (no full-chromosome downloads), and writes `data/`. Expected summary: 22 variants, 22 found, 2,504 people, about 30% cases per group. With the fixed seed, the CSVs are identical to the ones in this repo.

### 2. Run the notebook

Environment: Python 3.11 with numpy, pandas, scipy, matplotlib and scikit-learn (the course `biot6900` environment works).

1. Keep `data/` next to the notebook.
2. Confirm `USE_DEMO = False` in the setup cell.
3. Run **Kernel → Restart & Run All**. It runs in under a minute and writes everything in `figures/`.
4. Section 2 should print **source: REAL** with these checksums: genotypes `bbd9bb17`, variants `2a0cf215`, score `a8e3a9d6`, meta `b582b129`.

## Key results

| Ancestry | n | Mean PRS (z) | % above cohort top-10% cutoff | AUC (95% CI) | Predicted vs observed risk |
|---|---|---|---|---|---|
| EUR | 503 | −0.36 | 1.2 | 0.619 (0.566 to 0.672) | 0.286 vs 0.286 |
| AFR | 661 | +1.16 | 35.4 | 0.549 (0.503 to 0.597) | 0.510 vs 0.300 |
| EAS | 504 | −0.27 | 0.6 | 0.589 (0.533 to 0.638) | 0.297 vs 0.308 |
| SAS | 489 | −0.66 | 0.0 | 0.595 (0.546 to 0.650) | 0.249 vs 0.292 |
| AMR | 347 | −0.38 | 2.3 | 0.653 (0.593 to 0.710) | 0.284 vs 0.291 |

- **Harmonization:** 22/22 variants found and used; 12 flipped, 0 dropped, 0 not found; 3 palindromic SNPs flagged.
- **Incremental value:** cross-validated AUC 0.629 (age only) → 0.649 (age + PRS), +0.020.

## Interpretation notes

- **Portability.** The score was built from European GWAS hits that tag causal variants through European LD patterns. Shorter LD blocks and different allele frequencies in other ancestries weaken that tagging and shift the score's distribution.
- **One cutoff is not fair.** With equal true risk in every group, a single top-10% cutoff flagged 35.4% of AFR and 0% of SAS. This is entirely an allele-frequency effect, and applied clinically it would over-refer one group and exclude others.
- **AUC and calibration disagree.** AFR's AUC is only moderately lower, but its absolute risk is over-predicted by 21 points. Both metrics need to be reported per ancestry.
- **AMR's higher AUC is noise.** AMR is the smallest group, and its confidence interval overlaps EUR's.
- **Real vs built-in.** The distribution shift, cutoff fractions and AFR calibration failure come from real genotypes. The size of the AUC drop depends partly on the assumed ρ.
- **Clinical use.** The score adds little beyond age. At most, it could help prioritize surveillance within already high-risk groups in European-ancestry patients. It should not be used for population screening or in other ancestries without real-data validation.

## Troubleshooting notes

- **`conda install` killed on a login node:** the solver ran out of memory. Install on a compute node or via `sbatch` instead, and use `--solver=libmamba`.
- **Remote VCF fetch fails:** bcftools needs libcurl support. Alternatively, download the needed `ALL.chrN...vcf.gz` files plus their `.tbi` indexes and point `build_inputs.sh` at the local paths.

## References

- Jia G et al. *JNCI Cancer Spectrum* 2020;4(3):pkaa021 (PGS000725)
- Klein AP et al. *Nature Communications* 2018;9:556 (underlying GWAS)
- Martin AR et al. *Nature Genetics* 2019;51:584–591 (PRS portability and disparities)
- 1000 Genomes Project Consortium. *Nature* 2015;526:68–74
- Lambert SA et al. *Nature Genetics* 2021;53:420–425 (PGS Catalog)
