# BIOT6900 Portfolio Deliverable 2: Genomics & Precision Medicine Case Study (PDAC)

**Name:** Deepthi Dineshkumar
**Disease:** Pancreatic ductal adenocarcinoma (PDAC)
**Full write-up:** `Deepthi_Dineshkumar_PD2_Report.pdf` (8 pages)

## Overview

This folder packages Labs 5 to 8 into one precision-medicine argument and connects it back to the
Part I targets (XDH, MUC16, CEL, plus CFTR). **No primary analysis was re-run.** `pd2_figures.py`
reads the outputs already committed in the Module 5 to 8 folders and redraws them as the four
required figures; `pd2_build_report.py` assembles the PDF from those figures plus the tables.
The only additions are presentational: figures redrawn for legibility from the committed outputs.

## The argument in one paragraph

The case study assembles Labs 5 to 8 around the Part I shortlist. Variant interpretation in
CYP2C19 could not reach a patient-level call: population allele frequencies, not genotypes, and
pooled calibration dominated by DPYD. A European-derived PDAC PRS (PGS000725) discriminates
weakly (overall AUC 0.574) and transfers badly across ancestries (35.4% of AFR versus 0.0% of SAS
above one pooled top-decile cutoff). The knowledge graph recovers KRAS from Part I rank 69 to
rank 1 on topology and drug annotation alone, and ranks XDH second. **The measurement that would
decide whether that convergence describes a real target, which cells express it, has not been
made in PDAC.**

## Section-to-lab map

| Report section | Source | Key committed outputs used |
| --- | --- | --- |
| 1. Disease, target, question | Part I (Modules 2-4) | `Portfolio/prioritized_targets_pdac.csv` |
| 2. Variant interpretation and PGx | Lab 5 | Module 5 notebook and report (CYP2C19, pooled panel) |
| 3. Polygenic risk modelling | Lab 6 | Module 6 notebook outputs (PGS000725) |
| 4. Knowledge graph | Lab 7 | `module7/biomarker_ranking.csv`, `my_knowledge_graph.png` |
| 5. Single-cell context | Lab 8 | `module8/figures/`, `results_summary.json` |
| 6-7. Synthesis, limitations | all | - |

## Interpretation notes

### Section 2, variants (Tables 2 and 3)
- Reports the Lab 5 CYP2C19 interpretation as submitted: 80 variants (38 drug response, 33 VUS,
  9 B/LB, 0 P/LP), both anchors recovered (*2 AF 0.168, *17 AF 0.200).
- MAX_AF 0.25 retained both anchors and removed 2 of 80; a 0.001 rare-disease cutoff would have
  removed 19 variants, 14 of them drug response or P/LP.
- CADD selected for coverage (68 of 80 scored) over REVEL and AlphaMissense (54 each, missense
  only). Table 2 is the pooled five-gene CADD calibration from the Module 5 report; the 25-30
  band was excluded as unstable and the 30-50 band's infinite LR+ reflects no benign variants.
- Coverage gap: *7 (rs72558186, c.819+2T>A) received no CADD, REVEL or AlphaMissense score.
  Missing scores are not benign evidence.
- Table 3 is the rs41291556 (*8, W120R) discordance from the Module 5 report, with the same
  pattern at *22 (R186P) and *25 (F448L).
- No diplotype or metabolizer phenotype was determined; the eight VUS that shifted toward
  pathogenic are computational only and are not clinical reclassifications.

### Section 3, PRS (Figure 1 and Table 4)
- **Phenotypes were simulated**; this is stated in the section's first sentence. Panels (a) and (b)
  of Figure 1 use genotypes and published weights only and are outcome-independent; panels (c)
  and (d) depend on the simulation.
- Figure 1a draws Gaussian summaries of the observed per-group means and standard deviations
  rather than the raw histogram, so that all five groups and the cutoff stay legible in one panel.
  The underlying histogram is in the Module 6 report.

### Section 4, knowledge graph (Figure 2 and Table 5)
- **Correction to the Lab 7 report:** it states the driver ranks as 25th, 24th and 26th. The
  committed `biomarker_ranking.csv` gives **BRCA2 23rd (0.415), TP53 24th (0.404), SMAD4 25th
  (0.383)**. The report uses the CSV values.
- Known-gene recovery: KRAS rises from Module 2 rank 69 to graph rank 1, entirely from pathway
  membership and two drug edges, not from stronger omics evidence.
- Discordance examined: MUC16 is Part I rank 1 but graph rank 11, because it has no drug edge and
  p_ab is only 0.151. The demotion measures curation coverage, not biology.
- Dated annotation: the daraxonrasib approved/investigational correction moved KRAS from 2nd to
  1st. PARP1 (16th) and TYMS, RRM1, TOP1, EGFR (tied at 0.500) hold druggability 1.0 with no
  Module 2 score, so their rank rests entirely on drug status.
- Figure 2 is the committed `my_knowledge_graph.png`, cropped to its content bounding box and
  placed at full page width. The original full-resolution PNG is committed alongside it for
  label-level inspection.

### Section 5, single-cell (Figure 3)
- **Lab 8 was run on breast cancer with ERBB2 (HER2), not on PDAC.** The report states this at the
  top of the section, reports the result as it stands, and leaves the PDAC localisation question
  from Part I open.
- Figure 3b and 3c are redrawn from the committed Lab 8 numbers. 4b replaces the original violin
  plot, whose category labels overlapped and were illegible at print size. NK (n=50) and DC (n=36)
  reached 10-14% but are annotated as too small to interpret rather than plotted as bars.
- The measurement still required is cell-type-specific expression of XDH, MUC16, CEL and CFTR in
  PDAC, annotated with pancreatic markers (KRT19/EPCAM ductal and malignant, PRSS1/CPA1/CTRB1
  acinar, INS/GCG islet, COL1A1/ACTA2 fibroblast and CAF, PTPRC/CD68 immune, PECAM1 endothelial).

## Files in this folder

| File | Description |
| --- | --- |
| `dineshkumar_deepthi_portfolio2.pdf` | The 12-page case study submitted on Canvas |
| `pd2_figures.py` | Builds Figures 1-4 from the committed Lab 5-8 outputs |
| `pd2_build_report.py` | Assembles the PDF from those figures plus the tables |
| `figures/fig1_prs.png` | PGS000725: distribution, selection rate, calibration, AUC |
| `figures/fig2_graph.png` | PDAC knowledge graph, 80 nodes and 101 edges |
| `figures/fig3_singlecell.png` | ERBB2 single-cell panels (Lab 8, breast cancer) |
| `README.md` | This file |

**Inputs it reads** (paths relative to the repo root):

| Input | Produced by |
| --- | --- |
| `module7/biomarker_ranking.csv` | Module 7 notebook (committed) |
| `module7/my_knowledge_graph.png` | Module 7 notebook (committed) |
| `module8/figures/1_8_target_umap.png` | Module 8 notebook (committed) |
| `Portfolio/prioritized_targets_pdac.csv` | Portfolio Deliverable 1 |

PRS values in Figure 1 are transcribed from the committed Module 6 notebook output tables; the
numbers are listed at the top of `pd2_figures.py` with comments saying which table each came from.

## How to reproduce

1. Clone the repo and activate the course environment. The scripts need only `pandas`, `numpy`,
   `matplotlib`, `pillow` and `reportlab`:

   ```bash
   git clone https://github.com/Deepthi-Dineshkumar/biot6900.git
   cd biot6900
   conda activate biot6900
   ```

2. Rebuild the figures, then the PDF:

   ```bash
   python Portfolio2/pd2_figures.py
   python Portfolio2/pd2_build_report.py
   ```

If an input is missing the script stops and names the file it could not find.

## What did not work, or was not possible

- **No PDAC single-cell or spatial analysis was performed.** Lab 8 was run on breast cancer with
  ERBB2. The PDAC run that would answer the Part I localisation question was not completed before
  the deadline. Section 5 reports the breast result as a methods control and states the gap, and
  the cell-type column of the Table 4 evidence matrix is empty for every Part I target as a
  result. This is the single largest limitation of the deliverable.
- **No diplotypes or metaboliser phenotypes** could be assigned in Section 2, because the Lab 5
  pipeline held population allele frequencies rather than individual genotypes.
- **PRS phenotypes and age were simulated**, so all discrimination, calibration and incremental
  value figures are conditional on the simulation and are not evidence about real patients.
- **Three palindromic SNPs in PGS000725 remain strand-ambiguous**; allele matching alone cannot
  resolve them.
- **No variant-level interpretation of XDH, MUC16 or CEL** was performed; Lab 5 covered the
  pharmacogene panel only. This is listed as future work in Section 7.
- **The AlphaGenome regulatory analysis was skipped** in Lab 5, leaving the regulatory assessment
  of *7 and *2 incomplete.
- **MUC16 still has no structural assessment** (14,507 residues, beyond the 2,700-residue limit of
  the Module 4 workflow), carried forward from Part I as a data gap.
