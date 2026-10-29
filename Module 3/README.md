# BIOT6900 Module 3: Target Triage with Protein Language Models

**Name:** Deepthi Dineshkumar
**Disease:** Pancreatic Ductal Adenocarcinoma (PDAC), CPTAC PDAC cohort via LinkedOmics (continued from Module 2)
**Full write up:** `Module 3 Assignment Report.pdf`

## Overview
In Module 2, PDAC genes were ranked using RNA expression, protein abundance, and somatic mutation
evidence. In this module sequence based tractability prior built was added from the ESM C protein language model embeddings. It asks whether each target resembles proteins that already have clinical drug precedence. Module 2 evidence was combined with the prior to re-rank the targets, and the re-ranked list is the hand-off for structure based pocket assessment in Week 4.

## 1. Targets loaded and mapped, leakage guard applied

- **Input:** the top 15 targets from `module2/targets_pdac.csv`, the Module 2 output, which is already committed to this repo.
- **ID mapping:** Module 2 file stores Ensembl gene IDs, but `proteome.tsv` uses gene symbols.
- Further an additional step was added that converts the IDs to symbols, then matched each target to `proteome.tsv` while keeping the original embedding row positions.
- **Unmapped genes:** **0 of 15.** All 15 targets mapped.
- **Data pack:** 20,190 proteome proteins and a 1,201 protein training panel. Each protein is
  represented by a 960 dimensional ESM C 300M embedding.
- **Leakage guard:** Was applied before training. It removed **0 panel proteins**,
  because none of the 15 targets appear in the training panel. The guard still runs in the
  notebook, and its (empty) removal list is printed and saved when the notebook is run.

## 2. Prior trained and validated with a family aware split

I trained the any-clinical prior on the 1,201 protein panel and evaluated it with five fold cross
validation, using both a random split and a family aware split.

| Split | AUROC | Average precision (AP) |
|---|---|---|
| Random split (StratifiedKFold) | 0.848 | 0.727 |
| **Family split (GroupKFold on `family_group`)** | **0.732** | **0.548** |
| Chance level | 0.500 | 0.334 (prevalence of positives) |

**The family split values are the ones to trust.** In a random split, close paralogs from the same
protein family land on both sides of the split. Their embeddings are nearly identical, so the model
gets credit for recognizing family membership rather than for learning general features of
tractability. Holding out whole families tests whether the prior works on families it has never
seen labeled. That matters for this list, because several of my targets come from families with
little or no drugging history, such as mucins, keratins, and EPS8 like proteins. The random split
would have hidden a drop of 0.116 in AUROC and 0.179 in AP.

**Validation scope:** I only cross validated the any-clinical prior. I did not validate the small
molecule and antibody priors separately.

## 3. Modality priors, nearest drug targets, and anchor variant scores

### Modality priors
Each target received a small molecule prior (`p_sm`) and an antibody prior (`p_ab`). The report's
Figure 1 plots all 15 targets on these two axes.

- **Small molecule favored:** CEL (p_sm 0.716, p_ab 0.009) and XDH (p_sm 0.503, p_ab 0.003). Both
  are enzymes with defined catalytic pockets.
- **Antibody favored:** MUC5B and NPTX2 (p_ab about 0.50). Both are extracellular proteins without
  a classic enzyme pocket.
- **Mixed:** LRP1B (p_sm 0.284, p_ab 0.472).

### Nearest drug targets
For every target, the nearest proteins with clinical precedence (by cosine similarity of the
embeddings) are listed in `targets_pdac_w3.csv`. Some examples:

| Target | Nearest drugged neighbors | What it suggests |
|---|---|---|
| CFTR | GUCY2C, SCN7A, CASR (cosine 0.95) | Placed among multi-pass membrane drug targets |
| CEL | GPA33, MMP9, TNKS2 | Neighbors with different functions, so similarity does not prove a shared binding mechanism |

### Anchor variant scores
| Gene | Variant | Score | Within protein percentile |
|---|---|---|---|
| APOE | R176C (ε2) | -8.559 | 11.8 |
| APOE | C130R (ε4) | +7.208 | 100.0 |
| TREM2 | R47H | -5.156 | 47.7 |
| TP53 | R175H | -5.974 | 19.3 |

Positions use UniProt numbering, which includes the 18 residue signal peptide of APOE.

## 4. Re-ranked table exported

The re-ranked table is exported to **`targets_pdac_w3.csv`**. It uses a weight of **w = 0.6** for
the Module 2 evidence and 0.4 for the tractability prior.

| Rank | Gene | Module 2 rank | p_sm | p_ab | Final score |
|---|---|---|---|---|---|
| 1 | CEL | 4 | 0.716 | 0.009 | 0.88 |
| 2 | MUC16 | 1 | 0.00008 | 0.151 | 0.76 |
| 3 | LRP1B | 5 | 0.284 | 0.472 | 0.71 |
| 4 | COL11A1 | 3 | 0.100 | 0.0003 | 0.65 |
| 5 | BRSK2 | 6 | 0.267 | 0.0001 | 0.64 |
| 6 | XDH | 9 | 0.503 | 0.003 | 0.63 |
| 7 | CFTR | 2 | 0.044 | 0.001 | 0.61 |
| 8 | ACADL | 7 | 0.153 | 0.00002 | 0.55 |
| 9 | NOS1 | 8 | 0.209 | 0.0007 | 0.53 |
| 10 | NPTX2 | 11 | 0.054 | 0.501 | 0.52 |

The CSV contains all 15 targets, along with their nearest drugged neighbors.

## 5. Biological interpretation

### Re-ranking
- **The prior and Module 2 are close to independent** (Spearman ρ = -0.14, p = 0.62, n = 15). The
  prior adds new information, but it reshuffles the list rather than confirming it.
- **Biggest drop: CFTR (rank 2 to 7).** Both of its priors are close to zero, even though ivacaftor
  is an approved small molecule drug that acts on CFTR. CFTR is not in the training panel, so the
  prior never saw that label. Ivacaftor also binds unusual sites within the transmembrane domains
  rather than a classic enzyme pocket.
- **Biggest rises (+3 each): CEL (4 to 1), XDH (9 to 6), and FN1 (15 to 12).** XDH is a useful
  sanity check: it is the validated target of allopurinol and febuxostat, and the prior identifies
  it as small molecule tractable without having seen its label. FN1's p_sm of 0.482 is less
  convincing, because fibronectin's clinical uses are antibody based (EDB domain antibody fusions).

### Why w = 0.6
1. The Module 2 score is the disease specific evidence. The prior knows nothing about PDAC.
2. The family split AUROC of 0.73 is informative but far from decisive.
3. The top 15 Module 2 scores are compressed into a narrow range (0.932 to 0.979), so the prior
   is useful for breaking near ties.

**Sensitivity check:** CFTR is the most weight sensitive target. It ranks 3rd at w = 0.8 and 12th at
w = 0.4, which shows that the prior's main effect on this list is its judgment of CFTR.

### Modality
CEL is a lipid digesting enzyme with a catalytic site and a substrate binding pocket, so its strong
small molecule preference makes biological sense. Being secreted does not automatically make it a
good antibody target. MUC5B, a large extracellular mucin with no catalytic pocket, favors
antibodies instead. Together these suggest that structure and function, along with localization,
drive the modality call. MUC16's low antibody prior is a limitation, because antibody based
approaches against MUC16 already exist, and sequence truncation may explain the low score.

## 6. Variants, limitations, and Week 4 choice

### Variants
Sequence model scores do not predict disease risk. APOE ε2 scored strongly negative even though it
lowers Alzheimer's risk, and ε4 scored positive even though it raises risk. The score measures how
well a substitution fits learned evolutionary patterns, not whether it causes harm. TREM2 R47H and
TP53 R175H are clearly disease relevant but did not receive extreme scores. Percentiles show a
variant's rank within its own protein, not the probability of disease.

### Limitations
- Embeddings suggest tractability but cannot confirm an accessible binding pocket, or show that
  targeting a protein would benefit PDAC.
- Predictions may reflect historical research bias in which proteins have been drugged.
- Very long proteins (MUC16, MUC5B, LRP1B) are embedded from truncated sequences.
- Only the any-clinical prior was validated, not the separate small molecule and antibody priors.
- Module 2 compared against separate normal cohorts, so strong signals may reflect changes in
  tissue composition rather than tumor dependencies.

### Week 4 targets: CEL and XDH
- **CEL** ranked first overall (final 0.88, p_sm 0.716). Its catalytic pocket is a clear hypothesis
  to test through pocket geometry, accessibility, and local structural confidence.
- **XDH** ranked sixth and is a validated small molecule target. It gives a second enzyme for
  assessing catalytic pockets and the influence of cofactors and domain arrangement.

A good pocket would support structural tractability, but functional experiments would still be needed to decide whether inhibiting them is appropriate in PDAC.

## 7. Files in this folder

| File | Description |
|---|---|
| `BIOT6900_Module3_Assignment.ipynb` | Completed notebook covering target mapping, the leakage guard, prior training and validation, modality priors, nearest drug targets, anchor variant scoring, and re-ranking. Runs top to bottom with Restart & Run All |
| `targets_pdac_w3.csv` | Re-ranked target table with Module 2 rank, modality priors, final score, and nearest drugged neighbors. This is the hand-off to Week 4 |
| `Module 3 Assignment Report.pdf` | Written report covering re-ranking, modality, validation, anchor variants, limitations, and the Week 4 choice |
| `README.md` | This file |

## 8. How to reproduce

1. Clone the repo and open the Module 3 folder:

        git clone https://github.com/Deepthi-Dineshkumar/biot6900.git
        cd "biot6900/Module 3"

2. Install the required packages:
   
        pip install pandas numpy scipy scikit-learn matplotlib mygene jupyter

4. Place the course data pack (`proteome.tsv`, the embedding files, and the training panel) in location. The data pack is not committed to this repo.
5. `../module2/targets_pdac.csv` exists already in the module 2 folder.
6. Start Jupyter from inside `Module 3/`, open `BIOT6900_Module3_Assignment.ipynb`, and choose Kernel, then Restart & Run All.

**Expected results:** 15 of 15 targets mapped, 0 panel proteins removed, family split AUROC 0.732 and AP 0.548, CEL ranked first with a final score of 0.88, and `targets_pdac_w3.csv` written to this folder.
