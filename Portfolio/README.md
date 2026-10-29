# BIOT6900 Portfolio Deliverable 1: Multi-Omics Target Discovery Package (PDAC)

**Name:** Deepthi Dineshkumar

**Disease:** Pancreatic Ductal Adenocarcinoma (PDAC)

**Full write up:** `Deepthi_Dineshkumar_PD1_Report.pdf`

## Overview

This folder packages Modules 2 to 4 into one whole argument, from multi-omics evidence, through druggability triage, to structural assessment and a ranked shortlist. **No primary lab analyses were rerun.** `pd1_figures.py` reads the outputs already committed (or regenerated) in the Module 2, 3 and 4 folders and summarises them as the four required figures and the final prioritization table. The only additions are a Benjamini-Hochberg FDR adjustment of the Module 2 p values and explicitly stated evidence thresholds, used for the evidence Venn (Figure 1c).

Every target is judged on three criteria that run through the whole report:

| Criterion | Question | Evidence comes from |
|---|---|---|
| C1, direction | Is the gene increased in tumour at both RNA and protein? This supports investigating an inhibitory approach; functional experiments are needed to show the tumour depends on it. | Module 2 |
| C2, modality | Is a modality plausible: best prior (`p_sm` or `p_ab`) of at least 0.4, or prior clinical investigation for this gene? | Module 3 |
| C3, structure | Does the best pocket have the Module 4 verdict `keep` (lining pLDDT of 70 or more, internal PAE of 10 Å or less) **and** an fpocket druggability score of 0.5 or more? | Module 4 (CEL and XDH only) |

## 1. Final shortlist

| Rank | Target | Route | C1 | C2 | C3 | Why |
|---|---|---|---|---|---|---|
| 1 | **XDH** | Small molecule (repurposing) | Yes (RNA +2.79, protein +0.84) | Yes (p_sm 0.503; allopurinol, febuxostat) | Yes (pocket 2, score 0.979, covers 92% of the 2E1Q SAL site) | Only gene meeting all three criteria |
| 2 | **MUC16** | Antibody or ADC | Yes (RNA +3.70, protein +1.12) | Yes, by prior clinical investigation (p_ab only 0.151, truncated sequence) | Not assessed (14,507 residues) | Strongest disease evidence; missing criterion is a data gap |
| 3 | **CEL** | Hold until direction is resolved | No (RNA −6.35, protein −2.55) | Yes (p_sm 0.716) | Yes, uncertain (pocket 56, score 0.925, no ligand to check) | Strong pocket, but decreased in bulk tumour; cellular origin and therapeutic direction unresolved |

**Watch list:** MUC5B (antibody prior 0.55, but a secreted mucin would act as a soluble antigen sink) and FN1 (prior 0.48, but its clinical precedence is isoform specific). All 15 triaged genes, with their C1 to C3 flags and a one line rationale, are in `prioritized_targets_pdac.csv`.

**Tie-break rule.** When genes meet the same number of criteria, a gene missing a criterion because of a *data gap* (no structure yet) ranks above one missing a criterion because of a *biological gate* (lower in tumour). More modelling can fill a data gap; a pocket cannot settle whether inhibiting a protein that is decreased in tumour would help. This is why MUC16 ranks above CEL.

## 2. Interpretation notes

### Integration (Figure 1)

- **1,880 genes** survived the three-layer join. The limit was the protein and mutation overlap (only 1,882 of 3,362 mutated genes were measured as protein), not the RNA layer.

- **1,190 of 1,880 genes (63%)** were concordant between RNA and protein. Among the 153 genes with strong evidence at both layers, 140 (92%) were concordant.

- **The mutation layer is sparse.** 82% of the 1,880 genes are mutated in a single tumour out of 105, and only 14 reach 5%. The genomic percentile is mostly tied, so RNA and protein drive the ranking. **MUC16 is the only gene meeting all three selected evidence thresholds within the 1,880-gene intersection.** That is a statement about these thresholds, not about evidence in general.

- **Known genes recovered:** MUC16 (rank 1), COL11A1 (3) and KRT19 (10) at the top; drivers CDKN2A (32), KRAS (69), TP53 (381) and SMAD4 (530) lower down; RNF43 lost because it is missing from the proteomics file.

- **Discordant gene examined:** TP53 (RNA −0.05, protein +0.80, mutated in 83% of tumours), consistent with stabilised mutant p53.

- **Non-obvious hit:** XDH, which is not a known PDAC target but is up in tumour at both layers and ends up ranked first.

- **Clinical precedence:** the anti-MUC16 ADC phase I study (Liu et al., 2016; NCT01335958) included pancreatic cancer. This shows the modality has been investigated, not that it benefits PDAC; responses were reported only in ovarian cancer.

### Venn diagrams

- **Figure 1b, coverage:** gene IDs present in each layer file in `module2/data/`. Circle areas are not to scale, because the layers differ in size by more than ten times.

- **Figure 1c, evidence:** within the 1,880, RNA |log2FC| of 1 or more and protein |log2FC| of 0.5 or more, each at Benjamini-Hochberg FDR below 0.05, and mutation in at least 5% of tumours. **The FDR is adjusted across every gene tested in each layer** (60,579 RNA and 11,057 protein tests), not only the 1,880 retained genes. Adjusting within the 1,880 instead would give 262 RNA-evidence genes rather than 261; nothing else in the figure or ranking changes. These thresholds are set at the top of `pd1_figures.py` (`RNA_LFC`, `PROT_LFC`, `FDR`, `MUT_FREQ`). The region counts for both diagrams are saved in `venn_counts.csv`.

### Structure figure (Figure 2)

Models are drawn as a C-alpha trace viewed down the axis of the best kept pocket, with fpocket lining atoms as spheres and, for XDH, the 2E1Q SAL contact residues ringed. The script colours residues by **per-residue pLDDT** (AlphaFold bands) when the AlphaFold downloads `module4/structures/AF-<acc>.pdb` exist. If only the fpocket copy is available (fpocket overwrites the B-factor column), it falls back to the **Module 4 confidence regions** (PAE-confident domain, stretches below pLDDT 50, everything else). The figure in the submitted PDF uses the per-residue pLDDT version. The script prints which mode it used.

## 3. Files in this folder

| File | Description |
|---|---|
| `Deepthi_Dineshkumar_PD1_Report.pdf` | The 6-page package submitted on Canvas |
| `pd1_figures.py` | Builds all figures and tables below from the Module 2 to 4 outputs |
| `figures/fig1_multiomics_integration.png` | RNA vs protein scatter of the 1,880 genes, coverage Venn, evidence Venn |
| `figures/fig2_structures.png` | CEL and XDH AlphaFold models with kept pockets |
| `figures/fig3_prioritization.png` | Evidence matrix and package call for the 15 triaged genes |
| `prioritized_targets_pdac.csv` | Final table: Module 2 and 3 ranks, fold changes, priors, structure call, C1 to C3 flags, package call, rationale |
| `venn_counts.csv` | Region counts behind Figures 1b and 1c |
| `README.md` | This file |

**Inputs it reads** (paths relative to the repo root):

| Input | Produced by |
|---|---|
| `module2/targets_pdac.csv` | Module 2 notebook (committed) |
| `module2/data/pdac_transcriptomics.tsv`, `pdac_proteomics.tsv`, `pdac_mutation.tsv` | `module2/prep_pdac.py` (not committed, regenerate it) |
| `Module 3/targets_pdac_w3.csv` | Module 3 notebook (committed) |
| `module4/structure_assessment_pdac.csv` | Module 4 notebook (committed) |
| `module4/structures/AF-P19835*`, `AF-P47989*` and their `_out/pockets/` folders | Module 4 notebook and fpocket (not committed, regenerate them) |

The XDH SAL site residues and the CEL/XDH confidence regions are copied into the top of `pd1_figures.py` from the printed Module 4 output, with comments saying where they came from.

## 4. How to reproduce

1. Clone the repo:

        git clone https://github.com/Deepthi-Dineshkumar/biot6900.git
        cd biot6900

2. Activate the course environment. The script needs only `pandas`, `numpy` and `matplotlib`:

        conda activate biot6900

3. Rebuild the inputs that are not committed:

   - **Module 2 layers:** download the five CPTAC files listed in `module2/README.md` into `module2/`, then run `python prep_pdac.py` from inside `module2/`. This writes `module2/data/`.
   - **Module 4 structures:** open `module4/BIOT6900_Module4_Assignment.ipynb` and run all cells. This downloads the AlphaFold models into `module4/structures/` and runs fpocket (see the fpocket note in the Module 4 README if you are on Windows).

4. Run the script from the repo root (it finds the other folders on its own, so the working directory does not matter):

        python Portfolio/pd1_figures.py

**Expected output:**

    Fig 1: 1880 genes, 1190 concordant; coverage centre 1880; evidence RNA&protein 153 (140 concordant); all three: ['MUC16']
    Fig 2: structures coloured by per-residue pLDDT
    Fig 3 / CSV: shortlist = XDH, MUC16, CEL

If an input is missing, the script stops and lists exactly which file it could not find. If the AlphaFold models are missing, Figure 2 shows a note in place of that panel and the rest still runs.

## 5. What didn't work or was not possible

- **MUC16 has no structure assessment.** At 14,507 residues it is far beyond the 2,700-residue limit of the Module 4 workflow. It is ranked on Modules 2 and 3 plus clinical precedence, and the report flags this as a data gap.

- **CEL has no experimental ligand check**, because 1F6W has no bound small molecule. Its pocket call stays "possible but uncertain".

- **The Module 3 tractability priors missed two known drug targets** (CFTR and MUC16), so prior clinical investigation was allowed to satisfy C2. The precedence used for each gene is listed in `PRECEDENCE` in the script and in the `clinical_precedence` column of the CSV.

- **The package calls for MUC5B and FN1 are judgement calls**, not automatic outputs. They are listed in `CALLS` in the script with their reasons, so they can be checked or changed.

