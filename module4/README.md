# BIOT6900 Module 4: Structure Based Druggability of PDAC Targets

**Name:** Deepthi Dineshkumar

**Disease:** Pancreatic Ductal Adenocarcinoma (PDAC), (carried forward from Modules 2 and 3)
**Full write up:** `Deepthi_Dineshkumar_Module4_Report.pdf`

## Overview
Module 3 suggested that CEL and XDH look small molecule tractable based on sequence alone. This module tests whether their predicted 3D structures actually contain credible drug binding pockets.
For each target AlphaFold models were used, their confidence scores, experimental structures, and fpocket. This analysis tests structural tractability only. It does not show whether inhibiting either protein would benefit PDAC.

## 1. Targets chosen and why

| Gene | Protein | UniProt | Module 2 rank | Module 3 rank | p_sm | Experimental reference |
|---|---|---|---|---|---|---|
| CEL | Carboxyl ester lipase | P19835 | 4 | 1 | 0.716 | 1F6W chain A |
| XDH | Xanthine dehydrogenase/oxidase | P47989 | 9 | 6 | 0.503 | 2E1Q chain A |

- **CEL** ranked first after the Module 3 re-ranking and has the highest small molecule prior. It is a compact serine hydrolase, so whether its catalytic pocket is credible was the open question.
- **XDH** is the **positive control**. Its active site is the validated target of two approved drugs (allopurinol and febuxostat), so a working pipeline should find a strong pocket there.

**Inputs:** AlphaFold models and PAE files for both targets (AF-P19835 and AF-P47989) from the AlphaFold Protein Structure Database, and experimental structures 1F6W, 2E1Q, and 1TUP from the RCSB PDB.

## 2. pLDDT, PAE, and trustworthy regions

| Target | Mean pLDDT | Above 90 | Below 50 | Stretches below 50 (10+ residues) |
|---|---|---|---|---|
| CEL | 82.6 | 66.3% | 16.2% | 3 to 17; 555 to 592; 717 to 753 |
| XDH | 96.0 | 93.4% | 1.7% | 168 to 190 |

**Domains found by the PAE based finder:**
- **CEL:** one confident group, residues **21 to 552** (532 residues). This is the trusted region for pocket analysis. The low confidence N terminus overlaps the signal peptide, and the low confidence C terminal tail is the glycosylated, repeat rich mucin like region, which a single unglycosylated model cannot represent well.
- **XDH:** one confident group, residues **2 to 1333** (1,300 confident residues). Nearly the whole model is usable. Residues 168 to 190 need caution, but this loop is a known conformational switch (its cleavage converts the dehydrogenase form into the oxidase form), not a modelling error.

Because each protein produced only one confident group, inter domain PAE could not be calculated.`inter_pae` is NaN for both targets in the summary CSV. This is the expected output for a single group, not a failure.

## 3. Comparison with experimental structures

Residues were paired by Needleman Wunsch sequence alignment. Only identical residue pairs with
AlphaFold pLDDT of at least 70 were used for C alpha superposition. Both proteins formed a single
confident group, so no separate per domain superposition was needed.

| Target / reference | Confident identical pairs | RMSD | Largest deviations |
|---|---|---|---|
| CEL / 1F6W A | 527 | 1.35 Å | Residues 444 to 449, maximum 9.1 Å at residue 447 |
| XDH / 2E1Q A | 1,297 | 0.48 Å | 3 (5.6 Å), 537 (3.6 Å), 538 (3.4 Å), 1288 (2.7 Å) |

- **CEL:** agreement is below 1 Å almost everywhere, including the catalytic triad. The one large difference is a six residue surface loop far from the active site. Surface loops in 1F6W are known to be mobile, and the crystal has engineered surface mutations to help crystallization, so this looks like a conformational difference rather than a prediction error. 1F6W covers only the catalytic domain, so it cannot validate the C terminal tail.
  
- **XDH:** very close agreement overall. 2E1Q contains bound cofactors and the ligand SAL, while the AlphaFold model is a single chain with no bound components. The differences mostly reflect what the crystal has that the model lacks.

**Numbering and mutations:**
- 1F6W uses mature protein numbering, which is offset by CEL's 20 residue signal peptide. All
  positions in this README and the report use UniProt numbering, and the sequence alignment handles the conversion.
- 2E1Q is the E803V mutant of human XOR. Because the superposition keeps only identical residue pairs, position 803 is excluded automatically.

## 4. Pocket detection and confidence filtering

I ran fpocket on both AlphaFold models and applied these confidence rules:
- **Reject** a pocket if its lining residues have mean pLDDT below 70.
- **Suspect** if mean internal PAE is above 10 Å.
- **Keep** everything else.

| Target | Rank | Pocket | Druggability score | Mean pLDDT | PAE (Å) | Verdict |
|---|---|---|---|---|---|---|
| CEL | 1 | 56 | 0.925 | 93.7 | 2.1 | keep |
| CEL | 2 | 60 | 0.562 | 51.1 | 13.4 | reject |
| CEL | 3 | 3 | 0.193 | 98.5 | 1.2 | keep |
| CEL | 4 | 62 | 0.135 | 95.6 | 2.2 | keep |
| CEL | 5 | 11 | 0.093 | 98.2 | 1.2 | keep |
| XDH | 1 | 2 | 0.979 | 98.2 | 1.5 | keep |
| XDH | 2 | 1 | 0.933 | 96.6 | 2.1 | keep |
| XDH | 3 | 85 | 0.176 | 98.2 | 1.4 | keep |
| XDH | 4 | 3 | 0.074 | 96.1 | 1.8 | keep |
| XDH | 5 | 84 | 0.060 | 97.3 | 1.6 | keep |

### Experimental ligand site check (XDH only)
For each ligand bound to 2E1Q chain A, I mapped contact residues (4.5 Å cutoff) onto UniProt
numbering and checked how much of each site the detected pockets cover.

| Ligand (chain A) | Mapped site residues | Coverage by pocket 2 | Best matching pocket (rank) | Best coverage |
|---|---|---|---|---|
| DY1 A7003 | 27 | 74% | 2 (rank 1) | 74% |
| FAD A2003 | 39 | 0% | 1 (rank 2) | 90% |
| SAL A2006 | 12 | 92% | 2 (rank 1) | 92% |

1F6W contains no eligible bound small molecule, so this check could not be run for CEL.

## 5. Summary table exported

The per target summary is exported to **`structure_assessment_pdac.csv`**, with one row for each target. It records the trusted region, the experimental comparison (reference, residue pairs,RMSD), the inter domain PAE (NaN for both, as explained above), the best kept pocket and its confidence values, and the final druggability call.

| Target | Best kept pocket | Call |
|---|---|---|
| CEL | Pocket 56 (score 0.925) | Possible but uncertain |
| XDH | Pocket 2 (score 0.979) | Structure supports a small molecule pocket |

## 6. Druggability calls interpreted

### XDH: structure supports a small molecule pocket
Pocket 2 has the highest druggability score (0.979), high lining residue confidence, low PAE, and covers 92% of the residues that contact the bound ligand SAL. The best pocket matching each experimental site ranked first or second, so there is no sign that the scoring failed. As expected for the positive control, the pipeline found the known drug site. One caution is that SAL contacts residue 803, the engineered E803V position, which adds uncertainty when transferring this site to the wild type protein. DY1 and FAD are cofactors, so their sites should not be treated as interchangeable drug binding sites.

### CEL: possible but uncertain
Pocket 56 scores strongly and passes both confidence checks. The next kept pocket drops sharply to
0.193, so this is a single well defined site rather than a scatter of weak cavities. Pocket 60 shows
why geometry alone is not enough: it scored second but was rejected, because its lining residues
fall in the low confidence C terminal region (pLDDT 51.1, PAE 13.4 Å). The call stays "possible but
uncertain" for three reasons:
- There is no experimental ligand to check pocket 56 against.
- CEL's active site is partly covered by a surface loop until bile salt binds, and neither the model
  nor 1F6W shows that activated state.
- The model lacks glycosylation, which may affect pocket accessibility.

## 7. p53 anchor check, limitations, and integration

### p53 anchor check
The AlphaFold p53 model matched 1TUP chain B with an RMSD of 0.68 Å over 192 confident residues.
- **Zinc coordination:** C176, H179, C238, and C242 coordinate zinc at 2.1 to 2.3 Å.
- **DNA contacts:** R248 and R273 contact DNA (3.4 Å and 2.5 Å). R249 lies nearby at 7.2 Å.
- **Neither cutoff:** R175, G245, Y220, and R282 meet neither the zinc nor the DNA contact cutoff.

These results help explain how mutations can damage p53 through zinc dependent stability or DNA
binding. The 19th percentile Module 3 score for R175H shows that sequence based scores can miss
indirect destabilization near the zinc site.

### Limitations
- **Missing biology:** the models lack CEL glycosylation and XDH cofactors and partner proteins, which may change pocket accessibility.
- **Single conformation:** only one predicted conformation was assessed. Multiple structures or molecular dynamics would show whether the pockets stay stable.
- **Reference limits:** CEL was validated only against a catalytic domain construct, and the XDH reference carries an active site mutation.
- **Confidence is not binding:** high pLDDT does not confirm that a drug will bind or that the protein will be functionally affected.
- 
### Integration with Modules 2 and 3
Both targets were retained, and the structural results support the Module 3 small molecule
predictions. XDH has stronger structural and experimental support even though CEL had the higher Module 3 prior, so XDH is prioritized for binding studies. CEL needs further testing to see whether pocket 56 is functional, but weaker validation alone does not justify switching to an antibody approach. Both genes were **downregulated** in Module 2, so a druggable pocket does not show that inhibition would help in PDAC. Experiments in PDAC and normal pancreatic models are needed to decide whether either target should be inhibited or restored.

## 8. Files in this folder

| File | Description |
|---|---|
| `BIOT6900_Module4_Assignment.ipynb` | Completed notebook covering structure download, pLDDT and PAE analysis, domain finding, experimental superposition, pocket filtering, the ligand site check, and the p53 anchor check |
| `structure_assessment_pdac.csv` | Exported per target summary table |
| `Deepthi_Dineshkumar_Module4_Report.pdf` | Written report |
| `README.md` | This file |

Downloaded structures and fpocket output folders are not committed. They are regenerated when the
notebook is run (see the fpocket note below).

## 9. How to reproduce and what didn't work

1. Clone the repo and open the Module 4 folder:

        git clone https://github.com/Deepthi-Dineshkumar/biot6900.git
        cd "biot6900/Module 4"

2. Activate the course environment (`conda activate biot6900`).
3. Start Jupyter from inside this folder, open the notebook, and choose Kernel, then Restart & Run All. An internet connection is needed to download the AlphaFold and PDB files.

**fpocket would not install on Windows.** I ran Step 6 on the Northeastern Explorer cluster. I uploaded the two AlphaFold models, ran `fpocket -f` on each, and copied the `AF-P19835_out` and `AF-P47989_out` folders back into `structures/`. The Step 6 cell checks for these folders and reuses them when they exist, so the notebook still runs top to bottom on Windows. On Linux or macOS with fpocket installed, the cell runs fpocket directly.
