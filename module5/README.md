# BIOT6900 Module 5: Variant Interpretation in CYP2C19

**Name:** Deepthi Dineshkumar
**Gene:** CYP2C19
**Full write up:** `Deepthi_Dineshkumar_Module5_Report.pdf`
**AlphaGenome Atlas:** not run, `SKIP_ATLAS` was set to `True` (see Section 6).

## Overview
CYP2C19 metabolizes widely used drugs, including clopidogrel, several antidepressants, and proton
pump inhibitors. In this module, its ClinVar variants were fetched, audited which variants the
prediction tools can score, applied a frequency filter, calibrated CADD into ACMG evidence bands,
and used those bands to give the gene's variants of uncertain significance (VUS) a provisional
classification. CYP2C19 is a useful test case because its most important alleles are common
and are labelled "drug response" rather than pathogenic. Both of these features break
assumptions built into standard variant interpretation.

## 1. Gene fetched, cached, and quality checked

Variants were fetched from MyVariant.info and cached in `data/`. The fetch returned **80 CYP2C19
variants**:

| ClinVar label | Count |
|---|---|
| Drug response | 38 |
| VUS | 33 |
| Benign / likely benign (B/LB) | 9 |
| Pathogenic / likely pathogenic (P/LP) | 0 |

**Quality checks:**
1. **Anchor alleles present:** CYP2C19*2 (rs4244285, gnomAD AF 0.168) and CYP2C19*17
   (rs12248560, AF 0.200) were both returned.
2. **Fetch not truncated:** 80 variants is well below the fetch limit, and both anchors arrived,
   so the result is complete. The gene simply has fewer ClinVar entries than DPYD (426).
3. **Label breakdown inspected:** there are no P/LP variants, because almost all CYP2C19
   submissions are filed as drug response. This is what later prevents single gene calibration.

## 2. Coverage audit by consequence class

| Tool | Variants scored (of 80) | What it covers |
|---|---|---|
| **CADD** | **68** | Most classes, including 6 of 7 synonymous and 4 of 5 upstream variants |
| REVEL | 54 | Missense only |
| AlphaMissense | 54 | Missense only |

**Unscoreable classes:**
- **Structural deletions (8):** no tool scored them. CADD, REVEL, and AlphaMissense all work on
  single nucleotide variants, and these large deletions (for example
  `chr10:(?_96501534)_(96548712_?)del`, ClinVar 689097) have no rsID and no consequence call.
- **Splice donor:** rs72558186 (CYP2C19*7, c.819+2T>A) received no score from any tool, even
  though it is a known nonfunctional allele.
- **Other non missense classes:** REVEL and AlphaMissense left 26 variants unscored, including
  synonymous, upstream, intronic, intergenic, stop gained, and stop lost variants.

**CADD** was chosen as the predictor because it has the broadest coverage of the variant types that
matter for CYP2C19.

## 3. Frequency filter

**Threshold: MAX_AF = 0.25.** A rare disease cutoff does not suit CYP2C19, because its clinically
important alleles are common (*2 at 0.168 and *17 at 0.200). The 0.25 threshold kept both anchors
and removed only **2 of 80 variants**. When I inspected them, neither was labelled P/LP or drug
response. In this analysis, frequency is a screening step, not evidence that a variant has
normal function.

**The 0.001 contrast:** a rare disease cutoff of 0.001 would have removed **19 variants**,
**14 of them labelled drug response or P/LP**. These include *2, *3, *9, and the regulatory
variants rs12769205 and rs12777823. A strict rarity filter would throw away the best
characterized pharmacogenetic evidence while keeping poorly understood rare variants.

## 4. Calibration: the small gene route

**Single gene calibration failed.** After filtering, CYP2C19 has **0 P/LP and 5 B/LB** labelled
variants. Calibration returned one populated band with LR+ = 0.00 and four empty bands. This
reflects the lack of positive examples, not strong benign evidence. Drug response labels are not
pathogenicity labels, so the pathogenic side is empty by construction.

**Route taken: pooling with the other four pharmacogenes.** This gives 616 variants, of which
162 are labelled (47 P/LP, 115 B/LB), with MAX_AF = 0.25 applied throughout.

| CADD band | n | P/LP | B/LB | LR+ | ACMG evidence strength |
|---|---|---|---|---|---|
| 0 to 10 | 72 | 1 | 71 | 0.03 | Strong, benign |
| 10 to 20 | 33 | 0 | 33 | 0.00 | Strong, benign |
| 20 to 25 | 17 | 8 | 9 | 2.17 | Supporting, pathogenic |
| 25 to 30 | 6 | 4 | 2 | 4.89 | Not used (too few variants) |
| 30 to 50 | 34 | 34 | 0 | inf | Not reliable (no benign variants) |

Labelled variants contributed by each gene: DPYD 88 B/LB + 39 P/LP, VKORC1 11 + 7, CYP2C9 7 + 1,
CYP2C19 5 + 0, TPMT 4 + 0.

**Sample size engaged with:**
- **DPYD dominates.** DPYD supplies 127 of the 162 labelled variants, while CYP2C19 supplies only
  5 benign ones. The likelihood ratios are therefore not CYP2C19 specific. Applying them assumes
  CADD performs the same across different enzymes, which the data cannot test.
- **The 25 to 30 band rests on 6 variants**, so its LR+ of 4.89 is unstable.
- **The 30 to 50 band's infinite LR+** comes from having zero benign variants, not from strong evidence.
- **Safeguard added:** any VUS in a band with fewer than 10 labelled variants stays a VUS.

## 5. VUS classification and export

All 33 VUS were matched to a calibrated band and exported to **`pgx_classified_CYP2C19.csv`**,
along with their scores, band, evidence strength, and provisional call.

| Provisional outcome | Count |
|---|---|
| Shifted toward benign | 16 |
| Shifted toward pathogenic | 8 |
| Held as VUS (band has fewer than 10 labelled variants) | 3 |
| Remains VUS (no score, structural deletions) | 6 |
| **Total** | **33** |

All 8 pathogenic shifts rely on the CADD 20 to 25 band (supporting strength), which comes from
pooled, DPYD dominated calibration. These are **provisional computational calls, not clinical
reclassifications.** Unscored variants are reported separately and are never treated as benign.

## 6. Biological interpretation

### Thresholds
Pharmacogenes break the "rare equals damaging" logic. Common CYP2C19 alleles such as *2 change
drug metabolism in a large share of the population, so frequency filters must be loose enough
to keep them.

### Coverage gaps and the regulatory genome
- **Missing scores are not benign evidence.** CYP2C19*7 disrupts a canonical splice donor site and
  is a known nonfunctional allele, but no tool scored it. SpliceAI or Pangolin could predict the
  donor loss. A minigene assay, or RT-PCR and RNA-seq in a suitable model, could confirm abnormal
  transcripts and show whether nonsense mediated decay reduces transcript levels.
- **CYP2C19*2 shows the same blind spot.** It is annotated as synonymous with a CADD score of only
  5.686, yet it causes abnormal splicing.
- **AlphaGenome was not run.** Its splicing outputs could have linked *2 and *7 to their mechanism,
  and its regulatory outputs could have helped interpret noncoding variants. No Atlas results were
  generated, so these mechanisms are not claimed as model findings.

### Discordant variant: rs41291556 (CYP2C19*8, p.Trp120Arg)
This variant was flagged because REVEL and AlphaMissense differed by 0.50, above the 0.4 threshold.

| Tool / field | Value | Reads as |
|---|---|---|
| CADD phred | 25.70 | Top ~0.3% genome wide |
| REVEL | 0.440 | Below the usual 0.5 threshold, tolerated |
| AlphaMissense | 0.9402 | Likely pathogenic |
| gnomAD AF | 0.001571 | Rare |
| ClinVar | Drug response | Functionally characterized, not a P/LP assertion |

REVEL was trained to separate disease causing missense variants from neutral ones, which may not
capture effects on drug metabolism. Replacing a tryptophan with a charged arginine could plausibly
disrupt the enzyme, which fits the AlphaMissense and CADD scores. The same pattern appears for
R186P (*22) and F448L (*25), where REVEL scored 0.209 and 0.254 and AlphaMissense scored 0.800 and
0.798. This suggests a systematic calibration difference, but the scores cannot prove the cause or
the structural mechanism. It is another reason not to rely on REVEL alone.

## 7. Limits and next experiment

### What this analysis cannot establish
- **No genotypes or phase:** the data are population allele frequencies, so diplotypes and
  metabolizer phenotypes cannot be assigned, and it is unknown whether W120R and *17 occur on the
  same chromosome.
- **Phenoconversion:** co-medications, liver function, enzyme expression, and inflammation all
  change CYP2C19 activity, so these calls describe possible variant effects, not an individual's
  drug response.
- **ClinVar drug response labels** do not show how often carriers experience a clinical outcome.

## 8. Files in this folder

| File | Description |
|---|---|
| `BIOT6900_Module5_Part3_Assignment.ipynb` | Completed notebook (runs top to bottom with Restart & Run All) |
| `pgx_classified_CYP2C19.csv` | The 33 VUS with scores, band, evidence strength, and provisional call |
| `Deepthi_Dineshkumar_Module5_Report.pdf` | Report answering the four Q6 questions |
| `README.md` | This file |

The `data/` folder is in `.gitignore` and is not committed. The notebook recreates it from
MyVariant.info when it is run.

## 9. How to reproduce, and fallbacks

1. Clone the repo, open the Module 5 folder, and activate the course environment.
2. Open the notebook and choose Kernel, then Restart & Run All. An internet connection is needed
   for the MyVariant.info fetch. Leave `SKIP_ATLAS = True` unless you have an AlphaGenome API key.

**Fallbacks and issues:**
- **AlphaGenome skipped:** no API key, so the regulatory analysis is incomplete (see Section 6).
- **Pooled calibration:** single gene calibration was impossible, so the five pharmacogenes were
  pooled (see Section 4).
- **Small band safeguard:** VUS in bands with fewer than 10 labelled variants stay VUS (3 held back).
- **Unscoreable variants:** 6 structural deletion VUS are reported as "remains VUS (no score)", not as benign.
- **`RequestsDependencyWarning`:** a warning about `charset_normalizer` appeared on import. It did
  not affect the fetch, and the anchor checks passed.
- **pandas `TypeError` in Q5:** `.map()` on the categorical band column returned a categorical
  that could not be compared with an integer. I fixed this by casting the mapped columns with
  `.astype(object)` and `pd.to_numeric`.
