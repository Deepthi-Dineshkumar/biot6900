# BIOT 6900 · Module 5: Variant Interpretation on CYP2C19

**Name:** Deepthi Dineshkumar

**Gene:** CYP2C19 (no gene was assigned on Canvas, so CYP2C19 was chosen from the five options)

**AlphaGenome Atlas:** NOT run, and `SKIP_ATLAS` was set to `True`.
Section 2 of the report addresses what the Atlas could in principle have shown for
rs72558186 (CYP2C19*7) and rs4244285 (CYP2C19*2).

## Files

| File | Description |
|---|---|
| `BIOT6900_Module5_Part3_Assignment.ipynb` | Completed notebook (Restart & Run All, clean top to bottom) |
| `pgx_classified_CYP2C19.csv` | Q5 export, 33 VUS with evidence and provisional call |
| `Deepthi_Dineshkumar_Module5_Report.pdf` | Report answering the four Q6 questions |

`data/` is in `.gitignore` and is not committed. The notebook recreates it from MyVariant.info.

## What didn't work, and what was done about it

- **CYP2C19 cannot be calibrated on its own.** After the frequency filter the gene has
  **0 P/LP and 5 B/LB** labelled variants, so the first calibration returned a single
  populated band with `LR+ = 0.00` and four empty bands. The cause is structural: almost
  all CYP2C19 ClinVar submissions are filed as *drug response*, which is not a
  pathogenicity assertion, so the pathogenic side is empty by construction. The second
  route the assignment permits was taken instead: **pooling with the other four
  pharmacogenes** (616 variants, 162 labelled — 47 P/LP, 115 B/LB). The cost is disclosed
  in the notebook and argued in Section 1 of the report: DPYD supplies 127 of the 162
  labelled variants, so the likelihood ratios are not CYP2C19-specific.
- **Two bands are not usable even after pooling.** The 25–30 band rests on 6 variants, and
  the 30–50 band reports `LR+ = inf` because its benign count is zero. A rule was added to
  Q5 so that any VUS whose band rests on fewer than 10 labelled variants stays a VUS;
  3 variants were held back this way.
- **6 VUS could not be scored by any tool.** They are large structural deletions
  (e.g. `chr10:(?_96501534)_(96548712_?)del`, ClinVar 689097) with no rsID and no
  consequence call. CADD, REVEL and AlphaMissense are all SNV-level tools. These are
  reported separately in Q5 as "remains VUS (no score)", not as benign.
- **rs72558186 (CYP2C19*7) returned no score from any tool**, despite being a canonical
  splice-donor loss-of-function allele. This is the named unscored variant discussed in
  Section 2 of the report.
- **Fetch size.** The query returned 80 variants for CYP2C19. Both anchor alleles
  (rs4244285 *2 and rs12248560 *17) were present, so the fetch is not truncated, but the
  number is smaller than for DPYD (426) and is the reason the single-gene calibration fails.
- **A `RequestsDependencyWarning`** about `charset_normalizer` appeared on import. It did
  not affect the fetch; all downloads completed and the anchor checks passed.
- **A pandas `TypeError`** occurred in Q5 on the first run, because `.map()` on a categorical
  band column returns a categorical that cannot be compared to an integer. Fixed by casting
  the mapped columns with `.astype(object)` and `pd.to_numeric`.
