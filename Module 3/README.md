\# BIOT 6900 · Module 3: Target Triage with Protein Language Models



\*\*Name:\*\* Deepthi Dineshkumar

\*\*Disease:\*\* Pancreatic ductal adenocarcinoma (PDAC), CPTAC-PDAC cohort via LinkedOmics (continued from Module 2)



\## Files in this folder



| File | Description |

|---|---|

| `BIOT6900\_Module3\_Assignment.ipynb` | Completed assignment notebook (runs top to bottom with Restart \& Run All) |

| `targets\_pdac\_w3.csv` | Re-ranked target table, the hand-off to Week 4 |

| `BIOT6900\_Module3\_Report\_Dineshkumar.pdf` | Report answering the five Step 8 questions |


\*\*Input:\*\* `module2/targets\_pdac.csv`, the Module 2 output (already committed to the repo).



\## Notes and things that didn't work as expected



\- \*\*Gene IDs:\*\* my Module 2 file stores Ensembl gene IDs rather than gene symbols. I added a mapping step so the IDs match `proteome.tsv`. All 15 top targets mapped, so there were no unmapped genes.

\- \*\*Leakage guard:\*\* 0 panel proteins were removed, because none of my targets appear in the 1,201-protein training panel. The guard is still applied.

\- \*\*Validation scope:\*\* only the any-clinical prior was cross-validated. The small-molecule and antibody priors were not validated separately.

\- \*\*Known model limitations noted in the report:\*\* CFTR receives a very low prior despite having approved small-molecule drugs. Very long proteins (MUC16, MUC5B, LRP1B) are embedded from truncated sequences.





