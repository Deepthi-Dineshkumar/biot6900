# Module 9 — Clinical NLP PDAC

**Course:** BIOT 6900 · AI-Powered Biotech Workflows · Computational Lab 9 (Fall 2026)
**Notebook:** `BIOT6900_Module9_Notebook3_Graded.ipynb`

## Disease
Pancreatic ductal adenocarcinoma (PDAC) in hospitalized patients, with in-hospital death as the outcome. PDAC continues the disease I have used across this course (Modules 2–7).

| Role | Concept | Synonyms used | Why |
|---|---|---|---|
| Signal (cause) | metastatic_disease | metastatic disease, liver metastases, peritoneal carcinomatosis, stage iv | Baseline disease burden, known from staging before admission |
| Signal (cause) | cachexia | cachexia, cachectic, significant weight loss, sarcopenia | Tumour-driven wasting that predicts poor survival |
| Signal (cause) | cholangitis | cholangitis, ascending cholangitis, biliary sepsis | Biliary obstruction from the tumour causes the deterioration |
| Consequence (leak) | vasopressors | vasopressors, norepinephrine, levophed | Given because the patient is in shock |
| Consequence (leak) | intubation | intubated, mechanical ventilation, on the vent | Started because of respiratory failure |
| Consequence (leak) | comfort_care | comfort care, comfort measures only, hospice, dnr/dni | Decided once death is expected |

## How to re-run
1. Activate the `biot6900` environment (only `pandas`, `scikit-learn`, `matplotlib` needed; no downloads — all data are synthetic and generated in the notebook).
2. Open the notebook and run `Kernel → Restart & Run All`.
3. The cohort is set by `SEED = 2545448` in Part A (my student ID 002545448; leading zeros dropped because Python 3 rejects them in integer literals). Changing it changes the cohort and the numbers below.

## Key results
| Feature set | CV AUROC |
|---|---|
| Everything | 0.888 ± 0.051 |
| Declared leakage removed | 0.753 ± 0.053 |
| All leakage removed (+ `renal_replacement`) | 0.675 ± 0.064 |

- **Wildcards:** `renal_replacement` = **leakage** (dialysis for organ failure, documented only in the final note, ranked 3rd in importance, mortality 50.7% vs 13.2%); `lactate_not_clearing` = **signal** (admission-time hypoperfusion with a modest, plausible effect: mortality 24.2% vs 19.0%).
- **Naive vs assertion-aware** (`metastatic_disease`): 176 vs 96 — naive matching counts negated mentions like "No evidence of metastatic disease."

## Additional Robustness and Model Comparison Analysis

Three additional analyses were conducted to further assess the reliability of the results. First, the assertion engine was tested using PDAC-specific edge cases, and all cases were classified correctly. Second, a gradient boosting model was evaluated as an alternative to logistic regression, producing scores of 0.852, 0.701, and 0.601, respectively, following the same overall performance trend as the logistic regression model. Finally, a TF-IDF model was trained using the raw clinical notes. Although it achieved a high score of 0.993 when all notes were included, the most influential terms were associated with "required renal replacement therapy," suggesting potential information leakage. When restricted to the admission note alone, the score dropped to 0.531, highlighting the impact of information leakage on model performance.


## Interpretation notes
- About 55% of the model's apparent skill above chance was leakage. The all-features AUROC of 0.89 was driven largely by `comfort_care`, which in PDAC almost announces the outcome.
- The hidden leak (`renal_replacement`) alone cost 0.078 AUROC, more than half as much as my three declared leaks together (0.136). Checking only the leaks you planned for is not enough.
- The honest model (0.675) still has some real signal, mainly age, cachexia and lactate.
- Using only the admission note removes the leakage from the text model entirely, which shows the real fix is a prediction-time cutoff.
- Limitations: the synthetic cohort has perfectly clean note timing (signals only at admission, consequences only at discharge), a fixed three notes per patient (no documentation bias), and only 400 patients (81 deaths; CV SD ≈ ±0.05–0.06). Real charts would need note timestamps and a defined prediction time.

