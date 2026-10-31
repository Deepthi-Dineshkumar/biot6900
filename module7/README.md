# BIOT6900 Module 7: Knowledge Graph for Pancreatic Ductal Adenocarcinoma (PDAC)

**Name:** Deepthi Dineshkumar

**Disease:** Pancreatic ductal adenocarcinoma (PDAC), CPTAC-PDAC cohort via LinkedOmicsKB (continued from Modules 2–6)

**Database:** Neo4j AuraDB Free

**Full write-up:** `Deepthi_Dineshkumar_Module7_Report.pdf`

## Overview

This module turns the separate files built so far into one connected, queryable graph. Genes and their multi-omics evidence come from Module 2, druggability from the Module 3 protein-language-model prior, and variants from the Module 5 pharmacogenomic work. Phenotypes, pathways, drugs and related diseases were curated from public sources, and **every edge carries a `source` and a `confidence`**.

The graph answers the opening question for PDAC (which drugs act on genes associated with it), extends it to shared-symptom repurposing, ranks candidate biomarkers by centrality fused with my own evidence, and groups the graph into communities.

## Where each layer came from

| Layer | Source | Notes |
|---|---|---|
| Genes + evidence | **Module 2** `targets_pdac.csv` | RNA log2FC, protein log2FC, mutation frequency, multi-evidence score, rank |
| Druggability | **Module 3** `targets_pdac_w3.csv` | ESM C tractability priors `p_sm`, `p_ab` |
| Variants | **Module 5** `DPYD.tsv`, `CYP2C19.tsv`, `pgx_classified_CYP2C19.csv` | ClinVar significance and review stars (via MyVariant.info), CADD, gnomAD AF, and my provisional VUS calls |
| Phenotypes | HPO | 6 terms with HPO IDs |
| Pathways | Reactome | 15 pathways |
| Drugs | ChEMBL / DrugBank, FDA labels, NCCN-level PDAC regimens | 15 drugs, with status and indication (status checked against current FDA approvals) |
| Drug metabolism | CPIC guidelines | `METABOLIZED_BY` edges |
| Related diseases | OMIM / HPO / literature | Cystic fibrosis, chronic pancreatitis, MODY8 (CEL-MODY) |

### How the Module 5 variants connect

The Module 5 variants are **pharmacogenomic** (CYP2C19, pooled with DPYD), not PDAC risk variants, so they enter through **drug metabolism** rather than disease association. Neither pharmacogene has an `ASSOCIATED_WITH` edge to PDAC. The two routes are:

- **DPYD → fluorouracil.** DPYD clears fluorouracil, the backbone of FOLFIRINOX and NALIRIFOX, and fluorouracil acts on TYMS, a PDAC drug target. The four CPIC-actionable variants (`rs3918290` \*2A, `rs55886062` \*13, `rs67376798` D949V, `rs75017182` HapB3) were all present in the Module 5 fetch.
- **CYP2C19 → proton pump inhibitors.** PPIs are used alongside pancreatic enzyme replacement for PDAC-related steatorrhea, and CYP2C19 is the main enzyme clearing them. The graph holds \*2, \*3, \*7, \*8, \*17 and one of my provisional "towards pathogenic" VUS.

## Schema

**Node labels:** `Gene`, `Disease`, `Phenotype`, `Pathway`, `Drug`, `Variant` (all also `:Entity` with a unique `name`).

**Relationship types:**

| Relationship | Meaning |
|---|---|
| `(Gene)-[:ASSOCIATED_WITH]->(Disease)` | evidence that the gene relates to the disease |
| `(Disease)-[:HAS_PHENOTYPE]->(Phenotype)` | clinical presentation (HPO) |
| `(Gene)-[:PARTICIPATES_IN]->(Pathway)` | Reactome membership |
| `(Drug)-[:TARGETS]->(Gene)` | mechanism of action |
| `(Variant)-[:OCCURS_IN]->(Gene)` | Module 5 variants |
| `(Drug)-[:METABOLIZED_BY]->(Gene)` | **added** — pharmacokinetics, so Module 5 variants reach therapy |
| `(Drug)-[:TREATS]->(Phenotype)` | **added** — symptom-level treatment |
| `(Gene)-[:PREDICTS_RESPONSE_TO]->(Drug)` | **added** — biomarker-directed therapy (BRCA2 → olaparib) |

**Confidence scale on every edge:** 1.0 curated and definitive (approved drug mechanism, Mendelian gene–disease, CPIC level A), 0.8–0.9 curated annotation (Reactome, HPO, ClinVar ≥2 stars) 0.5–0.7 moderate (investigational drugs, ClinVar conflicting or low-star, non-Mendelian risk genes), 0.2–0.4 weak (pharmacology-only links, adjunct therapy, my provisional calls). Gene–PDAC edges carry my Module 2 multi-evidence score as their confidence, so the graph inherits that calibration rather than inventing a new one.

## Results

**Graph size:** 80 nodes (30 genes, 15 pathways, 15 drugs, 10 variants, 6 phenotypes, 4 diseases) and 101 relationships. Every relationship has a source and a confidence (checked in the notebook).

**Identifier harmonization.** Module 2 keys genes on Ensembl IDs while every other source uses HGNC symbols. Genes are keyed on one symbol with the Ensembl ID kept as a property, and nodes merge rather than overwrite, so a gene arriving from several sources stays one node. The notebook asserts that every gene node has exactly one Ensembl ID and that no two symbols share an ID. Two genes are shared across diseases and were merged, not duplicated:

- **CFTR** → PDAC, cystic fibrosis, chronic pancreatitis
- **CEL** → PDAC, MODY8 (CEL-MODY)

**Opening question (2 hops):** 12 drug–gene–PDAC rows, each showing the evidence and confidence on both edges. **Shared-symptom extension (5 hops):** 4 rows (2 CFTR modulators × 2 related diseases, with the shared symptoms collected per row), all running through CFTR via cystic fibrosis and chronic pancreatitis. Three full paths are printed node by node with the evidence on each hop.

**Top candidate biomarkers** (fused = 0.4 × PageRank percentile + 0.4 × Module 2 evidence + 0.2 × druggability):

| Rank | Gene | Fused | PageRank | Module 2 score (rank) | Druggability |
|---|---|---|---|---|---|
| 1 | KRAS | 0.885 | 0.0164 | 0.811 (69) | 1.000 (approved drug) |
| 2 | XDH | 0.841 | 0.0164 | 0.952 (9) | 0.504 (Module 3 prior) |
| 3 | CFTR | 0.760 | 0.0164 | 0.978 (2) | 0.044 (Module 3 prior) |
| 4 | CEL | 0.668 | 0.0061 | 0.962 (4) | 0.716 (Module 3 prior) |
| 5 | MUC5B | 0.626 | 0.0061 | 0.942 (13) | 0.546 (Module 3 prior) |

The notebook also prints a sensitivity check across four weightings; KRAS leads all four, and XDH, CFTR and CEL stay in the top 5 in all of them.

**Drug-status correction.** Daraxonrasib was first annotated as investigational, reflecting its trial status when the drug layer was curated. The FDA approved it on 26 August 2026 for metastatic pancreatic adenocarcinoma after at least one prior systemic therapy, so the annotation was corrected to approved and the notebook re-run. That raises the KRAS tractability proxy from 0.5 to 1.0 and moves KRAS from second to first (fused 0.885). All committed artifacts reflect the corrected run.

**Communities:** 8, modularity 0.626. They separate cleanly into a PDAC tumour-biology core, a pancreatic-insufficiency cluster (CFTR, CEL, PRSS1, SPINK1, CTRC with cystic fibrosis, chronic pancreatitis and MODY8), a fluoropyrimidine cluster (DPYD variants with fluorouracil, capecitabine, TYMS, RRM1), a CYP2C19/PPI cluster, a KRAS–EGFR signalling cluster, and small DNA-repair, purine-catabolism and TP53/CDKN2A clusters.

## Files in this folder

| File | Description |
|---|---|
| `BIOT6900_Module7_PDAC.ipynb` | Completed notebook; runs top to bottom with Restart & Run All |
| `my_knowledge_graph.png` | Graph coloured by node kind, laid out by community, with a legend |
| `nodes.csv` | All 80 nodes with their properties and community assignment |
| `edges.csv` | All 101 edges with type, source and confidence |
| `biomarker_ranking.csv` | Ranked candidate-biomarker table from §5 |
| `Deepthi_Dineshkumar_Module7_Report.pdf` | Report answering the four Part F questions |
| `data/` | Inputs copied from my earlier modules (see below) |
| `README.md` | This file |

`data/` holds `targets_pdac.csv` and `ensembl_symbol_mapping.csv` (Module 2/3), `targets_pdac_w3.csv` (Module 3), and `DPYD.tsv`, `CYP2C19.tsv`, `pgx_classified_CYP2C19.csv` (Module 5). They are copies, so this folder runs on its own.

## How to reproduce

1. Create a free Neo4j AuraDB instance at `console.neo4j.io` and download the credentials file. On AuraDB Free the **username is the instance ID**, not `neo4j`.
2. Save that file in this folder as `neo4j_credentials.env`. It is git-ignored and is **not** committed; the notebook reads it rather than holding any password in a cell.
3. Install the packages:

        pip install neo4j networkx matplotlib pandas requests

4. Open `BIOT6900_Module7_PDAC.ipynb` and choose Kernel, then Restart & Run All.

**Expected results:** "Connected to Neo4j"; 80 nodes and 101 relationships loaded, matching the Python build; 0 relationships missing source or confidence; 12 rows from the 2-hop query and 4 from the 5-hop query; 3 full paths printed; KRAS first in the biomarker table, XDH second; 8 communities with modularity 0.626; and `my_knowledge_graph.png`, `nodes.csv`, `edges.csv` and `biomarker_ranking.csv` written to this folder.

`RESET = True` in §2 clears the database before loading, so a re-run always reproduces exactly this graph instead of adding to an earlier one.

## Things that did not work, and choices made

- **Patient layer skipped.** Real PDAC patient genotype or outcome data were unavailable, and inventing a synthetic cohort would add nodes without adding evidence, so `PATIENTS` is left empty and the stratification query is not run.
- **Pharmacology-only genes carry weak edges.** TYMS, RRM1, TOP1, EGFR and PARP1 are not in my Module 2 table, so their link to PDAC rests only on a PDAC drug acting on them. Those edges are given confidence 0.3 and labelled as pharmacology-only, which is also why they rank mid-table rather than at the top.
- **CYP2C19 appears in my Module 2 table** (rank 493) but is deliberately not linked to PDAC, because differential expression in tumour tissue is not evidence that it drives the disease. It is in the graph as a pharmacogene only.
- **One provisional VUS edge is deliberately weak.** `rs371241592` enters at confidence 0.3, because my Module 5 call came from pooled, DPYD-dominated CADD calibration and is not a clinical reclassification.
- **Drug status is a dated claim.** Daraxonrasib moved from investigational to approved during this project (FDA, 26 August 2026), which changed the top of the biomarker ranking. Approval status and PDAC indication would be better stored as separate, dated annotations than collapsed into one status field.

