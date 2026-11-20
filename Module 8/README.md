# Module 8 — Single-Cell & Spatial Omics (Breast Cancer)

**Name:** Deepthi Dineshkumar · **Target gene:** ERBB2 (HER2) · **Course:** BIOT 6900, Module 8 Computational Lab

## Files

| File | What it is |
| --- | --- |
| `BIOT6900_Module8_Assignment_Starter.ipynb` | Completed notebook: every TODO filled, all cells run top to bottom with outputs saved |
| `dineshkumar_deepthi_module8_report.pdf` | The 4–5 page report linking which cells express ERBB2 to where they sit |
| `figures/` | All report figures (PNG), score tables (CSV) and `results_summary.json` with every key number |
| `README.md` | This file |

## Key results

ERBB2 is expressed almost only by tumor cells, and in the Visium section it sits evenly through the tumor core and margin.

- **Data:** 20,000 cells from 10 donors (CELLxGENE Census 2025-11-08; 10x 3′ v3, 3′ v2, 5′ v2), 18,388 after QC. Spatial: 10x Visium breast cancer Block A Section 1, 3,798 spots.
- **Integration:** scVI chosen (batch mixing 0.968, bio conservation 0.516, overall 0.697) over Harmony (0.912 / 0.495 / 0.662) and no integration (0.872 / 0.489 / 0.642). The scVI–Harmony gap is small (0.035).
- **Which cells:** ERBB2 is in 36% of 8,678 tumor cells vs 2.5–3.8% of T cells, macrophages and endothelial cells. It ranges from 93% (HTAPP-313) to 3% (pal_Patient 0178) of tumor cells across donors.
- **Where:** 2,376 tumor, 727 immune and 695 stroma spots. ERBB2 is high in tumor core (mean 1.48, 99% positive) and margin (1.44, 97%), with no difference (p = 0.12). Tumor–immune neighborhood z = −56.2 (immune-excluded-like); immune–stroma z = +8.2.

## What did not run, and what I changed

Every step ran to completion. These are the deviations from the starter, and why:

- **GPU step ran on CPU:** the assigned courses-gpu node had a Tesla P100 (CUDA capability 6.0), which the installed PyTorch 2.8 did not support ("no kernel image is available"). I hid the GPU (`CUDA_VISIBLE_DEVICES=""`) and trained scVI on CPU for 60 epochs instead of 200.
- **Cell-type label bug fixed:** the starter's keyword mapping matched "t cell" inside "malignant cell", "myofibroblast cell" and "mast cell", putting 8,526 malignant cells in the T group. A correction cell after 1.1b reassigns them, confirmed by marker genes. The CD68 suggestion printed by 1.1b came before this fix, so ERBB2 was kept.
- **Census access:** opened with `census_version="2025-11-08"` and Explorer's HTTP proxy passed to TileDB (`vfs.s3.proxy_*`), because compute nodes reach S3 only through the proxy.
- **Environment:** packages were installed in a virtual environment on the node's local `/tmp` disk, because the shared home storage was too slow. NumPy-2 incompatibilities in h5py, bottleneck, numexpr and transformers were fixed by upgrading them.
- **Added analysis:** section 1.8b splits ERBB2 by donor to show patient-to-patient differences.

## How to re-run

1. On Explorer, open Open OnDemand → JupyterLab (courses-gpu, 4 CPUs, 32 GB, about 6 hours).
2. In a terminal, build the environment on the node's local disk and register the kernel:

   ```bash
   source ~/m8_setup.sh   # creates /tmp/$USER-m8, installs packages, registers "Python 3 (Module 8)"
   ```

3. Open the notebook, choose **Kernel → Change Kernel → Python 3 (Module 8)**, then **Kernel → Restart Kernel and Run All Cells**. The full run takes about 30–40 minutes on CPU, including the Census download.
4. Figures and `results_summary.json` are written to `figures/`.
