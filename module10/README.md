# Module 10 — Medical Imaging AI: Build It, Then Break It (Computational Lab 10)

**BIOT 6900 · AI-Powered Biotech Workflows · Deepthi Dineshkumar**

**Question:** a pretrained DenseNet121 gets high accuracy on tumour-vs-normal colon histology easily. Can I *prove* whether to believe it? I fine-tune a clean model (Part A), probe it with Grad-CAM and occlusion (Part B), design my own confounder, a **scanner / slide-label barcode tag** on tumour tiles (Part C), and test whether the retrained model learned the shortcut, using the clean model as a control (Part D). The write-up is `dineshkumar_deepthi_module10_report.pdf` (Part E).

**Verdict: the confounded model learned "barcode → tumour".** The clean control largely ignored the tag (11.6% flip rate) while the confounded model flipped 100% of tagged normals and lost 97% of its tumour recall when the tag was removed.

## Files
| File | What it is |
|---|---|
| `BIOT6900_Module10_Graded_Confounder.ipynb` | Graded notebook, Parts A–E, executed top to bottom (Restart & Run All) |
| `dineshkumar_deepthi_module10_report.pdf` | 4-page report (Part E) |
| `figures/` | Every figure used in the notebook and report (`00_tiles`, `A_*`, `B1–B3`, `C1–C2`, `D1–D3`) |
| `results/module10_results.json` | Every number reported in the notebook and report |
| `data/` | **Not committed**, rebuilt automatically, see below |

## How to reproduce
1. **Environment:** the course `biot6900` conda env, selected as the notebook kernel.
   `pip install "monai[pillow]" "medmnist>=3.0" scikit-learn matplotlib` and `pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu`
   (Versions used: torch 2.15.0+cpu, MONAI 1.6.2.)
2. **Data:** the notebook looks for `data/pathmnist_tumour_vs_normal_subset.npz`. The course subset file was **not available on Canvas**, and the 64-px PathMNIST download repeatedly timed out, so this run **rebuilt the subset from the 28-px PathMNIST release** (`DL_SIZE = 28`, ~0.2 GB, downloaded into `data/medmnist/`). If the subset file is missing, the data cell does this automatically: 300 tiles per class for training, 200 per class for validation, 300 per class for test, sampled with seed 2545448. If you have the Canvas file, drop it into `data/` and the download is skipped (numbers will then differ slightly, because those tiles are 64-px-derived).
3. **Run:** open the notebook, select the `biot6900` kernel, then `Kernel → Restart & Run All`. On a laptop CPU it takes **about 4 minutes**, plus a one-time ~30 MB download of ImageNet DenseNet121 weights.
4. **Seed:** `2545448` (my student ID without leading zeros) is fixed. Small CPU/library differences in the last decimal are expected.

## Pipeline
- **Data:** PathMNIST (from NCT-CRC-HE-100K). Class 6, normal colon mucosa → label 0; class 8, colorectal adenocarcinoma epithelium → label 1. 600 train / 400 val / 600 test, balanced, resized to 128 px, ImageNet-normalised.
- **Model:** MONAI DenseNet121, ImageNet-pretrained, fine-tuned with Adam (lr 1e-4) for 3 epochs, batch 32.
- **Explainability:** Grad-CAM on `class_layers.relu`, with a min-max fix so high = important. Hand-built occlusion uses a 16-px grey patch at stride 8. Agreement and location are *measured* (evidence share on tissue or on the tag; Pearson r between maps), not eyeballed.
- **Confounder:** a fixed off-white barcode label, 44×18 px (4.8% of the tile), in the top-left corner, drawn on every tumour tile in training and validation. A fallback `CONFOUNDER = "stain_cast"` exists but was not needed.
- **Battery (run on both models):**
  - tumour recall with vs without the tag;
  - flip rate of correctly classified normals when the tag is added;
  - Grad-CAM and occlusion share on the tag, averaged over 8 tiles.

## Key results
| | Clean model (control) | Confounded model |
|---|---|---|
| Validation accuracy | 0.940 | **1.000** (tagged validation set) |
| Test accuracy, untagged tiles | 0.932 (sens 0.947, spec 0.917) | **0.512** (detects 8/300 tumours) |
| Tumour recall with tag → without tag | 0.987 → 0.947 | 1.000 → **0.027** |
| Normals called tumour: no tag → tag added | 0.083 → 0.190 | 0.003 → **1.000** |
| Flip rate normal → tumour | 0.116 | **1.000** (gap 0.884) |
| Mean Grad-CAM share on tag (tag = 4.8% of tile) | 16% | 15% |
| Mean occlusion share on tag | 23% | **66%** |

Iteration log: one attempt (scanner tag 44×18 px) was enough. The clean control's flip rate stayed below the 30% "inconclusive" threshold, so no re-design was needed.

## Interpretation notes
- **Clean model (Part B):** Grad-CAM and occlusion both keep their evidence on tissue (80–98%, matching the tissue area), and the maps change with the target class. The two methods agree only partly on *where* within the tissue (r = −0.03 to 0.62), partly because Grad-CAM uses a coarse 4×4 map. So the model looks tissue-driven, but I can't say which structures it relies on. It has no detectable background or border shortcut.
- **The control is the proof:** both models saw identical tagged tiles. The clean model barely reacted, so the tag is not inherently confusing. The gap between the two models (0.88 in flip rate, 0.93 in recall drop) comes from what the confounded model learned in training.
- **The clean control was not perfectly indifferent:** tagged normals called tumour rose from 8.3% to 19.0%. An unfamiliar high-contrast pattern nudges any model, and the control measures that baseline so it isn't mistaken for learning.
- **The confounded model abandoned morphology:** without the tag it calls almost nothing tumour. A constant, perfectly predictive pattern is easier to fit than tissue (shortcut learning). This 4.8% tag was learned *more* completely than Notebook 2's larger ink mark (recall without ink 0.19).
- **Grad-CAM alone would have missed it:** Grad-CAM's tag share was the same for both models (15% vs 16%). Occlusion (66% vs 23%) and the behavioural tests separated them. This is why the lab uses two methods and measures them.
- **Why accuracy misleads:** the worst model had the best validation score (1.000 vs 0.940), because the confound was in validation too. Deployed at a site without the tag, it would miss 97% of cancers.
- **Limitations:** the artifact is synthetic and perfectly correlated with the label. The data are 28-px tiles, upsampled, so Grad-CAM is coarse. There is one seed, 3 epochs and 600 training tiles. The thresholds (tissue mask < 0.80, "points at it" > 2× area, 30% cut-offs) are my own choices.
