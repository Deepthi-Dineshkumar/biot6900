"""
BIOT6900 Portfolio Deliverable 1: PDAC

No primary lab analyses are rerun here. The script reads outputs already produced in
Modules 2-4 and summarises them for the Deliverable 1 report. The only additions are a
Benjamini-Hochberg FDR adjustment of the Module 2 p values (across every gene tested in
each layer) and explicitly stated evidence thresholds, used for Figure 1c. Outputs:

Run from anywhere:  python Portfolio/pd1_figures.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------------------
# Paths (relative to the repo root, so the script works on any machine)
# --------------------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
REPO = HERE.parent
M2 = REPO / "module2"
M3 = REPO / "Module 3"
M4 = REPO / "module4"
STRUCT = M4 / "structures"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

INPUTS = {
    "Module 2 ranked table": M2 / "targets_pdac.csv",
    "Module 2 RNA layer": M2 / "data" / "pdac_transcriptomics.tsv",
    "Module 2 protein layer": M2 / "data" / "pdac_proteomics.tsv",
    "Module 2 mutation layer": M2 / "data" / "pdac_mutation.tsv",
    "Module 3 re-ranked table": M3 / "targets_pdac_w3.csv",
    "Module 4 structure summary": M4 / "structure_assessment_pdac.csv",
}
missing = [f"{k}: {v}" for k, v in INPUTS.items() if not v.exists()]
if missing:
    raise FileNotFoundError("Missing inputs (see README, 'How to reproduce'):\n  " + "\n  ".join(missing))

# Evidence thresholds for the second Venn diagram (Figure 1c)
RNA_LFC, PROT_LFC, FDR, MUT_FREQ = 1.0, 0.5, 0.05, 0.05

# Results carried over from Module 4 (printed in BIOT6900_Module4_Assignment.ipynb)
STRUCTURES = {
    "CEL": {"acc": "P19835", "pockets": {56: "#D62728"},
            "trusted": (21, 552), "low": [(3, 17), (555, 592), (717, 753)]},
    "XDH": {"acc": "P47989", "pockets": {2: "#D62728", 1: "#7B4EA3"},
            "trusted": (2, 1333), "low": [(168, 190)]},
}
# SAL contact residues in 2E1Q chain A, UniProt numbering (Module 4 ligand-site check)
XDH_SAL_SITE = [803, 874, 877, 881, 915, 1009, 1010, 1011, 1012, 1015, 1079, 1080]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})


# --------------------------------------------------------------------------------------
# Load Module 2-4 outputs
# --------------------------------------------------------------------------------------
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


m2 = pd.read_csv(INPUTS["Module 2 ranked table"])
m2 = m2.sort_values("score", ascending=False).reset_index(drop=True)
m2["m2_rank"] = np.arange(1, len(m2) + 1)
# BH-FDR across every gene tested in each layer (the full source tables), not only the 1,880
layer_tables = {name: pd.read_csv(INPUTS[f"Module 2 {name} layer"], sep="\t")
                for name in ("RNA", "protein", "mutation")}
for name, col in (("RNA", "rna_q"), ("protein", "prot_q")):
    tab = layer_tables[name]
    m2[col] = m2.gene.map(dict(zip(tab.gene, bh_fdr(tab.pval.fillna(1)))))

m3 = pd.read_csv(INPUTS["Module 3 re-ranked table"]).rename(columns={"rank": "m3_rank"})
m4 = pd.read_csv(INPUTS["Module 4 structure summary"])
symbol = dict(zip(m3.original_gene_id, m3.gene))   # Ensembl -> symbol for the triaged 15
LABELS = {"ENSG00000133703": "KRAS", "ENSG00000141510": "TP53", "ENSG00000147889": "CDKN2A",
          "ENSG00000141646": "SMAD4", **symbol}

layers = {name: set(tab.gene) for name, tab in layer_tables.items()}


# --------------------------------------------------------------------------------------
# Figure 1: integration scatter + Venn diagrams
# --------------------------------------------------------------------------------------
def venn_regions(a, b, c):
    return {"100": len(a - b - c), "010": len(b - a - c), "001": len(c - a - b),
            "110": len(a & b - c), "101": len(a & c - b), "011": len(b & c - a), "111": len(a & b & c)}


def draw_venn(ax, sets, names, colours, title):
    """Unweighted three-set Venn (layer sizes differ by >10x, so areas are not to scale)."""
    r, centres = 1.0, [(-0.58, 0.35), (0.58, 0.35), (0.0, -0.62)]
    for (x, y), col in zip(centres, colours):
        ax.add_patch(Circle((x, y), r, fc=col, ec=col, alpha=0.22, lw=1.2))
    reg = venn_regions(*sets)
    pos = {"100": (-1.0, 0.6), "010": (1.0, 0.6), "001": (0, -1.15), "110": (0, 0.78),
           "101": (-0.62, -0.36), "011": (0.62, -0.36), "111": (0, 0.0)}
    for k, (x, y) in pos.items():
        ax.text(x, y, f"{reg[k]:,}", ha="center", va="center",
                fontsize=8.5 if k == "111" else 7.5, fontweight="bold" if k == "111" else "normal")
    label_pos = [(-1.45, 1.45, "right"), (1.45, 1.45, "left"), (0.0, -1.9, "center")]
    for (x, y, ha), name, s, col in zip(label_pos, names, sets, colours):
        ax.text(x, y, f"{name}\n(n={len(s):,})", ha=ha, va="center", fontsize=6.8, color=col,
                fontweight="bold")
    ax.set_xlim(-3.2, 3.2); ax.set_ylim(-2.25, 1.85); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, loc="left")
    return reg


fig = plt.figure(figsize=(7.4, 3.9))
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], hspace=0.3, wspace=0.05)
ax = fig.add_subplot(gs[:, 0])
conc = m2.concordant.astype(bool)
size = 6 + 260 * m2.mut_freq
ax.scatter(m2.rna_lfc[conc], m2.prot_lfc[conc], s=size[conc], c="#4C78A8", alpha=0.35, lw=0,
           label=f"concordant ({conc.sum():,})")
ax.scatter(m2.rna_lfc[~conc], m2.prot_lfc[~conc], s=size[~conc], c="#F58518", alpha=0.45, lw=0,
           label=f"discordant ({(~conc).sum():,})")
top = m2.m2_rank <= 15
ax.scatter(m2.rna_lfc[top], m2.prot_lfc[top], s=size[top], facecolors="none", edgecolors="black", lw=0.8,
           label="Module 2 top 15")
ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
offsets = {"CEL": (8, -2), "XDH": (18, -4), "MUC16": (8, 12), "CFTR": (18, -6), "COL11A1": (6, 10),
           "KRT19": (-34, 8), "LRP1B": (10, -6), "KRAS": (40, -30), "TP53": (-24, 16), "CDKN2A": (34, -14),
           "SMAD4": (-34, -12), "MUC5B": (52, -2), "ACADL": (-36, -6)}
for gid, lab in LABELS.items():
    if lab not in offsets:
        continue
    row = m2.loc[m2.gene == gid].iloc[0]
    col = "#B5442B" if lab in ("KRAS", "TP53", "CDKN2A", "SMAD4") else "black"
    ax.annotate(lab, (row.rna_lfc, row.prot_lfc), xytext=offsets[lab], textcoords="offset points",
                fontsize=6.5, color=col, fontweight="bold" if lab in ("CEL", "XDH") else "normal",
                arrowprops=dict(arrowstyle="-", lw=0.4, color="grey"))
for f in (0.02, 0.10, 0.50):
    ax.scatter([], [], s=6 + 260 * f, c="lightgrey", label=f"mutated in {f:.0%} of tumours")
ax.legend(fontsize=5.8, loc="lower right", frameon=False, handletextpad=0.3, labelspacing=0.5)
ax.set_xlabel("RNA log$_2$FC (tumour vs normal)"); ax.set_ylabel("Protein log$_2$FC (tumour vs normal)")
ax.set_title(f"a  {len(m2):,} genes surviving the three-layer join", loc="left")

cov = draw_venn(fig.add_subplot(gs[0, 1]), [layers["RNA"], layers["protein"], layers["mutation"]],
                ["RNA", "Protein", "Mutation"], ["#4C78A8", "#54A24B", "#B5442B"],
                "b  Layer coverage (gene IDs)")
ev_sets = [set(m2.gene[(m2.rna_q < FDR) & (m2.rna_lfc.abs() >= RNA_LFC)]),
           set(m2.gene[(m2.prot_q < FDR) & (m2.prot_lfc.abs() >= PROT_LFC)]),
           set(m2.gene[m2.mut_freq >= MUT_FREQ])]
ev = draw_venn(fig.add_subplot(gs[1, 1]), ev_sets, ["RNA evidence", "Protein evidence", "Recurrent mutation"],
               ["#4C78A8", "#54A24B", "#B5442B"], "c  Evidence overlap within the 1,880")
fig.savefig(FIG / "fig1_multiomics_integration.png", bbox_inches="tight")
plt.close(fig)

all_three = [LABELS.get(g, g) for g in ev_sets[0] & ev_sets[1] & ev_sets[2]]
rp = ev_sets[0] & ev_sets[1]
rp_conc = int(m2.set_index("gene").loc[list(rp), "concordant"].sum())
pd.DataFrame([{"diagram": "coverage", **cov}, {"diagram": "evidence", **ev}]).to_csv(HERE / "venn_counts.csv", index=False)
print(f"Fig 1: {len(m2)} genes, {conc.sum()} concordant; coverage centre {cov['111']}; "
      f"evidence RNA&protein {len(rp)} ({rp_conc} concordant); all three: {all_three}")


# --------------------------------------------------------------------------------------
# Figure 2: structure models with kept pockets
# --------------------------------------------------------------------------------------
def read_atoms(path):
    rows = []
    for line in open(path):
        if line.startswith("ATOM"):
            rows.append((line[12:16].strip(), int(line[22:26]), float(line[30:38]),
                         float(line[38:46]), float(line[46:54]), float(line[60:66] or 0)))
    return pd.DataFrame(rows, columns=["atom", "resnum", "x", "y", "z", "b"])


def load_model(acc):
    """Prefer the AlphaFold download (pLDDT in the B-factor column); fall back to fpocket's copy."""
    af = STRUCT / f"AF-{acc}.pdb"
    if af.exists():
        return read_atoms(af), True
    out = STRUCT / f"AF-{acc}_out" / f"AF-{acc}_out.pdb"   # fpocket overwrites B-factors
    if out.exists():
        return read_atoms(out), False
    return None, False


def view_rotation(coords, toward):
    """Rotation that points `toward` at the viewer (+z) and spreads the rest along x."""
    z = toward / np.linalg.norm(toward)
    proj = coords - np.outer(coords @ z, z)
    _, _, vt = np.linalg.svd(proj - proj.mean(0), full_matrices=False)
    x = vt[0] - (vt[0] @ z) * z; x /= np.linalg.norm(x)
    return np.vstack([x, np.cross(z, x), z]).T


AF_BANDS = [(90, "#0053D6", "pLDDT > 90"), (70, "#65CBF3", "70-90"), (50, "#FFDB13", "50-70"),
            (-1, "#FF7D45", "< 50")]


def residue_colours(ca, info, has_plddt):
    if has_plddt:
        return [next(c for cut, c, _ in AF_BANDS if v > cut) for v in ca.b]
    lo, hi = info["trusted"]
    cols = []
    for r in ca.resnum:
        if any(a <= r <= b for a, b in info["low"]):
            cols.append("#FF7D45")
        else:
            cols.append("#0053D6" if lo <= r <= hi else "#B0B0B0")
    return cols


fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))
plddt_mode = None
for ax, (gene, info) in zip(axes, STRUCTURES.items()):
    atoms, has_plddt = load_model(info["acc"])
    pocket_files = [STRUCT / f"AF-{info['acc']}_out" / "pockets" / f"pocket{pid}_atm.pdb" for pid in info["pockets"]]
    if atoms is None or not all(f.exists() for f in pocket_files):
        ax.text(0.5, 0.5, f"{gene}: model or fpocket output not found in\n{STRUCT}\nRun the Module 4 notebook first.",
                ha="center", va="center", transform=ax.transAxes); ax.axis("off"); continue
    plddt_mode = has_plddt
    ca = atoms[atoms.atom == "CA"].reset_index(drop=True)
    xyz = ca[["x", "y", "z"]].to_numpy()
    centre = xyz.mean(0)
    pockets = {}
    for pid in info["pockets"]:
        p = read_atoms(STRUCT / f"AF-{info['acc']}_out" / "pockets" / f"pocket{pid}_atm.pdb")
        pockets[pid] = p
    main = pockets[list(info["pockets"])[0]]
    R = view_rotation(xyz - centre, main[["x", "y", "z"]].to_numpy().mean(0) - centre)
    v = (xyz - centre) @ R
    cols = residue_colours(ca, info, has_plddt)
    zmin, zmax = v[:, 2].min(), v[:, 2].max()
    depth = lambda zz: (zz - zmin) / (zmax - zmin)
    segs = [(i, (v[i, 2] + v[i + 1, 2]) / 2) for i in range(len(v) - 1)
            if np.linalg.norm(xyz[i + 1] - xyz[i]) < 4.5]
    for i, zz in sorted(segs, key=lambda s: s[1]):
        d = depth(zz)
        ax.plot(v[i:i + 2, 0], v[i:i + 2, 1], color=cols[i], lw=0.5 + 1.6 * d, alpha=0.25 + 0.75 * d,
                solid_capstyle="round", zorder=1 + d)
    for pid, colour in info["pockets"].items():
        pv = (pockets[pid][["x", "y", "z"]].to_numpy() - centre) @ R
        ax.scatter(pv[:, 0], pv[:, 1], s=14 + 22 * depth(pv[:, 2]), c=colour, edgecolors="white",
                   lw=0.3, alpha=0.95, zorder=5)
    if gene == "XDH":
        site = ca.resnum.isin(XDH_SAL_SITE).to_numpy()
        ax.scatter(v[site, 0], v[site, 1], s=34, facecolors="none", edgecolors="black", lw=0.9, zorder=6)
    ax.set_aspect("equal"); ax.axis("off")
    row = m4.set_index("gene").loc[gene]
    ax.set_title(f"{'a' if gene == 'CEL' else 'b'}  {gene} ({info['acc']}), mean pLDDT {row.mean_pLDDT}",
                 loc="left")

pocket_meta = {"CEL": "CEL pocket 56: score 0.925, lining pLDDT 93.7",
               "XDH2": "XDH pocket 2: score 0.979, pLDDT 98.2 (SAL site)",
               "XDH1": "XDH pocket 1: score 0.933, pLDDT 96.6 (FAD site)"}
if plddt_mode:
    conf = [Patch(color=c, label=l) for _, c, l in AF_BANDS]
else:
    conf = [Patch(color="#0053D6", label="PAE-confident domain"),
            Patch(color="#FF7D45", label="pLDDT < 50 stretch"),
            Patch(color="#B0B0B0", label="outside PAE-confident domain")]
marks = [Line2D([], [], marker="o", ls="", color="#D62728", label=pocket_meta["CEL"]),
         Line2D([], [], marker="o", ls="", color="#D62728", label=pocket_meta["XDH2"]),
         Line2D([], [], marker="o", ls="", color="#7B4EA3", label=pocket_meta["XDH1"]),
         Line2D([], [], marker="o", ls="", mfc="none", mec="black", label="2E1Q SAL contact residues")]
fig.legend(handles=conf + marks, loc="lower center", ncol=3, fontsize=6.3, frameon=False,
           bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(FIG / "fig2_structures.png", bbox_inches="tight")
plt.close(fig)
print("Fig 2: " + ("models not found, panels show a note" if plddt_mode is None else
       f"structures coloured by {'per-residue pLDDT' if plddt_mode else 'Module 4 confidence regions'}"))


# --------------------------------------------------------------------------------------
# Figure 3 + CSV: prioritization across all three analyses
# --------------------------------------------------------------------------------------
# Decision criteria (see report Section 5):
#   C1 direction  - increased in tumour at both RNA and protein. This supports investigating an
#                   inhibitory approach; it does not show that the tumour depends on the gene.
#   C2 tractable  - best modality prior >= 0.4, or prior clinical investigation for this gene
#   C3 structure  - best pocket has Module 4 verdict "keep" (lining pLDDT >= 70, internal PAE <= 10 A)
#                   AND fpocket druggability score >= 0.5. Only CEL and XDH were assessed.
POCKET_SCORE = 0.5
PRECEDENCE = {"XDH": "allopurinol, febuxostat", "MUC16": "phase I anti-MUC16 ADC incl. PDAC (investigated, benefit unproven)",
              "CFTR": "ivacaftor (potentiator)", "FN1": "EDB-domain antibody fusions"}
CALLS = {   # analyst call after applying C1-C3; rationale is printed into the CSV
    "XDH": ("Advance", 1, "Small molecule: up in tumour, enzyme-like prior, top pocket covers 92% of the 2E1Q SAL site"),
    "MUC16": ("Advance", 2, "Antibody/ADC: Module 2 #1; only gene meeting all three selected evidence thresholds "
                             "in the 1,880-gene intersection; prior clinical investigation (phase I ADC incl. PDAC), "
                             "not proven PDAC benefit; structure pending"),
    "CEL": ("Conditional", 3, "Strong kept pocket, but decreased in bulk tumour tissue; cellular origin and "
                              "therapeutic direction remain unresolved"),
    "MUC5B": ("Watch", None, "Antibody prior 0.55 but secreted gel-forming mucin: soluble-antigen sink"),
    "FN1": ("Watch", None, "Prior 0.48 but precedence is isoform (EDB) specific; gene-level data cannot resolve"),
}
tri = m3.merge(m2[["gene", "m2_rank", "rna_lfc", "prot_lfc", "mut_freq", "rna_q", "prot_q"]],
               left_on="original_gene_id", right_on="gene", suffixes=("", "_ens")).drop(columns="gene_ens")
tri = tri.merge(m4[[c for c in ("gene", "best_pocket_score", "best_pocket_verdict", "druggability_call")
                    if c in m4.columns]], on="gene", how="left")
tri["direction"] = np.where((tri.rna_lfc > 0) & (tri.prot_lfc > 0), "up",
                            np.where((tri.rna_lfc < 0) & (tri.prot_lfc < 0), "down", "mixed"))
tri["best_modality"] = np.where(tri.p_sm >= tri.p_ab, "small molecule", "antibody")
tri["clinical_precedence"] = tri.gene.map(PRECEDENCE).fillna("")
tri["C1_direction"] = tri.direction.eq("up")
tri["C2_tractable"] = (tri[["p_sm", "p_ab"]].max(axis=1) >= 0.4) | tri.clinical_precedence.ne("")
if "best_pocket_verdict" in tri.columns:
    kept = tri.best_pocket_verdict.astype(str).str.lower().eq("keep")
else:   # Module 4 exported only the best KEPT pocket, so its score already passed the confidence filter
    print("Note: no best_pocket_verdict column; assuming best_pocket_score is the best kept pocket")
    kept = tri.best_pocket_score.notna()
tri["C3_structure"] = (kept & tri.best_pocket_score.ge(POCKET_SCORE)).where(tri.druggability_call.notna(), other=np.nan)
tri["criteria_met"] = tri[["C1_direction", "C2_tractable"]].sum(axis=1) + tri.C3_structure.fillna(0).astype(int)


def default_call(r):
    if r.direction == "down":
        return "Down in tumour", "Decreased in bulk tumour tissue; cellular origin and therapeutic direction remain unresolved"
    return "Low tractability", "Up in tumour but no modality prior >= 0.4 and no precedence"


tri["package_call"], tri["shortlist_rank"], tri["rationale"] = zip(*[
    CALLS.get(r.gene, (*default_call(r)[:1], None, default_call(r)[1])) for r in tri.itertuples()])
tri["shortlist_rank"] = pd.to_numeric(tri.shortlist_rank)
order = {"Advance": 0, "Conditional": 1, "Watch": 2, "Low tractability": 3, "Down in tumour": 4}
tri = tri.sort_values(["shortlist_rank", "package_call", "m3_rank"],
                      key=lambda s: s.map(order) if s.name == "package_call" else s, na_position="last")
cols = ["shortlist_rank", "gene", "m2_rank", "m3_rank", "rna_lfc", "prot_lfc", "mut_freq", "direction",
        "p_sm", "p_ab", "best_modality", "clinical_precedence", "best_pocket_score", "druggability_call",
        "C1_direction", "C2_tractable", "C3_structure", "criteria_met", "package_call", "rationale"]
tri["shortlist_rank"] = tri.shortlist_rank.astype("Int64")
tri[cols].round(3).to_csv(HERE / "prioritized_targets_pdac.csv", index=False)

# --- matrix figure ---
CALL_COL = {"Advance": "#2E7D32", "Conditional": "#E69F00", "Watch": "#8DA0CB",
            "Low tractability": "#BDBDBD", "Down in tumour": "#D9D9D9"}
heads = ["Module 2\nrank", "RNA\nlog$_2$FC", "Protein\nlog$_2$FC", "Mutated\n(%)", "Module 3\nrank",
         "p$_{sm}$", "p$_{ab}$", "Structure\n(Module 4)", "Package call"]
n = len(tri)
widths = [0.85, 0.9, 0.9, 0.9, 0.9, 0.7, 0.7, 1.15, 1.75]
x0 = np.concatenate([[1.0], 1.0 + np.cumsum(widths)])
fig, ax = plt.subplots(figsize=(7.4, 0.23 * n + 0.8))
ax.set_xlim(0, x0[-1] + 0.05); ax.set_ylim(n - 0.4, -1.55); ax.axis("off")
cmap_lfc = plt.get_cmap("RdBu_r"); cmap_p = plt.get_cmap("Greens")
for j, h in enumerate(heads):
    ax.text(x0[j] + widths[j] / 2, -0.95, h, ha="center", va="center", fontsize=6.6, fontweight="bold")


def text_colour(fc):
    r, g, b = matplotlib.colors.to_rgb(fc)
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else "black"


for i, r in enumerate(tri.itertuples()):
    short = pd.notna(r.shortlist_rank)
    ax.text(0.95, i, f"{int(r.shortlist_rank)}. {r.gene}" if short else r.gene, ha="right", va="center",
            fontsize=7, fontweight="bold" if short else "normal")
    cells = [
        (f"{r.m2_rank}", "white"),
        (f"{r.rna_lfc:+.2f}", cmap_lfc(0.5 + np.clip(r.rna_lfc, -4, 4) / 8)),
        (f"{r.prot_lfc:+.2f}", cmap_lfc(0.5 + np.clip(r.prot_lfc, -2.5, 2.5) / 5)),
        (f"{100 * r.mut_freq:.1f}", plt.get_cmap("Oranges")(min(r.mut_freq / 0.1, 1) * 0.7)),
        (f"{r.m3_rank}", "white"),
        (f"{r.p_sm:.2f}", cmap_p(r.p_sm * 0.8)),
        (f"{r.p_ab:.2f}", cmap_p(r.p_ab * 0.8)),
        ("not assessed" if pd.isna(r.best_pocket_score) else f"pocket {r.best_pocket_score:.3f}",
         "#F2F2F2" if pd.isna(r.best_pocket_score) else ("#A5D6A7" if r.gene == "XDH" else "#FFE0A3")),
        (r.package_call, CALL_COL[r.package_call]),
    ]
    for j, (txt, fc) in enumerate(cells):
        ax.add_patch(plt.Rectangle((x0[j] + 0.02, i - 0.44), widths[j] - 0.04, 0.88, fc=fc, ec="none"))
        ax.text(x0[j] + widths[j] / 2, i, txt, ha="center", va="center", fontsize=6.3,
                color="grey" if txt == "not assessed" else text_colour(fc))
fig.savefig(FIG / "fig3_prioritization.png", bbox_inches="tight")
plt.close(fig)
print("Fig 3 / CSV: shortlist =", ", ".join(tri.dropna(subset=["shortlist_rank"]).gene))
print("Outputs written to", HERE)
