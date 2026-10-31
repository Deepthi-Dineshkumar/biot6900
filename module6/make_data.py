#!/usr/bin/env python3
"""
Build the four input files for a PDAC PRS:
  data/genotypes.csv  data/variants.csv  data/score.csv  data/meta.csv

Inputs (from build_inputs.sh):
  raw/<PGS>_hmPOS_GRCh37.txt.gz                  PGS Catalog harmonized score
  raw/genos.tsv                                  CHROM POS ID REF ALT GT...  (bcftools query)
  raw/samples.txt                                sample IDs, same order as GT columns
  raw/integrated_call_samples_v3.20130502.ALL.panel   1000G sample -> super-population

Genotypes, REF/ALT and score are REAL. Age and phenotype are SIMULATED (1000G has
no disease labels). Simulation assumptions are printed and written to
data/SIMULATION_NOTES.txt.
"""
import argparse, os
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--pgs", default="PGS000725")
p.add_argument("--raw", default="raw")
p.add_argument("--out", default="data")
p.add_argument("--seed", type=int, default=6900)
p.add_argument("--prevalence", type=float, default=0.30,
               help="case fraction per ancestry (case-control style, as in the class demo)")
p.add_argument("--or_per_sd", type=float, default=1.5,
               help="true OR per SD of genetic liability (PDAC PRS literature ~1.4-1.5)")
p.add_argument("--no_attenuation", action="store_true",
               help="sensitivity run: score tags causal signal equally well in every ancestry")
a = p.parse_args()
rng = np.random.default_rng(a.seed)
os.makedirs(a.out, exist_ok=True)

# ---------- 1. score (real) ----------
sc = pd.read_csv(os.path.join(a.raw, f"{a.pgs}_hmPOS_GRCh37.txt.gz"), sep="\t", comment="#", dtype=str)
chr_col = "hm_chr" if "hm_chr" in sc else "chr_name"
pos_col = "hm_pos" if "hm_pos" in sc else "chr_position"
sc = sc.dropna(subset=[chr_col, pos_col, "effect_allele", "effect_weight"])
sc["key"] = sc[chr_col].astype(str) + ":" + sc[pos_col].astype(str).str.replace(r"\.0$", "", regex=True)
rs = sc["hm_rsID"] if "hm_rsID" in sc else (sc["rsID"] if "rsID" in sc else pd.Series([""]*len(sc), index=sc.index))
sc["variant"] = np.where(rs.fillna("").str.startswith("rs"), rs, "var_" + sc["key"])
other = sc["other_allele"] if "other_allele" in sc else (sc["hm_inferOtherAllele"] if "hm_inferOtherAllele" in sc else pd.Series([""]*len(sc), index=sc.index))
score = pd.DataFrame({"variant": sc["variant"].values,
                      "effect_allele": sc["effect_allele"].str.upper().values,
                      "other_allele": other.fillna("").str.upper().values,
                      "weight": sc["effect_weight"].astype(float).values})
score = score.drop_duplicates("variant")
key2var = dict(zip(sc["key"], sc["variant"]))

# ---------- 2. genotypes + REF/ALT (real) ----------
samples = [s.strip() for s in open(os.path.join(a.raw, "samples.txt")) if s.strip()]
rows = []
with open(os.path.join(a.raw, "genos.tsv")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        key = f"{f[0]}:{f[1]}"
        if key not in key2var: continue
        gts = f[5:]
        dos = [g.count("1") if "." not in g else np.nan for g in gts]
        rows.append((key2var[key], f[3].upper(), f[4].upper(), dos))
seen, uniq = set(), []
for r in rows:                      # one record per score variant
    if r[0] not in seen: seen.add(r[0]); uniq.append(r)
variants = pd.DataFrame({"variant": [r[0] for r in uniq], "ref": [r[1] for r in uniq], "alt": [r[2] for r in uniq]})
geno = pd.DataFrame(np.array([r[3] for r in uniq], dtype=float).T, columns=variants["variant"])
geno.insert(0, "sample", samples)
geno = geno.fillna(geno.mean(numeric_only=True))   # 1000G phase 3 is fully phased; safety only

# ---------- 3. ancestry (real) ----------
panel = pd.read_csv(os.path.join(a.raw, "integrated_call_samples_v3.20130502.ALL.panel"), sep=r"\s+")
anc = panel.set_index("sample").loc[samples, "super_pop"].values

# ---------- 4. simulated age + phenotype ----------
# harmonize to the score's effect allele 
s = score.set_index("variant").loc[variants["variant"]]
eff = s["effect_allele"].values
keep = eff == variants["alt"].values
flip = (eff == variants["ref"].values) & ~keep
use = keep | flip
G = geno[variants["variant"]].to_numpy(float)
G[:, flip] = 2 - G[:, flip]
raw_prs = G[:, use] @ s["weight"].values[use]

# Tagging efficiency: correlation between the EUR-derived score and the true causal
# liability. Based on Martin et al. 2019 Nat Genet (PRS R^2 ~1.6x lower in SAS/AMR,
# ~2x lower in EAS, ~4.5x lower in AFR than EUR); rho = sqrt(1 / fold).
RHO = {"EUR": 1.00, "AMR": 0.79, "SAS": 0.79, "EAS": 0.71, "AFR": 0.47}
if a.no_attenuation: RHO = {k: 1.0 for k in RHO}

eur_sd = raw_prs[anc == "EUR"].std()
age = np.clip(rng.normal(65, 8, len(samples)), 40, 90)
beta_g, beta_age = np.log(a.or_per_sd), 0.06          # 0.06 log-OR per year of age
y = np.zeros(len(samples), int)
for pop, rho in RHO.items():
    m = anc == pop
    if not m.any(): continue
    z = (raw_prs[m] - raw_prs[m].mean()) / eur_sd       # centred within group: true risk is
    liab = rho * z + np.sqrt(max(0, 1 - rho**2)) * rng.normal(size=m.sum())  # not set by tag frequencies
    lin = beta_g * liab + beta_age * (age[m] - 65)
    lo, hi = -10, 10                                    # intercept giving the target prevalence
    for _ in range(60):
        b0 = (lo + hi) / 2
        if (1/(1+np.exp(-(b0+lin)))).mean() > a.prevalence: hi = b0
        else: lo = b0
    y[m] = (rng.uniform(size=m.sum()) < 1/(1+np.exp(-(b0+lin)))).astype(int)

meta = pd.DataFrame({"sample": samples, "ancestry": anc, "age": age.round(1), "phenotype": y})

# ---------- 5. write ----------
geno.to_csv(os.path.join(a.out, "genotypes.csv"), index=False)
variants.to_csv(os.path.join(a.out, "variants.csv"), index=False)
score.to_csv(os.path.join(a.out, "score.csv"), index=False)
meta.to_csv(os.path.join(a.out, "meta.csv"), index=False)

notes = f"""Score: {a.pgs} (PGS Catalog, GRCh37 harmonized) - {len(score)} variants; {len(variants)} found in 1000G phase 3
Genotypes: 1000 Genomes phase 3 (GRCh37), {len(samples)} people, biallelic SNPs, dosage = ALT count
Ancestry: 1000G super-population labels
SIMULATED: age ~ N(65, 8) clipped 40-90; phenotype from logistic model:
  logit = b0[group] + log({a.or_per_sd})*liability + 0.06*(age-65)
  liability = rho[group]*PRS_within_group + sqrt(1-rho^2)*noise   (PRS in EUR-SD units)
  rho = {RHO}  (tagging efficiency; Martin et al. 2019)
  b0[group] set so case fraction = {a.prevalence} in every group (equal true risk across groups)
Seed: {a.seed}
"""
open(os.path.join(a.out, "SIMULATION_NOTES.txt"), "w").write(notes)
print(notes)
print("case fraction by group:", meta.groupby("ancestry")["phenotype"].mean().round(3).to_dict())
