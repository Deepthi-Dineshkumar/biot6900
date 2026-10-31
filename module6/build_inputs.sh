#!/usr/bin/env bash
# Pull a PGS Catalog score + matching 1000 Genomes phase 3 genotypes, then build data/.
# Usage:  bash build_inputs.sh PGS000725      (needs bcftools, curl, wget, python3 + pandas)
set -euo pipefail
PGS=${1:-PGS000725}
G1K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
mkdir -p raw && cd raw

# 1. score (GRCh37 harmonized, to match 1000G phase 3)
wget -nc -q https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${PGS}/ScoringFiles/Harmonized/${PGS}_hmPOS_GRCh37.txt.gz
# 2. sample -> super-population panel
wget -nc -q ${G1K}/integrated_call_samples_v3.20130502.ALL.panel

# 3. positions to fetch (autosomes only)
zcat ${PGS}_hmPOS_GRCh37.txt.gz | grep -v '^#' | awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  {c=$h["hm_chr"]; p=$h["hm_pos"]; if (c ~ /^[0-9]+$/ && p ~ /^[0-9]+$/) print c"\t"p}' \
  | sort -k1,1n -k2,2n -u > regions.tsv
echo "positions to fetch: $(wc -l < regions.tsv)"

# 4. stream only those positions from the remote VCFs (tabix range requests, no full download)
LISTING=$(curl -s ${G1K}/)
: > genos.tsv
for c in $(cut -f1 regions.tsv | sort -nu); do
  VCF=$(echo "$LISTING" | grep -o "ALL\.chr${c}\.phase3[^\"<> ]*\.genotypes\.vcf\.gz" | head -1)
  echo "chr${c}: ${VCF}"
  [ -s samples.txt ] || bcftools query -l "${G1K}/${VCF}" > samples.txt
  awk -v c="$c" '$1==c' regions.tsv > reg_tmp.tsv
  bcftools view -R reg_tmp.tsv -m2 -M2 -v snps "${G1K}/${VCF}" -Ou \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' >> genos.tsv
done
rm -f reg_tmp.tsv
echo "SNP records fetched: $(wc -l < genos.tsv)"
cd ..

# 5. build the four CSVs
python3 make_data.py --pgs ${PGS}
