#!/bin/bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate
conda activate vcf2maf

# IO Vars
IN_SNP_VCF="${1}"
IN_INDEL_VCF="${2}"
MERGED_VCF="$(basename ${IN_SNP_VCF%%.vcf})_$(basename ${IN_INDEL_VCF}).gz"
FILTERED_VCF="${MERGED_VCF%%.vcf.gz}_filtered.vcf.gz"
OUT_MAF="${MERGED_VCF%%.vcf.gz}.maf"
TMP="${3}"

# References
REF_FASTA="${4}"
FILTER_VCF="${5}"

# Threading
MAX_THREADS="${6}"
VEP_THREADS=$(echo "${MAX_THREADS} - 2" | bc -l)
BCF_THREADS=$(echo "${MAX_THREADS} - 1" | bc -l)

# VEP variables
VEP="${HOME}/vep_hg19"
VEP_CACHE="${HOME}/vep_cache_hg19"
VEP_THREADS=4
VEP_BUFFER=5000
VEP_CACHE_VERSION=86

# Spiece variables
SPIECE="homo_sapiens"
NCBI_BUILD="GRCh37"

# Global functions
function index() {
  local VCF="${1}"

  bgzip --stdout "${VCF}" > "${VCF}.gz"
  tabix -p vcf "${VCF}.gz"
}

# Indexing vcf from varscan
index "${IN_SNP_VCF}"
index "${IN_INDEL_VCF}"

# Merging SNP and InDels
bcftools concat --allow-overlaps --rm-dups none --output "${MERGED_VCF}" --output-type z --threads "${BCF_THREADS}" "${IN_SNP_VCF}.gz" "${IN_INDEL_VCF}.gz"

bcftools index "${MERGED_VCF}"

bcftools filter --exclude 'ADP<=15|FORMAT/DP<=15|FORMAT/FREQ=1|FORMAT/PVAL<=0.05' --output "${FILTERED_VCF}" --output-type z --threads "${BCF_THREADS}" "${MERGED_VCF}"

bcftools index "${FILTERED_VCF}"

vcf2maf.pl --input-vcf "${FILTERED_VCF}" --output-maf "${OUT_MAF}" --tmp-dir "${TMP}" --vep-path "${VEP}" --vep-data "${VEP_CACHE}" --vep-forks "${VEP_THREADS}" --buffer-size "${VEP_BUFFER}" --cache-version "${VEP_CACHE_VERSION}" --species "${SPIECE}" --ncbi-build "${NCBI_BUILD}" --ref-fasta "${REF_FASTA}" --filter-vcf "${FILTER_VCF}"
