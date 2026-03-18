#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# ECOSSDB - Ecosystem Services Database - Installer
# ============================================================================
#
# Standalone install for ECOSSDB pipeline.
# When running within danaseq, use danaseq's install.sh instead (handles
# the ECOSSDB submodule automatically).
#
# Usage:
#   ./install.sh              # Set up conda env + parse CICES ontology
#   ./install.sh --env-only   # Only create conda environment
#   ./install.sh --skip-env   # Skip conda, just parse CICES if needed
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/envs"
SKIP_ENV=false
ENV_ONLY=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-env)  SKIP_ENV=true; shift ;;
        --env-only)  ENV_ONLY=true; shift ;;
        -h|--help)
            echo "Usage: ./install.sh [--env-only | --skip-env]"
            echo ""
            echo "  --env-only   Only create the conda environment"
            echo "  --skip-env   Skip conda env setup, just verify data files"
            echo ""
            echo "The ecossdb-core conda env is only needed for:"
            echo "  - Re-parsing the CICES xlsx (parse_cices.py requires openpyxl)"
            echo "  - Running ECOSSDB as a standalone Nextflow pipeline"
            echo ""
            echo "When running within danaseq, no separate env is needed — all"
            echo "pipeline scripts use Python stdlib only."
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Check conda/mamba
if ! command -v conda &>/dev/null && ! command -v mamba &>/dev/null; then
    echo "[ERROR] conda or mamba not found. Install Miniforge first."
    exit 1
fi
SOLVER=$(command -v mamba &>/dev/null && echo mamba || echo conda)

# 1. Conda environment (optional — only needed for standalone or CICES re-parsing)
if [[ "${SKIP_ENV}" == false ]]; then
    local_env="${ENV_DIR}/ecossdb-core"
    if [[ -d "${local_env}" ]]; then
        echo "[OK] ecossdb-core environment already exists"
    else
        echo "Creating ecossdb-core conda environment ..."
        ${SOLVER} env create -p "${local_env}" -f "${ENV_DIR}/ecossdb-core.yml" -y
        echo "[SUCCESS] ecossdb-core environment created at ${local_env}"
    fi
fi

if [[ "${ENV_ONLY}" == true ]]; then
    exit 0
fi

# 2. Verify CICES ontology (pre-parsed TSV should be in repo)
ontology_tsv="${SCRIPT_DIR}/db/ontology/cices_v5.2.tsv"
hierarchy_json="${SCRIPT_DIR}/db/ontology/es_hierarchy.json"
if [[ -f "${ontology_tsv}" && -f "${hierarchy_json}" ]]; then
    n_entries=$(wc -l < "${ontology_tsv}")
    echo "[OK] CICES 5.2 ontology: ${n_entries} entries"
else
    echo "[WARNING] CICES ontology not found. To generate it, run:"
    echo "  conda run -p ${ENV_DIR}/ecossdb-core python bin/parse_cices.py \\"
    echo "    --input /path/to/CICES_V5.2.xlsx \\"
    echo "    --output-tsv ${ontology_tsv} \\"
    echo "    --output-json ${hierarchy_json}"
fi

# 3. Verify mapping table
mapping="${SCRIPT_DIR}/db/mappings/es_gene_mapping.tsv"
if [[ -f "${mapping}" ]]; then
    n_mappings=$(($(wc -l < "${mapping}") - 1))
    echo "[OK] Gene-ES mapping table: ${n_mappings} mappings"
else
    echo "[WARNING] Mapping table not found at ${mapping}"
    echo "  Run: python bin/bootstrap_mapping.py --foam /path/to/FOAM-onto_rel1.tsv \\"
    echo "    --ontology ${ontology_tsv} --output ${mapping}"
fi

# 4. Verify SDG crosswalk
sdg_crosswalk="${SCRIPT_DIR}/db/ontology/sdg/cices_to_sdg.tsv"
if [[ -f "${sdg_crosswalk}" ]]; then
    n_links=$(($(wc -l < "${sdg_crosswalk}") - 1))
    echo "[OK] CICES-SDG crosswalk: ${n_links} links"
else
    echo "[WARNING] SDG crosswalk not found. Run: python bin/build_sdg_crosswalk.py"
fi

echo ""
echo "ECOSSDB is ready."
echo ""
echo "Standalone usage:"
echo "  nextflow run main.nf --annotations <path> --contig2bin <path> -profile test"
echo ""
echo "Within danaseq:"
echo "  Enabled by default (--run_ecossdb true). No separate setup needed."
