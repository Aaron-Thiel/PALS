#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# PALS Pipeline — Database Download Script
# ============================================================================
# Downloads all required databases for the PALS pipeline.
#
# Databases:
#   1. Bakta       — bacterial genome annotation (~30 GB full, ~1.4 GB light)
#   2. CheckM2     — genome quality assessment (~3.5 GB)
#   3. BUSCO       — genome completeness (lineage-specific, ~300 MB)
#   4. Kraken2     — taxonomic read classification (~8 GB standard-8)
#   5. Sourmash    — genus-level classification via GTDB (~2 GB)
#   6. EggNOG      — functional annotation with eggnog-mapper (~45 GB)
#   7. Pfam-A      — protein family HMMs for BiG-SCAPE 2 (~300 MB)
#
# antiSMASH databases are bundled inside the Docker image and do NOT need
# a separate download.
#
# Usage:
#   ./download_databases.sh [DATABASE_DIR]
#
# Default DATABASE_DIR: databases/ (relative to project root)
#
# DISCLAIMER: Download URLs, version numbers, and size estimates are provided
# as-is with no warranty of correctness or completeness. External databases are
# maintained by their respective projects and may change, move, or become
# unavailable at any time. Downloading and using these databases is at your own
# risk. Refer to each database's official documentation and license terms.
# ============================================================================

DB_BASE="${1:-databases}"
SKIP_ALL=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ============================================================================
# Helper functions
# ============================================================================

info()  { echo -e "${CYAN}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
err()   { echo -e "${RED}[ERROR]${NC} $*"; }

# Prompt Yes/No/All. Returns 0 for yes, 1 for no.
# Sets SKIP_ALL=true if user picks "All".
ask_download() {
    local name="$1"
    local size="$2"
    local dest="$3"

    if $SKIP_ALL; then
        return 0
    fi

    echo ""
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "  Database:  ${YELLOW}${name}${NC}"
    echo -e "  Size:      ~${size}"
    echo -e "  Location:  ${dest}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    while true; do
        read -rp "Download ${name}? [Y]es / [N]o / [A]ll: " answer
        case "${answer,,}" in
            y|yes)  return 0 ;;
            n|no)   return 1 ;;
            a|all)  SKIP_ALL=true; return 0 ;;
            *)      echo "  Please enter Y, N, or A." ;;
        esac
    done
}

check_command() {
    if ! command -v "$1" &>/dev/null; then
        err "'$1' is required but not found. Please install it first."
        exit 1
    fi
}

# ============================================================================
# Pre-flight checks
# ============================================================================

echo ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  PALS Pipeline — Database Downloader${NC}"
echo -e "${CYAN}============================================================${NC}"
echo -e "  Target directory: ${YELLOW}${DB_BASE}${NC}"
echo ""

check_command wget
check_command tar

mkdir -p "${DB_BASE}"

# ============================================================================
# 1. Bakta database
# ============================================================================

BAKTA_DIR="${DB_BASE}/bakta/db"
BAKTA_TYPE=""

download_bakta() {
    if [ -d "${BAKTA_DIR}" ] && [ -f "${BAKTA_DIR}/version.json" ]; then
        ok "Bakta database already exists at ${BAKTA_DIR} — skipping."
        return
    fi

    if ! ask_download "Bakta" "30 GB (full) / 1.4 GB (light)" "${BAKTA_DIR}"; then
        warn "Skipping Bakta database."
        return
    fi

    # Ask for full vs light
    while true; do
        read -rp "  Download [F]ull (~30 GB) or [L]ight (~1.4 GB) version? " bt
        case "${bt,,}" in
            f|full)  BAKTA_TYPE="full"; break ;;
            l|light) BAKTA_TYPE="light"; break ;;
            *)       echo "  Please enter F or L." ;;
        esac
    done

    mkdir -p "${DB_BASE}/bakta"

    if command -v bakta_db &>/dev/null; then
        info "Downloading Bakta ${BAKTA_TYPE} database via bakta_db..."
        if [ "${BAKTA_TYPE}" = "light" ]; then
            bakta_db download --output "${DB_BASE}/bakta" --type light
        else
            bakta_db download --output "${DB_BASE}/bakta" --type full
        fi
    else
        # Fallback: download via Zenodo
        info "bakta_db not found — downloading Bakta ${BAKTA_TYPE} database from Zenodo..."

        if [ "${BAKTA_TYPE}" = "light" ]; then
            BAKTA_URL="https://zenodo.org/records/10522951/files/db-light.tar.gz"
        else
            BAKTA_URL="https://zenodo.org/records/10522951/files/db.tar.gz"
        fi

        wget -q --show-progress -O "${DB_BASE}/bakta/db.tar.gz" "${BAKTA_URL}"
        info "Extracting Bakta database..."
        tar -xzf "${DB_BASE}/bakta/db.tar.gz" -C "${DB_BASE}/bakta/"
        rm -f "${DB_BASE}/bakta/db.tar.gz"
    fi

    ok "Bakta database installed at ${BAKTA_DIR}"
}

# ============================================================================
# 2. CheckM2 database
# ============================================================================

CHECKM2_DIR="${DB_BASE}/checkm2"
CHECKM2_FILE="${CHECKM2_DIR}/uniref100.KO.1.dmnd"

download_checkm2() {
    if [ -f "${CHECKM2_FILE}" ]; then
        ok "CheckM2 database already exists at ${CHECKM2_DIR} — skipping."
        return
    fi

    if ! ask_download "CheckM2" "3.5 GB" "${CHECKM2_DIR}"; then
        warn "Skipping CheckM2 database."
        return
    fi

    mkdir -p "${CHECKM2_DIR}"

    if command -v checkm2 &>/dev/null; then
        info "Downloading CheckM2 database via checkm2 CLI..."
        checkm2 database --download --path "${CHECKM2_DIR}"
    elif command -v docker &>/dev/null; then
        info "checkm2 not found locally — downloading via Docker container..."
        docker run --rm \
            -v "${CHECKM2_DIR}:/checkm2_db" \
            staphb/checkm2:latest \
            checkm2 database --download --path /checkm2_db
    else
        err "Neither checkm2 nor Docker is available."
        err "Install CheckM2 or Docker, then re-run."
        return 1
    fi

    ok "CheckM2 database installed at ${CHECKM2_DIR}"
}

# ============================================================================
# 3. BUSCO lineage database
# ============================================================================

BUSCO_DIR="${DB_BASE}/busco"
BUSCO_LINEAGE=""
BUSCO_LINEAGE_DIR=""

download_busco() {
    # Check if default lineage already exists
    if [ -d "${BUSCO_DIR}/lineages/lactobacillaceae_odb12" ]; then
        ok "BUSCO lineage lactobacillaceae_odb12 already exists — skipping."
        BUSCO_LINEAGE="lactobacillaceae_odb12"
        BUSCO_LINEAGE_DIR="${BUSCO_DIR}/lineages/${BUSCO_LINEAGE}"
        return
    fi

    if ! ask_download "BUSCO lineage" "300 MB" "${BUSCO_DIR}"; then
        warn "Skipping BUSCO database."
        return
    fi

    # Let user choose lineage
    echo ""
    echo "  Available BUSCO lineages:"
    echo "    [1] lactobacillaceae_odb12  — Lactobacillaceae family (default)"
    echo "    [2] lactobacillales_odb10   — Lactobacillales order (broader)"
    echo "    [3] bacilli_odb10           — Bacilli class"
    echo "    [4] bacteria_odb10          — All bacteria (most general)"
    echo ""

    while true; do
        read -rp "  Select lineage [1-4, default=1]: " choice
        case "${choice:-1}" in
            1) BUSCO_LINEAGE="lactobacillaceae_odb12"; break ;;
            2) BUSCO_LINEAGE="lactobacillales_odb10"; break ;;
            3) BUSCO_LINEAGE="bacilli_odb10"; break ;;
            4) BUSCO_LINEAGE="bacteria_odb10"; break ;;
            *) echo "  Please enter 1, 2, 3, or 4." ;;
        esac
    done

    BUSCO_LINEAGE_DIR="${BUSCO_DIR}/lineages/${BUSCO_LINEAGE}"

    if [ -d "${BUSCO_LINEAGE_DIR}" ]; then
        ok "BUSCO lineage ${BUSCO_LINEAGE} already exists — skipping."
        return
    fi

    mkdir -p "${BUSCO_DIR}/lineages"

    info "Downloading BUSCO lineage ${BUSCO_LINEAGE}..."

    # Try via busco CLI first
    if command -v busco &>/dev/null; then
        busco --download_path "${BUSCO_DIR}" --download "${BUSCO_LINEAGE}"
    elif command -v docker &>/dev/null; then
        info "busco not found locally — trying via Docker container..."
        docker run --rm \
            -v "${BUSCO_DIR}:/busco_db" \
            staphb/busco:latest \
            busco --download_path /busco_db --download "${BUSCO_LINEAGE}"
    else
        # Fallback: direct download from BUSCO data server
        info "Trying direct download from BUSCO data server..."
        BUSCO_URL="https://busco-data.ezlab.org/v5/data/lineages/${BUSCO_LINEAGE}.2024-01-08.tar.gz"

        if wget -q --show-progress -O "${BUSCO_DIR}/${BUSCO_LINEAGE}.tar.gz" "${BUSCO_URL}" 2>/dev/null; then
            info "Extracting BUSCO lineage..."
            tar -xzf "${BUSCO_DIR}/${BUSCO_LINEAGE}.tar.gz" -C "${BUSCO_DIR}/lineages/"
            rm -f "${BUSCO_DIR}/${BUSCO_LINEAGE}.tar.gz"
        else
            err "Could not download BUSCO lineage. Install busco or Docker."
            return 1
        fi
    fi

    ok "BUSCO lineage ${BUSCO_LINEAGE} installed at ${BUSCO_LINEAGE_DIR}"

    if [ "${BUSCO_LINEAGE}" != "lactobacillaceae_odb12" ]; then
        warn "The pipeline defaults to lactobacillaceae_odb12."
        warn "Update the -l flag in modules/qc/busco.nf to use ${BUSCO_LINEAGE}."
    fi
}

# ============================================================================
# 4. Kraken2 database
# ============================================================================

KRAKEN2_DIR="${DB_BASE}/kraken2"

download_kraken2() {
    if [ -f "${KRAKEN2_DIR}/hash.k2d" ]; then
        ok "Kraken2 database already exists at ${KRAKEN2_DIR} — skipping."
        return
    fi

    if ! ask_download "Kraken2" "8–70 GB (depends on variant)" "${KRAKEN2_DIR}"; then
        warn "Skipping Kraken2 database."
        return
    fi

    # Let user choose database variant
    echo ""
    echo "  Available Kraken2 databases (pre-built, from genome-idx.s3.amazonaws.com):"
    echo "    [1] Standard-8   — ~8 GB   RefSeq bacteria/archaea/viral/human, capped (default)"
    echo "    [2] Standard-16  — ~16 GB  Same taxa, higher resolution"
    echo "    [3] Standard     — ~70 GB  Full RefSeq, maximum sensitivity"
    echo "    [4] PlusPF-8     — ~8 GB   Standard + protozoa + fungi, capped"
    echo ""
    echo "  All variants include comprehensive bacterial RefSeq coverage."
    echo "  Standard-8 is sufficient for Lactobacillaceae-focused work."
    echo ""

    local KRAKEN2_URL=""
    while true; do
        read -rp "  Select database [1-4, default=1]: " choice
        case "${choice:-1}" in
            1)
                KRAKEN2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz"
                info "Selected: Standard-8 (~8 GB)"
                break
                ;;
            2)
                KRAKEN2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz"
                info "Selected: Standard-16 (~16 GB)"
                break
                ;;
            3)
                KRAKEN2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20241228.tar.gz"
                info "Selected: Standard (~70 GB)"
                break
                ;;
            4)
                KRAKEN2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20241228.tar.gz"
                info "Selected: PlusPF-8 (~8 GB)"
                break
                ;;
            *)
                echo "  Please enter 1, 2, 3, or 4."
                ;;
        esac
    done

    mkdir -p "${KRAKEN2_DIR}"

    info "Downloading Kraken2 pre-built database..."
    info "Source: https://benlangmead.github.io/aws-indexes/k2"
    info "This may take a while depending on your connection speed..."

    if wget -q --show-progress -O "${KRAKEN2_DIR}/kraken2_db.tar.gz" "${KRAKEN2_URL}" 2>/dev/null; then
        info "Extracting Kraken2 database..."
        tar -xzf "${KRAKEN2_DIR}/kraken2_db.tar.gz" -C "${KRAKEN2_DIR}/"
        rm -f "${KRAKEN2_DIR}/kraken2_db.tar.gz"
    else
        warn "Pre-built download failed. You can build a Kraken2 database manually:"
        warn "  kraken2-build --standard --db ${KRAKEN2_DIR}"
        warn "Or download a pre-built database from:"
        warn "  https://benlangmead.github.io/aws-indexes/k2"
        return 1
    fi

    ok "Kraken2 database installed at ${KRAKEN2_DIR}"
}

# ============================================================================
# 5. Sourmash GTDB database
# ============================================================================

SOURMASH_DIR="${DB_BASE}/sourmash"
SOURMASH_DB="${SOURMASH_DIR}/gtdb-reps-rs226-k31.dna.zip"
SOURMASH_LINEAGES="${SOURMASH_DIR}/gtdb-rs226.lineages.csv"

download_sourmash() {
    if [ -f "${SOURMASH_DB}" ] && [ -f "${SOURMASH_LINEAGES}" ]; then
        ok "Sourmash GTDB database already exists at ${SOURMASH_DIR} — skipping."
        return
    fi

    if ! ask_download "Sourmash GTDB RS226" "2 GB" "${SOURMASH_DIR}"; then
        warn "Skipping Sourmash database."
        return
    fi

    mkdir -p "${SOURMASH_DIR}"

    # Download GTDB RS226 representatives (k=31)
    if [ ! -f "${SOURMASH_DB}" ]; then
        info "Downloading GTDB RS226 representative genomes (k=31)..."
        wget -q --show-progress -O "${SOURMASH_DB}" \
            "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip"
    fi

    # Download GTDB RS226 lineages CSV
    if [ ! -f "${SOURMASH_LINEAGES}" ]; then
        info "Downloading GTDB RS226 lineages CSV..."

        # Try the .gz version first
        if wget -q --show-progress -O "${SOURMASH_LINEAGES}.gz" \
            "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv.gz" 2>/dev/null; then
            gunzip "${SOURMASH_LINEAGES}.gz"
        else
            # Try uncompressed
            wget -q --show-progress -O "${SOURMASH_LINEAGES}" \
                "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv"
        fi
    fi

    ok "Sourmash GTDB RS226 database installed at ${SOURMASH_DIR}"
}

# ============================================================================
# 6. EggNOG database
# ============================================================================

EGGNOG_DIR="${DB_BASE}/eggnog"

download_eggnog() {
    if [ -f "${EGGNOG_DIR}/eggnog.db" ]; then
        ok "EggNOG database already exists at ${EGGNOG_DIR} — skipping."
        return
    fi

    if ! ask_download "EggNOG" "45 GB" "${EGGNOG_DIR}"; then
        warn "Skipping EggNOG database."
        return
    fi

    mkdir -p "${EGGNOG_DIR}"

    if command -v download_eggnog_data.py &>/dev/null; then
        info "Downloading EggNOG database via download_eggnog_data.py..."
        download_eggnog_data.py --data_dir "${EGGNOG_DIR}" -y
    else
        # Try running inside the Docker container
        info "download_eggnog_data.py not found locally."
        info "Attempting download via eggnog-mapper Docker container..."

        if command -v docker &>/dev/null; then
            docker run --rm \
                -v "${EGGNOG_DIR}:/data/eggnog" \
                nanozoo/eggnog-mapper:2.1.13--c16a7d2 \
                download_eggnog_data.py --data_dir /data/eggnog -y
        else
            err "Neither download_eggnog_data.py nor Docker is available."
            err "Install eggnog-mapper or Docker, then re-run."
            return 1
        fi
    fi

    ok "EggNOG database installed at ${EGGNOG_DIR}"
}

# ============================================================================
# 7. Pfam-A HMM database (for BiG-SCAPE 2)
# ============================================================================

PFAM_DIR="${DB_BASE}/antismash/pfam/35.0"
PFAM_FILE="${PFAM_DIR}/Pfam-A.hmm"

download_pfam() {
    if [ -f "${PFAM_FILE}" ]; then
        ok "Pfam-A database already exists at ${PFAM_FILE} — skipping."
        return
    fi

    if ! ask_download "Pfam-A 35.0" "300 MB" "${PFAM_DIR}"; then
        warn "Skipping Pfam-A database."
        return
    fi

    check_command gunzip

    mkdir -p "${PFAM_DIR}"

    info "Downloading Pfam-A.hmm.gz from EBI (release 35.0)..."
    wget -q --show-progress -O "${PFAM_DIR}/Pfam-A.hmm.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"

    info "Extracting Pfam-A.hmm..."
    gunzip "${PFAM_DIR}/Pfam-A.hmm.gz"

    # Press HMM if hmmpress is available (optional but speeds up BiG-SCAPE)
    if command -v hmmpress &>/dev/null; then
        info "Running hmmpress on Pfam-A.hmm..."
        hmmpress "${PFAM_FILE}"
    else
        warn "hmmpress not found — BiG-SCAPE will press the HMM on first run (slower)."
    fi

    ok "Pfam-A 35.0 installed at ${PFAM_FILE}"
}

# ============================================================================
# Run downloads
# ============================================================================

download_bakta
download_checkm2
download_busco
download_kraken2
download_sourmash
download_eggnog
download_pfam

# ============================================================================
# Summary
# ============================================================================

echo ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  Download Summary${NC}"
echo -e "${CYAN}============================================================${NC}"

check_db() {
    local label="$1"
    local path="$2"
    if [ -e "$path" ]; then
        echo -e "  ${GREEN}✓${NC} ${label}  →  ${path}"
    else
        echo -e "  ${RED}✗${NC} ${label}  →  not found"
    fi
}

check_db "Bakta"      "${BAKTA_DIR}/version.json"
check_db "CheckM2"    "${CHECKM2_FILE}"
check_db "BUSCO"      "${BUSCO_LINEAGE_DIR:-${BUSCO_DIR}/lineages}"
check_db "Kraken2"    "${KRAKEN2_DIR}/hash.k2d"
check_db "Sourmash"   "${SOURMASH_DB}"
check_db "EggNOG"     "${EGGNOG_DIR}/eggnog.db"
check_db "Pfam-A"     "${PFAM_FILE}"
echo ""

echo -e "  ${YELLOW}Note:${NC} antiSMASH databases are expected to be bundled in the"
echo -e "        Docker image (nanozoo/antismash:8.0.0). If antiSMASH fails"
echo -e "        with database errors, see the antiSMASH docs for manual setup."
echo ""
echo -e "  Expected directory layout:"
echo -e "    ${DB_BASE}/"
echo -e "    ├── bakta/db/                  # Bakta annotation database"
echo -e "    ├── checkm2/                   # CheckM2 quality assessment"
echo -e "    ├── busco/lineages/            # BUSCO lineage datasets"
echo -e "    ├── kraken2/                   # Kraken2 taxonomy database"
echo -e "    ├── sourmash/                  # Sourmash GTDB sketches"
echo -e "    ├── eggnog/                    # EggNOG functional annotation"
echo -e "    └── antismash/pfam/35.0/       # Pfam-A HMMs (for BiG-SCAPE)"
echo ""
echo -e "  ${YELLOW}Disclaimer:${NC} URLs and sizes are approximate and provided as-is."
echo -e "  Databases are maintained by their respective projects and may change."
echo -e "  Downloading and using these databases is at your own risk."
echo ""
ok "Done."
