#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GEO_FILE="${ROOT_DIR}/mesh/unit_square.geo"
MSH_FILE="${ROOT_DIR}/mesh/unit_square.msh"

gmsh -2 "${GEO_FILE}" -format msh2 -o "${MSH_FILE}"
echo "Generated ${MSH_FILE}"
