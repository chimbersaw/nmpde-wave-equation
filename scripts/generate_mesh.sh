#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UNIT_GEO_FILE="${ROOT_DIR}/mesh/unit_square.geo"
UNIT_MSH_FILE="${ROOT_DIR}/mesh/unit_square.msh"
RADIAL_GEO_FILE="${ROOT_DIR}/mesh/radial_square.geo"
RADIAL_MSH_FILE="${ROOT_DIR}/mesh/radial_square.msh"

gmsh -2 "${UNIT_GEO_FILE}" -format msh2 -o "${UNIT_MSH_FILE}"
echo "Generated ${UNIT_MSH_FILE}"

gmsh -2 "${RADIAL_GEO_FILE}" -format msh2 -o "${RADIAL_MSH_FILE}"
echo "Generated ${RADIAL_MSH_FILE}"
