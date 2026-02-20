#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MESH_DIR="$ROOT_DIR/mesh"

if ! command -v gmsh >/dev/null 2>&1; then
  echo "gmsh not found in PATH" >&2
  exit 1
fi

if [ ! -d "$MESH_DIR" ]; then
  echo "mesh directory not found: $MESH_DIR" >&2
  exit 1
fi

for geo in "$MESH_DIR"/*.geo; do
  [ -e "$geo" ] || continue
  msh="${geo%.geo}.msh"
  echo "Generating $(basename "$msh")"
  gmsh -2 -format msh2 "$geo" -o "$msh"
done

echo "Mesh generation complete."
