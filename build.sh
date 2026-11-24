#!/usr/bin/env bash
# Build f1 — 2D incompressible Navier-Stokes, real-time, raylib.
#   ./build.sh                 -> Release con OpenMP
#   ./build.sh -DUSE_OPENMP=OFF -> sin OpenMP
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

cmake -B build -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build build -j

cat <<'EOF'

built ./f1
use:  ./f1
EOF
