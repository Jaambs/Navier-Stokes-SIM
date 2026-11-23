#!/usr/bin/env bash
# Build f1_fluid — 2D incompressible Navier-Stokes, real-time, raylib.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

echo "raylib $(pkg-config --modversion raylib)"
g++ -O3 -std=c++17 -march=native -DNDEBUG -Wall -Wextra -Wshadow main.cpp -o f1 \
    $(pkg-config --cflags --libs raylib) -lm -lpthread -ldl
cat <<'EOF'

built ./f1
use:  ./f1 
EOF
