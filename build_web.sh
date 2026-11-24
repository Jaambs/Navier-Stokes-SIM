#!/usr/bin/env bash
# Build for the browser with Emscripten. Produces web/index.{html,js,wasm}.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# emscripten and the wasm build of raylib live outside the project
export PATH="$HOME/emsdk/upstream/emscripten:$HOME/emsdk/upstream/bin:$PATH"
RAYLIB="$HOME/.cache/raylib_web/src"

SOURCES="src/main.cpp src/BC.cpp src/cylinder.cpp src/mouse.cpp
         src/solver.cpp src/vorticity.cpp src/multigrid.cpp"

mkdir -p web

em++ -O2 -std=c++17 $SOURCES -o web/index.html \
    -Iinclude -I"$RAYLIB" -DPLATFORM_WEB \
    "$RAYLIB/libraylib.web.a" \
    -sUSE_GLFW=3 -sASYNCIFY \
    -sINITIAL_MEMORY=67108864 -sALLOW_MEMORY_GROWTH=1 \
    -sGL_ENABLE_GET_PROC_ADDRESS \
    --shell-file minshell.html

echo "built web/  ->  python3 -m http.server 8000 --directory web"
