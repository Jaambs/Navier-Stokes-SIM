#!/usr/bin/env bash
# static linking compiler script, portable version of the .exe
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

RAYLIB=~/.cache/raylib_win/raylib-5.5_win64_mingw-w64
# raylib estatico para Windows: se baja solo la primera vez
if [ ! -f "$RAYLIB/lib/libraylib.a" ]; then
    mkdir -p ~/.cache/raylib_win
    curl -sSL -o /tmp/raylib.zip https://github.com/raysan5/raylib/releases/download/5.5/raylib-5.5_win64_mingw-w64.zip
    unzip -oq /tmp/raylib.zip -d ~/.cache/raylib_win
fi

SOURCES="src/main.cpp src/BC.cpp src/cylinder.cpp src/mouse.cpp
         src/solver.cpp src/vorticity.cpp src/multigrid.cpp"

/mnt/c/TDM-GCC-64/bin/g++.exe -O3 -std=c++17 -march=native -mwindows $SOURCES -o f1_fluid.exe \
    -I"$(wslpath -w ./include)" \
    -I"$(wslpath -w "$RAYLIB/include")" -L"$(wslpath -w "$RAYLIB/lib")" \
    -static -static-libgcc -static-libstdc++ \
    -lraylib -lopengl32 -lgdi32 -lwinmm -lshell32 -luser32 -lcomdlg32

echo "built ./f1_fluid.exe"
