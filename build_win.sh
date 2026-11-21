#!/usr/bin/env bash
# static linking compiler script, portable version of the .exe
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

RAYLIB=~/.cache/raylib_win/raylib-5.5_win64_mingw-w64
# raylib estático para Windows: se baja solo la primera vez
if [ ! -f "$RAYLIB/lib/libraylib.a" ]; then
    mkdir -p ~/.cache/raylib_win
    curl -sSL -o /tmp/raylib.zip https://github.com/raysan5/raylib/releases/download/5.5/raylib-5.5_win64_mingw-w64.zip
    unzip -oq /tmp/raylib.zip -d ~/.cache/raylib_win
fi

/mnt/c/TDM-GCC-64/bin/g++.exe -O2 -std=c++17 -mwindows main.cpp -o f1_fluid.exe \
    -I"$(wslpath -w "$RAYLIB/include")" -L"$(wslpath -w "$RAYLIB/lib")" \
    -static -static-libgcc -static-libstdc++ \
    -lraylib -lopengl32 -lgdi32 -lwinmm -lshell32 -luser32 -lcomdlg32

echo "built ./f1_fluid.exe"
