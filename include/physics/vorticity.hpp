#pragma once
#include <raylib.h>
#include "core/grid.hpp"

//exact vorticity on a cell corner
inline Real cornerVorticity(const FaceX& u, const FaceY& v, int i, int j, Real h) {
    return ((v(i, j) - v(i - 1, j)) - (u(i, j) - u(i, j - 1))) / h;
}

//Vorticity on cell centers (average of the four corners)
void computeVorticity(CellField& omega, const FaceX& u, const FaceY& v, Real h);

//Colour maps: vorticity blue -> white -> red, speed black -> cyan -> yellow -> red
Color vorticityToColor(Real w, Real max_w);
Color speedToColor(Real speed, Real max_speed);
