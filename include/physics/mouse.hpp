#pragma once
#include "core/grid.hpp"

// Drops dye and applies a force in a smooth blob around the mouse position.
void applyMouseBrush(FaceX& u, FaceY& v, CellField& dyeR, CellField& dyeB, CellField& dyeG,
                     Real mx, Real my, Real mdx, Real mdy,
                     Real r_color, Real b_color, Real g_color,
                     Real dt, Real h, Real force_strengh);
