#pragma once
#include "core/grid.hpp"

//----------------------- BOUNDARY CONDITIONS -----------------------------

// Dirichlet boundary condition p = 0 on physical walls.
// The wall sits half a cell between ghost and real cell, because pressure lives
// on cell centers: (p(0,j) + p(1,j)) / 2 = 0  =>  p(0,j) = -p(1,j).
void applyDirichletZero(CellField& p);

// Von Neumann: the derivative along the normal is 0 on the wall, so the ghost
// cell takes the same value as the border real cell.
void applyPressureNeumann(CellField& p);

// Velocity BCs: no trespassing + no sliding (ghost mirror).
void applyVelocityBC(FaceX& u, FaceY& v);
void applyVelocityUBC(FaceX& u);
void applyVelocityVBC(FaceY& v);

// Channel version: inlet u = U_in, outlet zero-gradient, walls no-slip when
// nu > 0 and free-slip when nu <= 0.
void applyChannelVelocityBC(FaceX& u, FaceY& v, Real U_in,
                            bool set_outlet_face, bool slip);

// Pressure BC for the channel: Neumann on inlet and walls, p = 0 at the outlet.
void applyChannelPressureBC(CellField& p);
