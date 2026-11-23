#pragma once
#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <raylib.h>
#include <raymath.h>
#include <cassert>
#include <algorithm>
#include "grid.hpp"
#include "mouse.hpp"

//----------------------- BOUNDARY CONDITIONS -----------------------------

// Dirichlet boundary condition p = 0 on physical walls 
// Wall is half a cell between ghost and real cell because preassure lives on cell centers
// (p(0,j)+p(1,j)) / 2 = 0 therefore p(0,j) = - p(1,j).
//with this, ghost cells have the same value as the border real cells
inline void applyDirichletZero(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0,j)         = -p(1,j);
        p(p.nx+1, j)   = -p(p.nx, j);
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = -p(i, 1);
        p(i, p.ny + 1) = -p(i, p.ny);
    }
}

//Von Neumann boundary condition
//The derivative with respect to time must be 0 on the walls, so the gradient between the ghost cells
//and the next real cell must be 0, therefore the ghost cell value must be equal to the real cell value
inline void applyPressureNeumann(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0,j)         = p(1,j);
        p(p.nx+1, j)   = p(p.nx, j);
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = p(i, 1);
        p(i, p.ny + 1) = p(i, p.ny);
    }
}

//Boundary conditions for the velocity, no tresspassing
inline void applyVelocityBC(FaceX& u, FaceY& v) {
    //No tresspassing meaning velocity = 0 on real walls
    for (int j = 1; j <= u.ny; ++j) {
        u(1, j)    = 0.0;
        u(u.nx, j) = 0.0;
    }
    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)    = 0.0;
        v(i, v.ny) = 0.0;
    }

    //No sliding condition, tangencial velocity on boundary = 0
    for (int i = 0; i <= u.nx; ++i) {
        u(i, 0)      = -u(i, 1);
        u(i, u.ny+1) = -u(i, u.ny);
    }
    for (int j = 0; j <= v.ny; ++j) {
        v(0, j)      = -v(1, j);
        v(v.nx+1, j) = -v(v.nx, j);
    }
}
//Boundary conditions for velocity on x and y separately
inline void applyVelocityUBC(FaceX& u) {
        for (int j = 1; j <= u.ny; ++j) {
        u(1, j)    = 0.0;
        u(u.nx, j) = 0.0;
    }
    for (int i = 0; i <= u.nx; ++i) {
        u(i, 0)      = -u(i, 1);
        u(i, u.ny+1) = -u(i, u.ny);
    }
}

inline void applyVelocityVBC(FaceY& v) {
    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)    = 0.0;
        v(i, v.ny) = 0.0;
    }
    for (int j = 0; j <= v.ny; ++j) {
        v(0, j)      = -v(1, j);
        v(v.nx+1, j) = -v(v.nx, j);
    }
}

// Velocity BC for channel
inline void applyChannelVelocityBC(FaceX& u, FaceY& v, double U_in,
                                   bool set_outlet_face, double nu) {
    bool free_slip = (nu <= 0.0);

    for (int j = 1; j <= u.ny; ++j) {
        u(1, j) = U_in;
        u(0, j) = U_in; 
    }
    for (int j = 1; j <= v.ny; ++j) {
        v(0, j) = -v(1, j); 
    }

    for (int j = 1; j <= u.ny; ++j) {
        if (set_outlet_face) {
            u(u.nx, j) = u(u.nx - 1, j);
        }
        u(u.nx + 1, j) = u(u.nx, j); 
    }
    for (int j = 1; j <= v.ny; ++j) {
        v(v.nx + 1, j) = v(v.nx, j);
    }

    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)        = 0.0;
        v(i, v.ny)     = 0.0;
        v(i, 0)        = 0.0;
        v(i, v.ny + 1) = 0.0;
    }

    for (int i = 0; i <= u.nx + 1; ++i) {
        if (free_slip) {
            u(i, 0)       = u(i, 1);
            u(i, u.ny +1) = u(i, u.ny);
        } else {
        u(i, 0)        = -u(i, 1);
        u(i, u.ny + 1) = -u(i, u.ny);
        }
    }
}

// Pressure BC for channel
inline void applyChannelPressureBC(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0, j)        =  p(1, j);      
        p(p.nx + 1, j) = -p(p.nx, j);   
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = p(i, 1);       
        p(i, p.ny + 1) = p(i, p.ny);    
    }
}