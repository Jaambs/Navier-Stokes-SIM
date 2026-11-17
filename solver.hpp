#pragma once
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
// Dirichlet boundary condition p = 0 on physical walls 
// Wall is half a cell between ghost and real cell because preassure lives on cell centers
// (p(0,j)+p(1,j)) / 2 = 0 therefore p(0,j) = - p(1,j).
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
//with this, ghost cells have the same value as the border real cells

//Generic Bilineal interpolation algorithm
template <typename Field>
inline double sample(const Field& f, double gx, double gy) {
    //we take only the valid cell range without ghost cells
    gx=std::clamp(gx, 1.0, static_cast<double>(f.nx));
    gy=std::clamp(gy, 1.0, static_cast<double>(f.ny));

    //find left inferior corner
    int i0= static_cast<int>(std::floor(gx));
    int j0= static_cast<int>(std::floor(gy));
    //protection against cases where gx or gy are equal to f.nx or f.ny so i0+1 ist greater than f.nx
    if (i0>=f.nx) i0 = f.nx-1;
    if (j0>=f.ny) j0 = f.ny-1;

    //relative distances
    double s= gx-i0;
    double t= gy-j0;

    //Interpolation
    return (1.0 - s) *  (1.0 - t) * f(i0,     j0    )
     +        s  *      (1.0 - t) * f(i0 + 1, j0    )
     + (1.0 - s) *             t  * f(i0,     j0 + 1)
     +        s  *             t  * f(i0 + 1, j0 + 1);
}

// ----------- ADVECTION -------------

//horizontal velocity advection (u on FaceX)
inline void advectU(FaceX& un, const FaceX& u, const FaceY& v, double dt, double h) {
    for (int i=1; i<=u.nx; ++i) {
        for (int j= 1; j <= u.ny; ++j) {
            //physical position and time step back
            double x = xu(i,h) - dt * u(i,j);
            double y = yu(j,h) - dt * vAtU(v, i, j);

            //convert x,y to grid index
            double gx = x / h + 1.0;
            double gy = y / h + 0.5;

            //Interpolacion bilineal
            un(i,j) = sample(u, gx, gy);
        }
    }
}

//vertical velocity advection (v on FaceY)
inline void advectV(FaceY& vn, const FaceX& u, const FaceY& v, double dt, double h) {
    for (int i=1; i<=v.nx; ++i) {
        for (int j= 1; j <= v.ny; ++j) {
            //physical position and time step back
            double x = xv(i, h) - dt * uAtV(u, i, j);
            double y = yv(j, h) - dt * v(i, j);

            //convert x,y to grid index
            double gx = x / h + 0.5;
            double gy = y / h + 1.0;

            //Interpolacion bilineal
            vn(i,j) = sample(v, gx, gy);
        }
    }
}

//scalar advection
inline void advectScalar(CellField& d, const CellField& d0,
                         const FaceX& u, const FaceY& v, double dt, double h) {
    for (int i = 1; i <= d.nx; ++i) {
        for (int j = 1; j <= d.ny; ++j) {
            double x = xc(i, h) - dt * uAtCenter(u, i, j);
            double y = yc(j, h) - dt * vAtCenter(v, i, j);

            // Convertir (x, y) a coordenadas de índice de centros
            double gx = x / h + 0.5;
            double gy = y / h + 0.5;

            d(i, j) = sample(d0, gx, gy);
        }
    }
}

// -------------- DIVERGENCE --------------

//Function to calculate the divergence inside the cell with center i+0.5, j+0.5
inline void divergence(CellField& div, const FaceX& u, const FaceY& v, double h) {
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            div(i, j) = (u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j)) / h;
        }
    }
}

// Function to calculate the max value of a field inside an scalar field
inline double maxAbs(const CellField& f) {
    double max_val = 0.0;
    for (int i = 1; i <= f.nx; ++i) {
        for (int j = 1; j <= f.ny; ++j) {
            max_val = std::max(max_val, std::abs(f(i, j)));
        }
    }
    return max_val;
}

// --------------- Poisson Solver -----------------

//Jacobi solver (will be unused just to compare with gauss-seidel)
inline void jacobi(CellField& p, const CellField&rhs, double h, int iters) {
    CellField p_new=p;
    double h2 = h*h;

    for (int k = 0; k < iters; ++k) {
        applyDirichletZero(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p_new(i,j) = 0.25 * (p(i - 1, j) + p(i + 1, j)
                                   + p(i, j - 1) + p(i, j + 1)
                                   - h2 * rhs(i, j));
            }
        }
        std::swap(p.a,p_new.a);
    }
    applyDirichletZero(p);
}

//GaussSeidel solver
inline void gaussSeidel(CellField&p, const CellField& rhs, double h, double iters) {
    double h2 = h * h;

    for (int k = 0; k <= iters; ++k) {
        applyDirichletZero(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p(i, j) = 0.25 * (p(i - 1, j) + p(i + 1, j)
                                + p(i, j - 1) + p(i, j + 1)
                                - h2 * rhs(i, j));
            }
        }
        applyDirichletZero(p);
    }
}