#pragma once
#include <algorithm>
#include <cmath>
#include "core/grid.hpp"
#include "physics/cylinder.hpp"
#include "physics/multigrid.hpp"

//----------------------------------------------------
// Generic bilinear interpolation (template: must live in the header)
template <typename Field>
inline Real sample(const Field& f, Real gx, Real gy) {
    //we take only the valid cell range without ghost cells
    gx = std::clamp(gx, 1.0f, static_cast<Real>(f.nx));
    gy = std::clamp(gy, 1.0f, static_cast<Real>(f.ny));

    //find left inferior corner
    int i0 = static_cast<int>(std::floor(gx));
    int j0 = static_cast<int>(std::floor(gy));
    //protection against gx or gy equal to f.nx / f.ny so i0+1 is out of range
    if (i0 >= f.nx) i0 = f.nx - 1;
    if (j0 >= f.ny) j0 = f.ny - 1;

    //relative distances
    Real s = gx - i0;
    Real t = gy - j0;

    //Interpolation
    return (1.0f - s) * (1.0f - t) * f(i0,     j0    )
     +        s       * (1.0f - t) * f(i0 + 1, j0    )
     + (1.0f - s) *         t      * f(i0,     j0 + 1)
     +        s       *      t      * f(i0 + 1, j0 + 1);
}

// ----------- ADVECTION -------------
void advectU(FaceX& un, const FaceX& u, const FaceY& v, Real dt, Real h);
void advectV(FaceY& vn, const FaceX& u, const FaceY& v, Real dt, Real h);
void advectScalar(CellField& d, const CellField& d0,
                  const FaceX& u, const FaceY& v, Real dt, Real h);

// -------------- DIVERGENCE --------------
void divergence(CellField& div, const FaceX& u, const FaceY& v, Real h);
Real maxAbs(const CellField& f);

// -------------- POISSON SOLVERS --------------
// Gauss-Seidel red-black (parallel)
void gaussSeidelChannelRB(CellField& p, const CellField& rhs, const CylinderMask& cyl,
                          Real h, int iters);

// -------------- PROJECTION --------------
void projectChannel(FaceX& u, FaceY& v, CellField& p, CellField& div, CylinderMask& cyl,
                    CellField& rhs, Real dt, Real h, Real U_in, int iters, bool slip);
void projectChannelMG(FaceX& u, FaceY& v, CellField& p, CellField& div, const CylinderMask& cyl,
                      MultigridSolver& mg, Real dt, Real h, Real U_in, bool slip);

// ------------- VISCOSITY DIFFUSION --------------
void diffuseChannelU(FaceX& u, const FaceX& u0,
                     Real nu, Real dt, Real h, Real U_in, int iters = 30);
void diffuseChannelV(FaceY& v, const FaceY& v0,
                     Real nu, Real dt, Real h, int iters = 30);

// ------------------------------- DEPRECATED / UNUSED ----------------
//Kept for comparison and for other purposes
void jacobi(CellField& p, const CellField& rhs, Real h, int iters);
void gaussSeidelNeumann(CellField& p, const CellField& rhs, Real h, Real iters);
void DiffuseU(FaceX& u, const FaceX& u0, Real nu, Real dt, Real h, int iters);
void DiffuseV(FaceY& v, const FaceY& v0, Real nu, Real dt, Real h, int iters);
void project(FaceX& u, FaceY& v, CellField& p, CellField& div, Real dt, Real h, int iter);
