#pragma once
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <raymath.h>
#include <algorithm>
#include "grid.hpp"
#include "mouse.hpp"
#include "BC.hpp"

//----------------------------------------------------

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
//GaussSeidel with Neumann condition
inline void gaussSeidelNeumann(CellField&p, const CellField& rhs, double h, double iters) {
    double h2 = h * h;

    for (int k = 0; k <= iters; ++k) {
        applyPressureNeumann(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p(i, j) = 0.25 * (p(i - 1, j) + p(i + 1, j)
                                + p(i, j - 1) + p(i, j + 1)
                                - h2 * rhs(i, j));
            }
        }
        applyPressureNeumann(p);
    }
}

// ------------------ PROJECTION ------------------
inline void project(FaceX& u, FaceY& v, CellField& p, CellField& div, double dt, double h, int iter) {
    //use divergence to calcuulate the rhs eq
    divergence(div, u, v, h);
    CellField rhs(Grid{div.nx, div.ny, h});
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            rhs(i, j) = div(i, j) / dt;
        }
    }
    //Solve with gaussSeidel
    gaussSeidelNeumann(p, rhs, h, iter);

    //substract the pressure gradient on the inner faces
    for (int i = 2; i <= u.nx - 1; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            u(i, j) -= dt * (p(i, j) - p(i - 1, j)) / h;
        }
    }
    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 2; j <= v.ny - 1; ++j) {
            v(i, j) -= dt * (p(i, j) - p(i, j - 1)) / h;
        }
    }

    //reimpose velocity BC on walls
    applyVelocityBC(u, v);
}

//Projection for channel
inline void projectChannel(FaceX& u, FaceY& v, CellField& p, CellField& div,
                           double dt, double h, double U_in, int iters = 80) {
    divergence(div, u, v, h);
    CellField rhs(Grid{div.nx, div.ny, h});
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            rhs(i, j) = div(i, j) / dt;
        }
    }

    // Gauss-Seidel con contorno de canal (Neumann en entrada/paredes, p=0 en salida)
    double h2 = h * h;
    for (int k = 0; k < iters; ++k) {
        applyChannelPressureBC(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p(i, j) = 0.25 * (p(i - 1, j) + p(i + 1, j)
                                + p(i, j - 1) + p(i, j + 1)
                                - h2 * rhs(i, j));
            }
        }
    }
    applyChannelPressureBC(p);

    // Restar gradiente en u: llegamos hasta i = u.nx (la salida es libre y p=0 la corrige)
    for (int i = 2; i <= u.nx; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            u(i, j) -= dt * (p(i, j) - p(i - 1, j)) / h;
        }
    }
    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 2; j <= v.ny - 1; ++j) {
            v(i, j) -= dt * (p(i, j) - p(i, j - 1)) / h;
        }
    }

    // Reimponer entrada, paredes y fantasmas de salida (sin machacar u(u.nx, j))
    applyChannelVelocityBC(u, v, U_in, false);
}

// ------------- Viscosity Diffusion ----------------

//Difussion for FaceX
inline void DiffuseU(FaceX& u, const FaceX& u0, double nu, double dt, double h, int iters) {
    if (nu <= 0.0) {
        u.a = u0.a;
        applyVelocityUBC(u);
        return;
    }
    double alpha = (nu * dt) / (h * h);
    double inv_denom = 1 / (1.0 + 4.0 * alpha);
    for (int k = 0; k <= iters; ++k) {
        for (int i = 2; i <= u.nx - 1; ++i) {
            for (int j = 1; j <= u.ny; ++j) {
                u(i, j) = (u0(i, j) + alpha * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1))) * inv_denom;
            }
        }
    }
    applyVelocityUBC(u);
}

//Difussion for FaceV
inline void DiffuseV(FaceY& v, const FaceY& v0, double nu, double dt, double h, int iters) {
    if (nu <= 0.0) {
        v.a = v0.a;
        applyVelocityVBC(v);
        return;
    }
    double alpha = (nu * dt) / (h * h);
    double inv_denom = 1 / (1.0 + 4.0 * alpha);
    for (int k = 0; k <= iters; ++k) {
        applyVelocityVBC(v);
        for (int i = 1; i <= v.nx; ++i) {
            for (int j = 2; j <= v.ny - 1; ++j) {
                v(i, j) = (v0(i, j) + alpha * (v(i + 1, j) + v(i - 1, j) + v(i, j + 1) + v(i, j - 1))) * inv_denom;
            }
        }
    }
    applyVelocityVBC(v);
}

//Diffusion for channel

inline void diffuseChannelU(FaceX& u, const FaceX& u0, const FaceY& /*v_dummy*/,
                            double nu, double dt, double h, double U_in, int iters = 30) {
    if (nu <= 0.0) { u.a = u0.a; return; }
    double alpha = (nu * dt) / (h * h);
    double inv_denom = 1.0 / (1.0 + 4.0 * alpha);

    for (int k = 0; k < iters; ++k) {
        for (int j = 1; j <= u.ny; ++j) {
            u(1, j)        = U_in;
            u(0, j)        = U_in;
            u(u.nx, j)     = u(u.nx - 1, j);
            u(u.nx + 1, j) = u(u.nx, j);
        }
        for (int i = 0; i <= u.nx + 1; ++i) {
            u(i, 0)        = -u(i, 1);
            u(i, u.ny + 1) = -u(i, u.ny);
        }

        for (int i = 2; i <= u.nx - 1; ++i) {
            for (int j = 1; j <= u.ny; ++j) {
                u(i, j) = (u0(i, j) + alpha * (u(i - 1, j) + u(i + 1, j)
                                             + u(i, j - 1) + u(i, j + 1))) * inv_denom;
            }
        }
    }
}

inline void diffuseChannelV(FaceY& v, const FaceY& v0,
                            double nu, double dt, double h, int iters = 30) {
    if (nu <= 0.0) { v.a = v0.a; return; }
    double alpha = (nu * dt) / (h * h);
    double inv_denom = 1.0 / (1.0 + 4.0 * alpha);

    for (int k = 0; k < iters; ++k) {
        for (int i = 1; i <= v.nx; ++i) {
            v(i, 1)    = 0.0;
            v(i, v.ny) = 0.0;
        }
        for (int j = 1; j <= v.ny; ++j) {
            v(0, j)        = -v(1, j);
            v(v.nx + 1, j) =  v(v.nx, j);
        }

        for (int i = 1; i <= v.nx; ++i) {
            for (int j = 2; j <= v.ny - 1; ++j) {
                v(i, j) = (v0(i, j) + alpha * (v(i - 1, j) + v(i + 1, j)
                                             + v(i, j - 1) + v(i, j + 1))) * inv_denom;
            }
        }
    }
}