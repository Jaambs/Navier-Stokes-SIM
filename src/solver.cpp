#include <omp.h>
#include <algorithm>
#include <cmath>
#include "physics/solver.hpp"
#include "physics/BC.hpp"

// ----------- ADVECTION -------------

void advectU(FaceX& un, const FaceX& u, const FaceY& v, Real dt, Real h) {
    #pragma omp parallel for schedule(static)
    for (int i=1; i<=u.nx; ++i) {
        for (int j= 1; j <= u.ny; ++j) {
            //physical position and time step back
            Real x = xu(i,h) - dt * u(i,j);
            Real y = yu(j,h) - dt * vAtU(v, i, j);

            //convert x,y to grid index
            Real gx = x / h + 1.0f;
            Real gy = y / h + 0.5f;

            //Interpolacion bilineal
            un(i,j) = sample(u, gx, gy);
        }
    }
}

void advectV(FaceY& vn, const FaceX& u, const FaceY& v, Real dt, Real h) {
    #pragma omp parallel for schedule(static)
    for (int i=1; i<=v.nx; ++i) {
        for (int j= 1; j <= v.ny; ++j) {
            //physical position and time step back
            Real x = xv(i, h) - dt * uAtV(u, i, j);
            Real y = yv(j, h) - dt * v(i, j);

            //convert x,y to grid index
            Real gx = x / h + 0.5f;
            Real gy = y / h + 1.0f;

            //Interpolacion bilineal
            vn(i,j) = sample(v, gx, gy);
        }
    }
}

void advectScalar(CellField& d, const CellField& d0,
                  const FaceX& u, const FaceY& v, Real dt, Real h) {
    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= d.nx; ++i) {
        for (int j = 1; j <= d.ny; ++j) {
            Real x = xc(i, h) - dt * uAtCenter(u, i, j);
            Real y = yc(j, h) - dt * vAtCenter(v, i, j);

            // Convertir (x, y) a coordenadas de indice de centros
            Real gx = x / h + 0.5f;
            Real gy = y / h + 0.5f;

            d(i, j) = sample(d0, gx, gy);
        }
    }
}

// -------------- DIVERGENCE --------------

void divergence(CellField& div, const FaceX& u, const FaceY& v, Real h) {
    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            div(i, j) = (u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j)) / h;
        }
    }
}

Real maxAbs(const CellField& f) {
    Real max_val = 0.0f;
    for (int i = 1; i <= f.nx; ++i) {
        for (int j = 1; j <= f.ny; ++j) {
            max_val = std::max(max_val, std::abs(f(i, j)));
        }
    }
    return max_val;
}

// --------------- Poisson Solver -----------------

//Gauss-Seidel Red and Black for paralelitation
void gaussSeidelChannelRB(CellField& p, const CellField& rhs, const CylinderMask& cyl,
                          Real h, int iters) {
    Real h2 = h * h;

    #pragma omp parallel num_threads(8)
    {
        for (int k = 0; k < iters; ++k) {
            #pragma omp single
            applyChannelPressureBC(p);
            //color 0 red, color 1 black
            for (int color = 0; color < 2; ++color) {
                #pragma omp for schedule(static)
                for (int i = 1; i <= p.nx; ++i) {
                    int j_start = 1 + ((i + 1 + color) & 1);
                    for (int j = j_start; j <= p.ny; j += 2) {
                        if (cyl.isSolid(i,j)) continue;

                        Real p_c = p(i, j);
                        Real p_l = cyl.isSolid(i - 1, j) ? p_c : p(i - 1, j);
                        Real p_r = cyl.isSolid(i + 1, j) ? p_c : p(i + 1, j);
                        Real p_b = cyl.isSolid(i, j - 1) ? p_c : p(i, j - 1);
                        Real p_t = cyl.isSolid(i, j + 1) ? p_c : p(i, j + 1);

                        p(i, j) = 0.25f * (p_l + p_r + p_b + p_t - h2 * rhs(i, j));
                    }
                }
            }
        }
        applyChannelPressureBC(p);
    }
}

// ------------------ PROJECTION ------------------

void projectChannel(FaceX& u, FaceY& v, CellField& p, CellField& div, CylinderMask& cyl,
                    CellField& rhs, Real dt, Real h, Real U_in, Real nu, int iters) {

    divergence(div, u, v, h);
    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            rhs(i, j) = div(i, j) / dt;
        }
    }

    // Gauss-Seidel RB
    gaussSeidelChannelRB(p, rhs, cyl, h, iters);

    // substract gradient in u
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

    // Reimpose BC
    applyChannelVelocityBC(u, v, U_in, false, nu);
}

void projectChannelMG(FaceX& u, FaceY& v, CellField& p, CellField& div, const CylinderMask& cyl,
                      MultigridSolver& mg, Real dt, Real h, Real U_in, Real nu) {

    divergence(div, u, v, h);
    #pragma omp parallel for schedule(static)
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            mg.rhs[0](i, j) = div(i, j) / dt;
        }
    }

    mg.solve(p, 3);

    // substract gradient in u
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

    // Reimpose BC
    applyChannelVelocityBC(u, v, U_in, false, nu);
    applyCylinderVelocityBC(u, v, cyl);
}

// ------------- Viscosity Diffusion ----------------

void diffuseChannelU(FaceX& u, const FaceX& u0,
                     Real nu, Real dt, Real h, Real U_in, int iters) {
    if (nu <= 0.0f) { u.a = u0.a; return; }
    Real alpha = (nu * dt) / (h * h);
    Real inv_denom = 1.0f / (1.0f + 4.0f * alpha);

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

void diffuseChannelV(FaceY& v, const FaceY& v0,
                     Real nu, Real dt, Real h, int iters) {
    if (nu <= 0.0f) { v.a = v0.a; return; }
    Real alpha = (nu * dt) / (h * h);
    Real inv_denom = 1.0f / (1.0f + 4.0f * alpha);

    for (int k = 0; k < iters; ++k) {
        for (int i = 1; i <= v.nx; ++i) {
            v(i, 1)    = 0.0f;
            v(i, v.ny) = 0.0f;
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

// ------------------- DEPRECATED / UNUSED ----------------

void jacobi(CellField& p, const CellField& rhs, Real h, int iters) {
    CellField p_new = p;
    Real h2 = h*h;

    for (int k = 0; k < iters; ++k) {
        applyDirichletZero(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p_new(i,j) = 0.25f * (p(i - 1, j) + p(i + 1, j)
                                   + p(i, j - 1) + p(i, j + 1)
                                   - h2 * rhs(i, j));
            }
        }
        std::swap(p.a, p_new.a);
    }
    applyDirichletZero(p);
}

void gaussSeidelNeumann(CellField& p, const CellField& rhs, Real h, Real iters) {
    Real h2 = h * h;

    for (int k = 0; k <= iters; ++k) {
        applyPressureNeumann(p);
        for (int i = 1; i <= p.nx; ++i) {
            for (int j = 1; j <= p.ny; ++j) {
                p(i, j) = 0.25f * (p(i - 1, j) + p(i + 1, j)
                                + p(i, j - 1) + p(i, j + 1)
                                - h2 * rhs(i, j));
            }
        }
        applyPressureNeumann(p);
    }
}

void DiffuseU(FaceX& u, const FaceX& u0, Real nu, Real dt, Real h, int iters) {
    if (nu <= 0.0f) {
        u.a = u0.a;
        applyVelocityUBC(u);
        return;
    }
    Real alpha = (nu * dt) / (h * h);
    Real inv_denom = 1 / (1.0f + 4.0f * alpha);
    for (int k = 0; k <= iters; ++k) {
        for (int i = 2; i <= u.nx - 1; ++i) {
            for (int j = 1; j <= u.ny; ++j) {
                u(i, j) = (u0(i, j) + alpha * (u(i + 1, j) + u(i - 1, j)
                                             + u(i, j + 1) + u(i, j - 1))) * inv_denom;
            }
        }
    }
    applyVelocityUBC(u);
}

void DiffuseV(FaceY& v, const FaceY& v0, Real nu, Real dt, Real h, int iters) {
    if (nu <= 0.0f) {
        v.a = v0.a;
        applyVelocityVBC(v);
        return;
    }
    Real alpha = (nu * dt) / (h * h);
    Real inv_denom = 1 / (1.0f + 4.0f * alpha);
    for (int k = 0; k <= iters; ++k) {
        applyVelocityVBC(v);
        for (int i = 1; i <= v.nx; ++i) {
            for (int j = 2; j <= v.ny - 1; ++j) {
                v(i, j) = (v0(i, j) + alpha * (v(i + 1, j) + v(i - 1, j)
                                             + v(i, j + 1) + v(i, j - 1))) * inv_denom;
            }
        }
    }
    applyVelocityVBC(v);
}

void project(FaceX& u, FaceY& v, CellField& p, CellField& div, Real dt, Real h, int iter) {
    //use divergence to calcuulate the rhs eq
    divergence(div, u, v, h);
    CellField rhs_local(Grid{div.nx, div.ny, h});
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            rhs_local(i, j) = div(i, j) / dt;
        }
    }
    //Solve with gaussSeidel
    gaussSeidelNeumann(p, rhs_local, h, iter);

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
