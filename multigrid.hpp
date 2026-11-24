#pragma once
#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include "grid.hpp"
#include "cylinder.hpp"

struct MultigridSolver {
    static constexpr int LEVELS = 3;

    Grid grids[LEVELS];
    std::vector<CellField> rhs;                 //Right-hand-side ec
    std::vector<CellField> sol;                 //Solution on the current level
    std::vector<CellField> res;                 //Residue r = rhs - Lap(sol)
    std::vector<std::vector<uint8_t>> solid;    //Solid mask on each level (cylinder mask rescaled)

    MultigridSolver(const Grid& g0, const CylinderMask& cyl) {

        //allocate the memory for each grid on each level
        for (int l = 0; l < LEVELS; ++l) {
            int factor = 1 << l;
            grids[l] = Grid{g0.nx / factor, g0.ny / factor, g0.h * Real(factor)};
            rhs.emplace_back(grids[l]);
            sol.emplace_back(grids[l]);
            res.emplace_back(grids[l]);
            solid.emplace_back((grids[l].nx + 2) * (grids[l].ny + 2), 0);
        }

        //Copy solid mask at level 0
        for (int i = 1; i <= grids[0].nx; ++i) {
            for (int j = 1; j <= grids[0].ny; ++j) {
                isSolid(0, i, j) = cyl.isSolid(i, j);
            }
        }

        //Propagate the solid mask to higher levels
        for (int l = 0; l < LEVELS - 1; ++l) {
            for (int ic = 1; ic <= grids[l + 1].nx; ++ic) {
                for (int jc = 1; jc <= grids[l + 1].ny; ++jc) {
                    int i0 = 2 * ic - 1, j0 = 2 * jc - 1;
                    uint8_t sum = isSolid(l, i0, j0) + isSolid(l, i0 + 1, j0)
                                + isSolid(l, i0, j0 + 1) + isSolid(l, i0 + 1, j0 +1);
                    isSolid(l + 1, ic, jc) = (sum >= 2) ? 1 : 0;
                }
            }
        }

    }
    uint8_t& isSolid(int l, int i, int j) {
        return solid[l][(grids[l].ny + 2) * i + j];
    }
    uint8_t isSolid(int l, int i, int j) const {
        return solid[l][(grids[l].ny + 2) * i + j];
    }

    //channel BC: Neumann floor, celing, entrance, dirichlet exit
    void applyBC(int l) {
        CellField& f = sol[l];
        for (int j = 1; j <= f.ny; ++j) {
            f(0, j)        = f(1, j);       //entrance dp/dx = 0
            f(f.nx + 1, j) = -f(f.nx, j);   //exit p = 0
        }
        for (int i = 1; i <= f.nx; ++i) {
            f(i, 0)       = f(i, 1);        //floor dp/dy = 0
            f(i, f.ny +1) = f(i, f.ny);     //roof dp/dy  = 0
        }
    }
    //Gauss-Seidel Smoother
    void smooth(int l, int iters) {
        CellField& f = sol[l];
        const CellField& b = rhs[l];
        Real h2 = grids[l].h * grids[l].h;

        for (int k = 0; k < iters; ++k) {
            applyBC(l);
            for (int i = 1; i <= f.nx; ++i) {
                for (int j = 1; j <= f.ny; ++j) {
                    if (isSolid(l, i, j)) continue;

                    Real p_c = f(i, j);
                    Real p_l = isSolid(l, i - 1, j) ? p_c : f(i - 1, j);
                    Real p_r = isSolid(l, i + 1, j) ? p_c : f(i + 1, j);
                    Real p_b = isSolid(l, i, j - 1) ? p_c : f(i, j - 1);
                    Real p_t = isSolid(l, i, j + 1) ? p_c : f(i, j + 1);

                    f(i, j) = Real(0.25) * (p_l + p_r + p_b + p_t - h2 * b(i, j));
                }
            }
        }
        applyBC(l);
    }

    //compute res = rhs - Lapl(sol)
    void computeResidue(int l) {
        const CellField& f = sol[l];
        const CellField& b = rhs[l];
        CellField& r = res[l];
        Real inv_h2 = Real(1.0) / (grids[l].h * grids[l].h);

        applyBC(l);
        for (int i = 1; i <= f.nx; ++i) {
            for (int j = 1; j <= f.ny; ++j) {
                if (isSolid(l, i, j)) {
                    r(i, j) = Real(0.0);
                    continue;
                }
                Real p_c = f(i, j);
                Real p_l = isSolid(l, i - 1, j) ? p_c : f(i - 1, j);
                Real p_r = isSolid(l, i + 1, j) ? p_c : f(i + 1, j);
                Real p_b = isSolid(l, i, j - 1) ? p_c : f(i, j - 1);
                Real p_t = isSolid(l, i, j + 1) ? p_c : f(i, j + 1);

                Real lap = (p_l + p_r + p_b + p_t - Real(4.0) * p_c) * inv_h2;
                r(i, j) = b(i, j) - lap;
            }
        }
    }

    //residue restriction
    void restrictResidue(int l) {
        const CellField& r_fine = res[l];
        CellField& b_coarse = rhs[l + 1];

        for (int ic = 1; ic <= grids[l + 1].nx; ++ic) {
            for (int jc = 1; jc <= grids [l + 1].ny; ++jc) {
                if (isSolid(l + 1, ic, jc)) {
                    b_coarse(ic, jc) = Real(0.0);
                    continue;
                }
                int i0 = 2 * ic - 1, j0 = 2 * jc - 1;
                Real sum_r = Real(0.0);
                int fluid_count = 0;

                if (!isSolid(l, i0, j0))         { sum_r += r_fine(i0, j0);         ++fluid_count; }
                if (!isSolid(l, i0 + 1, j0))     { sum_r += r_fine(i0 + 1, j0);     ++fluid_count; }
                if (!isSolid(l, i0, j0 + 1))     { sum_r += r_fine(i0, j0 + 1);     ++fluid_count; }
                if (!isSolid(l, i0 + 1, j0 + 1)) { sum_r += r_fine(i0 + 1, j0 + 1); ++fluid_count; }

                b_coarse(ic, jc) = (fluid_count > 0) ? (sum_r / Real(fluid_count)) : Real(0.0);
            }
        }
    }

    // Bilineal prologation of the cell (thick l+1 -> thin l)
    void prolongateAndAdd(int l) {
        applyBC(l + 1);
        const CellField& e_coarse = sol[l + 1];
        CellField& f_fine = sol[l];

        for (int ic = 1; ic <= grids[l + 1].nx; ++ic) {
            for (int jc = 1; jc <= grids[l + 1].ny; ++jc) {
                int i0 = 2* ic - 1, j0 = 2* jc - 1;

                for (int di = 0; di <= 1; ++di) {
                    for (int dj = 0; dj <= 1; ++dj) {
                        int ifine = i0 + di;
                        int jfine = j0 + dj;
                        if (isSolid(l, ifine, jfine)) continue;
                        
                        int in = (di == 0) ? ic - 1 : ic + 1;
                        int jn = (dj == 0) ? jc - 1 : jc + 1;

                        Real e00 = e_coarse(ic, jc);
                        Real e10 = isSolid(l + 1, in, jc) ? e00 : e_coarse(in, jc);
                        Real e01 = isSolid(l + 1, ic, jn) ? e00 : e_coarse(ic, jn);
                        Real e11 = isSolid(l + 1, in, jn) ? e00 : e_coarse(in, jn);

                        f_fine(ifine, jfine) += (Real(9.0) * e00 + Real(3.0) * e10 + 
                                    Real(3.0) * e01 + Real(1.0) * e11) / Real(16.0);
                    }
                }
            }
        }
    }
    //execute N V-cyles from the current warm start p state
    void solve(CellField& p, int v_cycles = 3, int nu1 = 3, int nu2 = 4, int coarse_iters = 20) {

        for (int c = 0; c < v_cycles; ++c) {
            //First we from gigger to thiner
            for (int l = 0; l < LEVELS - 1; ++l) {
                smooth(l, nu1);
                computeResidue(l);
                restrictResidue(l);
                std::fill(sol[l + 1].a.begin(), sol[l + 1].a.end(), Real(0.0));
            }

            //Now we finish on the second level
            smooth(LEVELS - 1, coarse_iters);

            //And now we go up again
            for (int l = LEVELS - 2; l >= 0; --l) {
                prolongateAndAdd(l);
                smooth(l, nu2);
            }
        }
        p.a = sol[0].a;
    }
};