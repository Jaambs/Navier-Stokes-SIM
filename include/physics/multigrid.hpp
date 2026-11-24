#pragma once
#include <cassert>
#include <cstdint>
#include <vector>
#include "core/grid.hpp"
#include "physics/cylinder.hpp"

// Geometric multigrid, V-cycle, for the pressure Poisson equation on a uniform
// cell-centered grid with Neumann walls and p = 0 at the outlet.
struct MultigridSolver {
    static constexpr int LEVELS = 3;

    Grid grids[LEVELS];
    std::vector<CellField> rhs;                 //Right-hand-side on each level
    std::vector<CellField> sol;                 //Solution on each level
    std::vector<CellField> res;                 //Residue r = rhs - Lap(sol)
    std::vector<std::vector<uint8_t>> solid;    //Solid mask on each level (cylinder rescaled)

    MultigridSolver(const Grid& g0, const CylinderMask& cyl);

    uint8_t& isSolid(int l, int i, int j) {
        assert(l >= 0 && l < LEVELS);
        assert(i >= 0 && i <= grids[l].nx + 1 && j >= 0 && j <= grids[l].ny + 1);
        return solid[l][(grids[l].ny + 2) * i + j];
    }
    uint8_t isSolid(int l, int i, int j) const {
        assert(l >= 0 && l < LEVELS);
        assert(i >= 0 && i <= grids[l].nx + 1 && j >= 0 && j <= grids[l].ny + 1);
        return solid[l][(grids[l].ny + 2) * i + j];
    }

    void applyBC(int l);                 //channel BC on level l
    void smooth(int l, int iters);       //Gauss-Seidel smoother
    void computeResidue(int l);          //res = rhs - Laplacian(sol)
    void restrictResidue(int l);         //restrict res[l] into rhs[l+1]
    void prolongateAndAdd(int l);        //prolongate sol[l+1] into sol[l]
    void solve(CellField& p, int v_cycles = 3, int nu1 = 3, int nu2 = 4, int coarse_iters = 20);
    void updateMask(const CylinderMask& cyl);
};
