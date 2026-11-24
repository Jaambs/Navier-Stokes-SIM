#pragma once
#include <omp.h>
#include <vector>
#include <cstdint>
#include "grid.hpp"

// We have to build a structure that impose the condition u = 0, v = 0, whenever the corresponding face is from the cylinder
// in order to impose the no-slip condition

struct CylinderMask {
    int nx, ny;
    bool active = true;
    Real xc_cyl, yc_cyl, diameter;
    std::vector<uint8_t> solid;          // (nx+2) * (ny+2)
    std::vector<uint8_t> solidU;         // (nx+3) * (ny+2)
    std::vector<uint8_t> solidV;         // (nx+2) * (ny+3)
    int strideC, strideU, strideV;

    CylinderMask(const Grid& g)
        : nx(g.nx), ny(g.ny),
        solid((g.nx + 2) * (g.ny + 2), 0),
        solidU((g.nx + 3) * (g.ny + 2), 0),
        solidV((g.nx + 2) * (g.ny + 3), 0),
        strideC(g.ny + 2), strideU(g.ny + 2), strideV(g.ny + 3)
    {
        // Cylinder dimensions
        diameter = (g.ny / 6.0f) * g.h;
        xc_cyl   = (g.nx / 5.0f) * g.h;
        yc_cyl   = (g.ny / 2.0f) * g.h;
        active   = (diameter > Real(0.0));
        if (!active) return;
        Real r2 = (0.5f * diameter) * (0.5f * diameter);

        // Marking cells that should be solid
        // Centers
        for (int i = 1; i <= g.nx; ++i) {
            for (int j = 1; j <= g.ny; ++j) {
                Real dx = xc(i, g.h) - xc_cyl;
                Real dy = yc(j, g.h) - yc_cyl;
                if (dx * dx + dy * dy < r2) {
                    isSolid(i, j) = 1;
                }
            }
        }
        // FaceX
        for (int i = 1; i <= g.nx + 1; ++i) {
            for (int j = 1; j <= g.ny; ++j) {
                if (isSolid(i - 1, j) || isSolid(i, j)) {
                    isSolidU(i, j) = 1;
                }
            }
        }
        // FaceY
        for (int i = 1; i <= g.nx; ++i) {
            for (int j = 1; j <= g.ny + 1; ++j) {
                if (isSolid(i, j - 1) || isSolid(i, j)) {
                    isSolidV(i, j) = 1;
                }
            }
        }
    }

    uint8_t& isSolid(int i, int j)       {return solid[strideC * i + j]; }
    uint8_t  isSolid(int i, int j) const {return solid[strideC * i + j]; }

    uint8_t& isSolidU(int i, int j)       {return solidU[strideU * i + j]; }
    uint8_t  isSolidU(int i, int j) const {return solidU[strideU * i + j]; }

    uint8_t& isSolidV(int i, int j)       {return solidV[strideV * i + j]; }
    uint8_t  isSolidV(int i, int j) const {return solidV[strideV * i + j]; }
};

// u = 0, v = 0 on all cylinder faces

inline void applyCylinderVelocityBC(FaceX& u, FaceY& v, const CylinderMask& cyl) {
    if (!cyl.active) return;
    for (int i = 1; i <= u.nx; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            if (cyl.isSolidU(i, j)) u(i, j) = 0.0f;
        }
    }

    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 1; j <= v.ny; ++j) {
            if (cyl.isSolidV(i, j)) v(i, j) = 0.0f;
        }
    }
}