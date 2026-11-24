#pragma once
#include <cstdint>
#include <vector>
#include "core/grid.hpp"

// Solid mask for the cylinder. A cell is solid if its center is inside the
// circle; a face is solid if either of its two neighbouring cells is solid.
struct CylinderMask {
    int nx, ny;
    bool active = true;
    Real xc_cyl, yc_cyl, diameter;
    std::vector<uint8_t> solid;          // (nx+2) * (ny+2)
    std::vector<uint8_t> solidU;         // (nx+3) * (ny+2)
    std::vector<uint8_t> solidV;         // (nx+2) * (ny+3)
    int strideC, strideU, strideV;

    CylinderMask(const Grid& g);

    uint8_t& isSolid(int i, int j)       {return solid[strideC * i + j]; }
    uint8_t  isSolid(int i, int j) const {return solid[strideC * i + j]; }

    uint8_t& isSolidU(int i, int j)       {return solidU[strideU * i + j]; }
    uint8_t  isSolidU(int i, int j) const {return solidU[strideU * i + j]; }

    uint8_t& isSolidV(int i, int j)       {return solidV[strideV * i + j]; }
    uint8_t  isSolidV(int i, int j) const {return solidV[strideV * i + j]; }
};

// u = 0, v = 0 on all cylinder faces (no-slip obstacle).
void applyCylinderVelocityBC(FaceX& u, FaceY& v, const CylinderMask& cyl);
