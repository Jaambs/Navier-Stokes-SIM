#pragma once
#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

//Function to drop dye and apply forces when moving the mouse
inline void applyMouseBrush(FaceX& u, FaceY& v, CellField& dyeR, CellField& dyeB, CellField& dyeG,
                            Real mx, Real my, Real mdx, Real mdy,
                            Real r_color, Real b_color, Real g_color, Real dt, Real h, Real force_strengh) {
    Real radius = 3*h * force_strengh * 0.1f;
    Real r2 = radius * radius;

    //Calculate a box around the mouse so we dont have to check every single box
    int i_min = std::max(2,        static_cast<int>((mx - radius) / h));
    int i_max = std::min(u.nx - 1, static_cast<int>((mx + radius) / h) + 2);
    int j_min = std::max(1,        static_cast<int>((my - radius) / h));
    int j_max = std::min(u.ny,     static_cast<int>((my + radius) / h) + 2);

    //Add horizontal velocity on FaceX
    for (int i = i_min; i <= i_max; ++i) {
        for (int j = j_min; j <= j_max; ++j) {
            Real dx = xu(i, h) - mx;
            Real dy = yu(j, h) - my;
            Real d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                Real w = std::exp (Real(-d2/(0.5f * r2))); //smooth profile
                u(i, j) += dt * force_strengh * mdx * w;
            }
        }
    }

    //Redefine the limits
    i_min = std::max(1,        static_cast<int>((mx - radius) / h));
    i_max = std::min(v.nx,     static_cast<int>((mx + radius) / h) + 2);
    j_min = std::max(2,        static_cast<int>((my - radius) / h));
    j_max = std::min(v.ny - 1, static_cast<int>((my + radius) / h) + 2);
    //Add vertical velocity

    for (int i = i_min; i <= i_max; ++i) {
        for (int j = j_min; j <= j_max; ++j) {
            Real dx = xv(i, h) - mx;
            Real dy = yv(j, h) - my;
            Real d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                Real w = std::exp (Real(-d2/(0.5f * r2))); //smooth profile
                v(i, j) += dt * force_strengh * mdy * w;
            }
        }
    }

    //Redefine the limits
    i_min = std::max(1,        static_cast<int>((mx - radius) / h));
    i_max = std::min(dyeR.nx,     static_cast<int>((mx + radius) / h) + 2);
    j_min = std::max(1,        static_cast<int>((my - radius) / h));
    j_max = std::min(dyeR.ny, static_cast<int>((my + radius) / h) + 2);

    //Add dye on cell centers
    for (int i = i_min; i <= i_max; ++i) {
        for (int j = j_min; j <= j_max; ++j) {
            Real dx = xc(i, h) - mx;
            Real dy = yc(j, h) - my;
            Real d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                Real w = std::exp (Real(-d2 / (0.5f * r2)));
                dyeR(i, j) = std::min(1.0f, dyeR(i, j) + 0.5f * r_color * w);
                dyeB(i, j) = std::min(1.0f, dyeB(i, j) + 0.5f * b_color * w);
                dyeG(i, j) = std::min(1.0f, dyeG(i, j) + 0.5f * g_color * w);
            }
        }
    }

}