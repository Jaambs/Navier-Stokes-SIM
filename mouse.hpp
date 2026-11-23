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
                            double mx, double my, double mdx, double mdy,
                            double r_color, double b_color, double g_color, double dt, double h, double force_strengh) {
    double radius = 3*h * force_strengh * 0.1;
    double r2 = radius * radius;

    //Calculate a box around the mouse so we dont have to check every single box
    int i_min = std::max(2,        static_cast<int>((mx - radius) / h));
    int i_max = std::min(u.nx - 1, static_cast<int>((mx + radius) / h) + 2);
    int j_min = std::max(1,        static_cast<int>((my - radius) / h));
    int j_max = std::min(u.ny,     static_cast<int>((my + radius) / h) + 2);

    //Add horizontal velocity on FaceX
    for (int i = i_min; i <= i_max; ++i) {
        for (int j = j_min; j <= j_max; ++j) {
            double dx = xu(i, h) - mx;
            double dy = yu(j, h) - my;
            double d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                double w = std::exp (-d2/(0.5 * r2)); //smooth profile
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
            double dx = xv(i, h) - mx;
            double dy = yv(j, h) - my;
            double d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                double w = std::exp (-d2/(0.5 * r2)); //smooth profile
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
            double dx = xc(i, h) - mx;
            double dy = yc(j, h) - my;
            double d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                double w = std::exp (-d2 / (0.5 * r2));
                dyeR(i, j) = std::min(1.0, dyeR(i, j) + 0.5 * r_color * w);
                dyeB(i, j) = std::min(1.0, dyeB(i, j) + 0.5 * b_color * w);
                dyeG(i, j) = std::min(1.0, dyeG(i, j) + 0.5 * g_color * w);
            }
        }
    }

}