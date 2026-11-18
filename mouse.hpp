#pragma once
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
                            double r_color, double b_color, double g_color, double dt, double h) {
    double radius = 3*h;
    double r2 = radius * radius;
    double force_strengh = 15.0;

    //Add horizontal velocity on the faces FaceX
    //We are going to check every single box, on an small grid it isnt really necesary create bounding parameters
    for (int i = 2; i <= u.nx-1; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            double dx = xu(i, h) - mx;
            double dy = yu(j, h) - my;
            double d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                double w = std::exp (-d2/(0.5 * r2)); //smooth profile
                u(i, j) += dt * force_strengh * mdx * w;
            }
        }
    }
    
    //Add vertical velocity

    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 2; j <= v.ny-1; ++j) {
            double dx = xv(i, h) - mx;
            double dy = yv(j, h) - my;
            double d2 = dx * dx + dy * dy;
            if (d2 < r2) {
                double w = std::exp (-d2/(0.5 * r2)); //smooth profile
                v(i, j) += dt * force_strengh * mdy * w;
            }
        }
    }

    //Add dye on cell centers
    for (int i = 1; i <= dyeR.nx; ++i) {
        for (int j = 1; j <= dyeR.ny; ++j) {
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