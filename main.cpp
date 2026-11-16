//Objective
// Build navier-stokes simulation using raylib for visualizing
// navier-stokes based on MAC with semi-lagragian advection Gauss-Seidel Diffusion and poisson / Gauss-Seidel presion
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <raylib.h>
#include <raymath.h>
#include <cassert>
#include "grid.hpp"
#include "solver.hpp"

int main(){
    Grid g{16, 16, 1.0 / 16.0};
    FaceX u(g);
    FaceY v(g);
    CellField div(g);

    //Test 1: u = x, v = 0 -> div must be equal to 1
    for (int i = 1; i <= u.nx; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            u(i, j) = xu(i, g.h);
        }
    }

    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 1; j <= v.ny; ++j) {
            v(i, j) = 0.0;
        }
    }

    divergence(div, u, v, g.h);
    
    // Measure the error regarding 1
    double max_err_test1 = 0.0;
    for (int i = 1; i <= div.nx; ++i) {
        for (int j = 1; j <= div.ny; ++j) {
            max_err_test1 = std::max(max_err_test1, std::abs(div(i, j) - 1.0));
        }
    }

    // Printing the results
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Test 1 (u=x, v=0)  -> div(5,5) = " << div(5, 5)
              << " | Error max respecto a 1.0: " << max_err_test1 << "\n";

    // Test 2 u = y, v = -x -> div aprox 0
    for (int i = 1; i <= u.nx; ++i) {
        for (int j = 1; j <= u.ny; ++j) {
            u(i, j) = yu(j, g.h);
        }
    }

    for (int i = 1; i <= v.nx; ++i) {
        for (int j = 1; j <= v.ny; ++j) {
            v(i, j) = -xv(i, g.h);
        }
    }

    divergence(div,u,v,g.h);
    std::cout << "Test 2 (u=y, v=-x) -> maxAbs(div) = " << maxAbs(div) << "\n";
    
    return 0;
}