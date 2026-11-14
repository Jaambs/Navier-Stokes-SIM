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
    Grid g{4,4,1.0};
    CellField p(g);
    FaceX u(g);
    FaceY v(g);

    for (int i = 1; i <= g.nx; ++i)
        for (int j = 1; j <= g.ny; ++j)
            p(i, j) = xc(i, g.h) + yc(j, g.h);

    for (int i = 1; i <= u.nx; ++i) 
        for (int j = 1; j <= u.ny; ++j)
            u(i, j) = xu(i, g.h);

    for (int i = 1; i <= v.nx; ++i)
        for (int j = 1; j <= v.ny; ++j) 
            v(i, j) = yv(j, g.h);

    //Now we print u(i,j), the results should be 0.0,1.0,2.0,3.0,4.0
    std::cout << "--- Malla u (xu) ---\n";
    for (int j = u.ny; j >= 1; --j) {
        for (int i = 1; i <= u.nx; ++i)
            std::cout << std::setw(5) << u(i, j) << " ";
        std::cout << "\n";
    }
    std::cout << "--- Malla v (xv) ---\n";
    for (int j = v.ny; j>=1; --j) {
        for (int i = 1; i <=v.nx; ++i)
            std::cout << std::setw(5) << v(i,j) << " ";
        std::cout <<"\n";
    }
    std::cout << "--- Malla p ---\n";
    for (int j = p.ny; j>=1; --j) {
        for (int i = 1; i <=p.nx; ++i)
            std::cout << std::setw(5) << std::fixed << std::setprecision(2) << p(i,j) << " ";
        std::cout <<"\n";
    }
    return 0;
}