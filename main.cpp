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
    Grid g{10,5,1.0};
    CellField dye0(g), dye1(g);
    FaceX u(g), un(g);
    FaceY v(g), vn(g);

    double U = 1.0;
    double dt = 0.5;

    //uniform velocity including ghost cells
    for (int i = 0; i<= u.nx + 1; ++i) {
        for (int j = 0; j<= u.ny + 1; ++j)
            u(i, j) = U;
    }
    //Test 1 advectU de constante ha de ser constante
    advectU(un, u, v, dt, g.h);
    std::cout << "Test 1 (u const = 1.0 tras advectU): " << un(5,3) << "\n\n";

    //Test 2: Mancha de tinta en x = 3.5 y 2.5 -- celda(4,3)
    for (int i = 1; i <= g.nx; ++i) {
        for (int j = 1; j <= g.ny; ++j) {
            double dx = xc(i, g.h) - 3.5;
            double dy = yc(j, g.h) - 2.5;
            dye0(i, j) = std::exp(-(dx * dx + dy * dy));
        }
    }
    auto printRow = [&](const CellField& f, const char* label) {
        std::cout << label << ": ";
        for (int i = 1; i <= g.nx; ++i)
            std::cout << std::fixed << std::setprecision(2) << f(i, 3) << " ";
        std::cout << "\n";
    };

    printRow(dye0, "Fila central (t = 0.0)");

    //Test 3 Avance de 4 pasos en dt = 0.5 (desplazamiento de dos celdas a la derecha)
    for (int step = 1; step <=4; ++step) {
        advectScalar(dye1, dye0, u, v, dt, g.h);
        std::swap(dye0.a, dye1.a); //actualización de dye0 para el sig paso
        if (step % 2 == 0)
            printRow(dye0, "Tras 1.0s mas    ");
    }
    
    return 0;
}