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
    int N = 32;
    Grid g{N, N, 1.0 / N};
        
    //Test for Advection. Compare numerical to analitical solution
    //in this case we take the sol p(x,y) = sen(pix)sen(piy)
    //so rhs = -2pi*pi*sen(pix)sen(piy)
    //objective to obtain p(x,y through the jacobi and gaussSeidel solver)
    CellField rhs(g);
    CellField p_exact(g);

    //fill the cells with analitical sol and the rhs
    for (int i = 1; i <= g.nx; ++i) {
        for (int j = 1; j <= g.ny; ++j) {
            double x = xc(i, g.h);
            double y = yc(j, g.h);
            p_exact(i, j) = std::sin(PI*x)*std::sin(PI*y);
            rhs(i, j) = -2*PI*PI*std::sin(PI*x)*std::sin(PI*y);
        }
    }
    //build a temp function to calculate the max error between teh two cells
    auto maxError = [&](const CellField& p) {
        double err = 0;
        for (int i = 1; i <= g.nx; ++i) {
            for (int j = 1; j <= g.ny; ++j) {
                err=std::max(err,std::abs(p(i, j)-p_exact(i, j)));
            }
        }
        return err;
    };

    std::cout << std::setprecision(5) << std::scientific;
    std::cout << "Error máximo con respecto al comportamiento exacto\n";
    std::cout << "Iteraciones || Jacobi || Gauss-Seidel\n";

    //compare the results of different number of interactions for jacobi and gauss-seidel
    for (int iters : {50, 100, 200, 400, 800, 1600}) {
        CellField p_jac(g);
        CellField p_Seidel(g);

        jacobi(p_jac, rhs, g.h, iters);
        gaussSeidel(p_Seidel, rhs, g.h, iters);

        std::cout << std::setw(11) << iters << "||"
                  << std::setw(12) << maxError(p_jac)
                  << std::setw(12) << maxError(p_Seidel) << "\n";
    }



    return 0;
}