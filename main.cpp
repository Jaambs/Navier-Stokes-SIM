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

//Data structure MAC Framework
struct Grid {
    int nx, ny;                 //Number of cells on x and y
    double h;                   //size of the square cells
};

//Scalar field with cell centers (pressure, dye, div)
struct CellField {
    int nx, ny;                 //grid.nx, grid.ny
    int stride;                 //distance in memory between (i, j) and (i+1, j) used to convert 2D into 1D for a vector
    std::vector<double> a;      //(nx+2)*(ny+2)

    CellField(const Grid& g)
        : nx(g.nx), ny(g.ny), stride(g.ny +2), a((g.nx + 2)*(g.ny + 2), 0.0) {} 
    //Read-Write operator
    double& operator()(int i, int j) {
        assert(i>= 0 && i <= nx + 1 && j >= 0 && j<= ny+1);
        return a[stride *i +j];
    }

    //Read-Only operator
    double operator() (int i, int j) const {
        assert(i>= 0 && i <= nx + 1 && j >= 0 && j<= ny+1);
        return a[stride *i +j];
    }
};

//Field for horizontal velocity
struct FaceX {
    int nx, ny;                 //nx + 1, ny
    int stride;                 //ny+2
    std::vector<double> a;      //nx+3 * ny

    FaceX(const Grid& g) 
        : nx(g.nx + 1), ny(g.ny), stride(g.ny + 2), a((g.nx + 3) * (g.ny + 2), 0.0) {}

    double& operator()(int i, int j) {
        assert(i >= 0 && i <= nx + 1 && j >= 0 && j <= ny + 1);
        return a[stride * i + j];
    }
    double operator()(int i, int j) const {
        assert(i >= 0 && i <= nx + 1 && j >= 0 && j <= ny + 1);
        return a[stride * i + j];
    }
};

//Field for vertical velocity
struct FaceY {
    int nx, ny;                 //nx, ny + 1
    int stride;                 //ny + 3
    std::vector<double> a;      //nx + 2 * ny + 3

    FaceY(const Grid& g) 
        : nx(g.nx), ny(g.ny + 1), stride(g.ny + 3), a((g.nx + 2) * (g.ny + 3), 0.0) {}

    double& operator()(int i, int j) {
        assert(i >= 0 && i <= nx + 1 && j >= 0 && j <= ny + 1);
        return a[stride * i + j];
    }
    double operator()(int i, int j) const {
        assert(i >= 0 && i <= nx + 1 && j >= 0 && j <= ny + 1);
        return a[stride * i + j];
    }
};
//Center of the cells coordinates
inline double xc(int i, double h) {return (i-0.5)*h;}
inline double yc(int j, double h) {return (j-0.5)*h;}
//Center of the vertical faces coordinates
inline double xu(int i, double h) {return (i-1)*h;}
inline double yu(int j, double h) {return (j-0.5)*h;}
//Center of the horizontal faces coordinates
inline double xv(int i, double h) {return (i-0.5)*h;}
inline double yv(int j, double h) {return (j-1)*h;}


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
    return 0;
}