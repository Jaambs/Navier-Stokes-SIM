#pragma once
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
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
    std::vector<double> a;      //nx+3 * ny + 2

    FaceX(const Grid& g) 
        : nx(g.nx + 1), ny(g.ny), stride(g.ny + 2), a((g.nx + 3) * (g.ny + 2), 0.0) {}

    double& operator()(int i, int j) {return a[stride * i + j];
    }
    double operator()(int i, int j) const {return a[stride * i + j];
    }
};

//Field for vertical velocity
struct FaceY {
    int nx, ny;                 //nx, ny + 1
    int stride;                 //ny + 3
    std::vector<double> a;      //nx + 2 * ny + 3

    FaceY(const Grid& g) 
        : nx(g.nx), ny(g.ny + 1), stride(g.ny + 3), a((g.nx + 2) * (g.ny + 3), 0.0) {}

    double& operator()(int i, int j) {return a[stride * i + j];
    }
    double operator()(int i, int j) const {return a[stride * i + j];
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

//Calculation of the velocity on vertical faces, horizontal faces and centre

//v evaluated on vertical face u(i,j)
inline double vAtU(const FaceY& v, int i, int j) {
    return 0.25 * (v(i-1,j)+v(i,j)+v(i-1,j+1) + v(i, j+1));
}
//u evaluated on horizontal face v(i,j)
inline double uAtV(const FaceX& u, int i, int j) {
    return 0.25 * (u(i,j-1)+u(i,j)+u(i+1,j-1) + u(i+1, j));
}
// u and v evaluated on the center of the cell (i,j)
inline double uAtCenter(const FaceX& u, int i, int j) {
    return 0.5 * (u(i,j) + u(i+1,j));
}
inline double vAtCenter (const FaceY& v, int i, int j) {
    return 0.5 * (v(i,j) + v(i,j+1));
}