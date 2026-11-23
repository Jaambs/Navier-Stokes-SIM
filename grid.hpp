#pragma once
#include <omp.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>


using Real = float;
//Data structure MAC Framework
struct Grid {
    int nx, ny;                 //Number of cells on x and y
    Real h;                   //size of the square cells
};

//Scalar field with cell centers (pressure, dye, div)
struct CellField {
    int nx, ny;                 //grid.nx, grid.ny
    int stride;                 //distance in memory between (i, j) and (i+1, j) used to convert 2D into 1D for a vector
    std::vector<Real> a;      //(nx+2)*(ny+2)

    CellField(const Grid& g)
        : nx(g.nx), ny(g.ny), stride(g.ny +2), a((g.nx + 2)*(g.ny + 2), 0.0f) {} 
    //Read-Write operator
    Real& operator()(int i, int j) {
        assert(i>= 0 && i <= nx + 1 && j >= 0 && j<= ny+1);
        return a[stride *i +j];
    }

    //Read-Only operator
    Real operator() (int i, int j) const {
        assert(i>= 0 && i <= nx + 1 && j >= 0 && j<= ny+1);
        return a[stride *i +j];
    }
};

//Field for horizontal velocity
struct FaceX {
    int nx, ny;                 //nx + 1, ny
    int stride;                 //ny+2
    std::vector<Real> a;      //nx+3 * ny + 2

    FaceX(const Grid& g) 
        : nx(g.nx + 1), ny(g.ny), stride(g.ny + 2), a((g.nx + 3) * (g.ny + 2), 0.0f) {}

    Real& operator()(int i, int j) {return a[stride * i + j];
    }
    Real operator()(int i, int j) const {return a[stride * i + j];
    }
};

//Field for vertical velocity
struct FaceY {
    int nx, ny;                 //nx, ny + 1
    int stride;                 //ny + 3
    std::vector<Real> a;      //nx + 2 * ny + 3

    FaceY(const Grid& g) 
        : nx(g.nx), ny(g.ny + 1), stride(g.ny + 3), a((g.nx + 2) * (g.ny + 3), 0.0f) {}

    Real& operator()(int i, int j) {return a[stride * i + j];
    }
    Real operator()(int i, int j) const {return a[stride * i + j];
    }
};

//Center of the cells coordinates
inline Real xc(int i, Real h) {return (i-0.5f)*h;}
inline Real yc(int j, Real h) {return (j-0.5f)*h;}
//Center of the vertical faces coordinates
inline Real xu(int i, Real h) {return (i-1)*h;}
inline Real yu(int j, Real h) {return (j-0.5f)*h;}
//Center of the horizontal faces coordinates
inline Real xv(int i, Real h) {return (i-0.5f)*h;}
inline Real yv(int j, Real h) {return (j-1)*h;}

//Calculation of the velocity on vertical faces, horizontal faces and centre

//v evaluated on vertical face u(i,j)
inline Real vAtU(const FaceY& v, int i, int j) {
    return 0.25f * (v(i-1,j)+v(i,j)+v(i-1,j+1) + v(i, j+1));
}
//u evaluated on horizontal face v(i,j)
inline Real uAtV(const FaceX& u, int i, int j) {
    return 0.25f * (u(i,j-1)+u(i,j)+u(i+1,j-1) + u(i+1, j));
}
// u and v evaluated on the center of the cell (i,j)
inline Real uAtCenter(const FaceX& u, int i, int j) {
    return 0.5f * (u(i,j) + u(i+1,j));
}
inline Real vAtCenter (const FaceY& v, int i, int j) {
    return 0.5f * (v(i,j) + v(i,j+1));
}