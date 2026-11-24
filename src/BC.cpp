#include "physics/BC.hpp"

//----------------------- BOUNDARY CONDITIONS -----------------------------

void applyDirichletZero(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0,j)         = -p(1,j);
        p(p.nx+1, j)   = -p(p.nx, j);
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = -p(i, 1);
        p(i, p.ny + 1) = -p(i, p.ny);
    }
}

void applyPressureNeumann(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0,j)         = p(1,j);
        p(p.nx+1, j)   = p(p.nx, j);
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = p(i, 1);
        p(i, p.ny + 1) = p(i, p.ny);
    }
}

void applyVelocityBC(FaceX& u, FaceY& v) {
    //No trespassing meaning velocity = 0 on real walls
    for (int j = 1; j <= u.ny; ++j) {
        u(1, j)    = 0.0f;
        u(u.nx, j) = 0.0f;
    }
    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)    = 0.0f;
        v(i, v.ny) = 0.0f;
    }

    //No sliding condition, tangencial velocity on boundary = 0
    for (int i = 0; i <= u.nx; ++i) {
        u(i, 0)      = -u(i, 1);
        u(i, u.ny+1) = -u(i, u.ny);
    }
    for (int j = 0; j <= v.ny; ++j) {
        v(0, j)      = -v(1, j);
        v(v.nx+1, j) = -v(v.nx, j);
    }
}

void applyVelocityUBC(FaceX& u) {
    for (int j = 1; j <= u.ny; ++j) {
        u(1, j)    = 0.0f;
        u(u.nx, j) = 0.0f;
    }
    for (int i = 0; i <= u.nx; ++i) {
        u(i, 0)      = -u(i, 1);
        u(i, u.ny+1) = -u(i, u.ny);
    }
}

void applyVelocityVBC(FaceY& v) {
    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)    = 0.0f;
        v(i, v.ny) = 0.0f;
    }
    for (int j = 0; j <= v.ny; ++j) {
        v(0, j)      = -v(1, j);
        v(v.nx+1, j) = -v(v.nx, j);
    }
}

void applyChannelVelocityBC(FaceX& u, FaceY& v, Real U_in,
                            bool set_outlet_face, Real nu) {
    bool free_slip = (nu <= 0.0f);

    for (int j = 1; j <= u.ny; ++j) {
        u(1, j) = U_in;
        u(0, j) = U_in;
    }
    for (int j = 1; j <= v.ny; ++j) {
        v(0, j) = -v(1, j);
    }

    for (int j = 1; j <= u.ny; ++j) {
        if (set_outlet_face) {
            u(u.nx, j) = u(u.nx - 1, j);
        }
        u(u.nx + 1, j) = u(u.nx, j);
    }
    for (int j = 1; j <= v.ny; ++j) {
        v(v.nx + 1, j) = v(v.nx, j);
    }

    for (int i = 1; i <= v.nx; ++i) {
        v(i, 1)        = 0.0f;
        v(i, v.ny)     = 0.0f;
        v(i, 0)        = 0.0f;
        v(i, v.ny + 1) = 0.0f;
    }

    for (int i = 0; i <= u.nx + 1; ++i) {
        if (free_slip) {
            u(i, 0)       = u(i, 1);
            u(i, u.ny +1) = u(i, u.ny);
        } else {
            u(i, 0)        = -u(i, 1);
            u(i, u.ny + 1) = -u(i, u.ny);
        }
    }
}

void applyChannelPressureBC(CellField& p) {
    for (int j = 1; j <= p.ny; ++j) {
        p(0, j)        =  p(1, j);      //entrance dp/dx = 0
        p(p.nx + 1, j) = -p(p.nx, j);   //exit p = 0
    }
    for (int i = 1; i <= p.nx; ++i) {
        p(i, 0)        = p(i, 1);       //floor dp/dy = 0
        p(i, p.ny + 1) = p(i, p.ny);    //roof dp/dy = 0
    }
}
