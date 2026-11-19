#include <algorithm>
#include <cmath>
#include "solver.hpp"
#include "grid.hpp"

//exact vorticity on a cell corner
inline double cornerVorticity (const FaceX& u, const FaceY& v, int i, int j, double h) {
    return ((v(i, j) - v(i - 1, j)) - (u(i, j) - u(i, j - 1))) / h;
}

//Vorticity calculation on center 
inline void computeVorticity(CellField& omega, const FaceX& u, const FaceY& v, double h) {
    for (int i = 1; i <= omega.nx; ++i) {
        for (int j = 1; j <= omega.ny; ++j) {
            double w00 = cornerVorticity(u, v, i,     j,     h);
            double w10 = cornerVorticity(u, v, i + 1, j,     h);
            double w01 = cornerVorticity(u, v, i,     j + 1, h);
            double w11 = cornerVorticity(u, v, i + 1, j + 1, h);
            omega(i, j) = 0.25 * (w00 + w01 + w10 + w11);
        }
    }
}

//Color map for the vorticity Blue -> negative, Red -> Positive, White -> 0
inline Color vorticityToColor(double w, double max_w) {
    double t = std::clamp(w / max_w, -1.0, 1.0);
    if (t >= 0) { //Blanco (255, 255, 255) a rojo (220, 30, 30)
        unsigned char gb = static_cast<unsigned char>(255 * (1 - 0.88 * t));
        unsigned char r = static_cast<unsigned char>(255 - 35.0 * t);
        return Color{r, gb, gb, 255};
    } else { //Blanco a azul (30, 70, 230)
        double a = -t;
        unsigned char r = static_cast<unsigned char>(255 * (1.0 - 0.88 * a));
        unsigned char g = static_cast<unsigned char>(255 * (1.0 - 0.72 * a));
        unsigned char b = static_cast<unsigned char>(255 - 25.0 * a);
        return Color {r, g, b, 255};
    }
}

//Color map for speed: black -> cyan -> yellow -> white
inline Color speedToColor(double speed, double max_speed) {
    double t = std::clamp(speed / max_speed, 0.0, 1.0);
    double r = 0.0, g = 0.0, b = 0.0;

    if (t < 0.25) {
        // Azul oscuro (0, 0, 0.5) -> Cian (0, 1, 1)
        double u = t / 0.25;
        r = 0.0;
        g = u;
        b = 0.5 + 0.5 * u;
    } else if (t < 0.5) {
        // Cian (0, 1, 1) -> Verde (0, 1, 0)
        double u = (t - 0.25) / 0.25;
        r = 0.0;
        g = 1.0;
        b = 1.0 - u;
    } else if (t < 0.75) {
        // Verde (0, 1, 0) -> Amarillo (1, 1, 0)
        double u = (t - 0.5) / 0.25;
        r = u;
        g = 1.0;
        b = 0.0;
    } else {
        // Amarillo (1, 1, 0) -> Rojo (0.9, 0, 0)
        double u = (t - 0.75) / 0.25;
        r = 1.0 - 0.1 * u;
        g = 1.0 - u;
        b = 0.0;
    }

    return Color{
        static_cast<unsigned char>(r * 255.0),
        static_cast<unsigned char>(g * 255.0),
        static_cast<unsigned char>(b * 255.0),
        255
    };
}
