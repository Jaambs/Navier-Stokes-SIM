#include <algorithm>
#include <cmath>
#include "physics/vorticity.hpp"

void computeVorticity(CellField& omega, const FaceX& u, const FaceY& v, Real h) {
    for (int i = 1; i <= omega.nx; ++i) {
        for (int j = 1; j <= omega.ny; ++j) {
            Real w00 = cornerVorticity(u, v, i,     j,     h);
            Real w10 = cornerVorticity(u, v, i + 1, j,     h);
            Real w01 = cornerVorticity(u, v, i,     j + 1, h);
            Real w11 = cornerVorticity(u, v, i + 1, j + 1, h);
            omega(i, j) = 0.25f * (w00 + w01 + w10 + w11);
        }
    }
}

//Color map for the vorticity Blue -> negative, Red -> Positive, White -> 0
Color vorticityToColor(Real w, Real max_w) {
    Real t = std::clamp(w / max_w, -1.0f, 1.0f);
    if (t >= 0) { //Blanco (255, 255, 255) a rojo (220, 30, 30)
        unsigned char gb = static_cast<unsigned char>(255 * (1 - 0.88f * t));
        unsigned char r = static_cast<unsigned char>(255 - 35.0f * t);
        return Color{r, gb, gb, 255};
    } else { //Blanco a azul (30, 70, 230)
        Real a = -t;
        unsigned char r = static_cast<unsigned char>(255 * (1.0f - 0.88f * a));
        unsigned char g = static_cast<unsigned char>(255 * (1.0f - 0.72f * a));
        unsigned char b = static_cast<unsigned char>(255 - 25.0f * a);
        return Color {r, g, b, 255};
    }
}

//Color map for speed: black -> cyan -> yellow -> red
Color speedToColor(Real speed, Real max_speed) {
    Real t = std::clamp(speed / max_speed, 0.0f, 1.0f);
    Real r = 0.0f, g = 0.0f, b = 0.0f;

    if (t < 0.25f) {
        // Azul oscuro (0, 0, 0.5) -> Cian (0, 1, 1)
        Real u = t / 0.25f;
        r = 0.0f;
        g = u;
        b = 0.5f + 0.5f * u;
    } else if (t < 0.5f) {
        // Cian (0, 1, 1) -> Verde (0, 1, 0)
        Real u = (t - 0.25f) / 0.25f;
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - u;
    } else if (t < 0.75f) {
        // Verde (0, 1, 0) -> Amarillo (1, 1, 0)
        Real u = (t - 0.5f) / 0.25f;
        r = u;
        g = 1.0f;
        b = 0.0f;
    } else {
        // Amarillo (1, 1, 0) -> Rojo (0.9, 0, 0)
        Real u = (t - 0.75f) / 0.25f;
        r = 1.0f - 0.1f * u;
        g = 1.0f - u;
        b = 0.0f;
    }

    return Color{
        static_cast<unsigned char>(r * 255.0f),
        static_cast<unsigned char>(g * 255.0f),
        static_cast<unsigned char>(b * 255.0f),
        255
    };
}
