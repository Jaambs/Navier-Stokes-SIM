//Objective
// Build navier-stokes simulation using raylib for visualizing
// navier-stokes based on MAC with semi-lagragian advection Gauss-Seidel Diffusion and poisson / Gauss-Seidel presion
#include <iostream>
#include <omp.h>
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
#include "mouse.hpp"
#include "vorticity.hpp"
#include "BC.hpp"
#include "cylinder.hpp"
#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

int main(){

    // ---- SIM PARAMETERS ----

    const int NW         = 256;
    const int NH         = 128;
    const int SCALE      = 6;
    const int WIN_HEIGHT = NH * SCALE;
    const int WIN_WIDTH  = NW * SCALE;
    const double dt      = 1.0 / 60.0;
    double sim_time      = 0.0;
    SetTargetFPS(60);
    double nu = 1.67e-3;
    double force = 15.0;
    const double U_in = 1.0;
    const double L_ref = 1.0 / 6.0;
    Grid g{NW, NH, 1.0 / NH};
    CylinderMask cyl(g);
    bool inject_dye = true;
    bool use_advection = true;
    bool use_diffusion = true;
    bool use_projection = true;

    // Physical size of the domain: NW*h x NH*h = 2 x 1 for a 2:1 channel
    const double DOMAIN_W = NW * g.h;
    const double DOMAIN_H = NH * g.h;
     
    // ---- Velocity Divergence and Pressure Fields ----

    FaceX u(g), un(g);
    FaceY v(g), vn(g);
    CellField div(g), p(g), omega(g), rhs(g);
    
    // ---- RGB Dye Fields ----

    CellField dyeR(g), dyeR_new(g);
    CellField dyeG(g), dyeG_new(g);
    CellField dyeB(g), dyeB_new(g);
    enum class ViewMode {DYE, VORTICITY, VELOCITY};
    ViewMode view = ViewMode::DYE;

    // ---- Strohual Measurements ----

    int i_probe = static_cast<int>((cyl.xc_cyl + 3.0 * cyl.diameter) / g.h) + 1;
    int j_probe = g.ny / 2 + 1;
    //frecuency calculations and csv export variables
    double prev_v = 0.0;
    double last_cross_time = 1.0;
    double strohual = 0.0;

    // ---- Inicialization of the program ----

    InitWindow(WIN_WIDTH, WIN_HEIGHT, "Navier-Stokes SIM (MAC 256x128)");

    Image Img = GenImageColor(NW, NH, BLACK);
    Texture2D tex = LoadTextureFromImage(Img);
    SetTextureFilter(tex, TEXTURE_FILTER_BILINEAR);

    // ----- START OF THE SIMULATION LOOP -----
    while (!WindowShouldClose()) {
        sim_time += dt;

        //Reset the fluid with R
        if (IsKeyPressed(KEY_R)) {
        std::fill(u.a.begin(), u.a.end(), U_in);
        std::fill(v.a.begin(), v.a.end(), 0.0);
        std::fill(p.a.begin(), p.a.end(), 0.0);
        std::fill(dyeR.a.begin(), dyeR.a.end(), 0.0);
        std::fill(dyeG.a.begin(), dyeG.a.end(), 0.0);
        std::fill(dyeB.a.begin(), dyeB.a.end(), 0.0);
        prev_v = 0.0;
        last_cross_time = -1.0;
        strohual = 0.0;
        }

        //Switch Dye modes
        if (IsKeyPressed(KEY_TAB)) {
            if (view == ViewMode::DYE)              view = ViewMode::VORTICITY;
            else if (view == ViewMode::VORTICITY)   view = ViewMode::VELOCITY;
            else                                    view = ViewMode::DYE;
        }
        //Turn equation terms
        if (IsKeyPressed(KEY_ONE))   use_advection  = !use_advection;
        if (IsKeyPressed(KEY_TWO))   use_diffusion  = !use_diffusion;
        if (IsKeyPressed(KEY_THREE)) use_projection = !use_projection;

        //change viscosity
        if (IsKeyPressed(KEY_O)) {
            nu = (nu <= 0.0) ? 1e-5 : nu * 1.5;
        }
        if (IsKeyPressed(KEY_P)) {
            nu = (nu <= 1e-7) ? 0.0 : nu / 1.5;
        }
        if (IsKeyPressed(KEY_I)) {
            nu = 0.0;
        }

        //Change force strengh
        if (IsKeyPressed(KEY_K)) {
            force = (force <= 1e-2) ? 5.0 : force * 1.5;
        }
        if (IsKeyPressed(KEY_L)) {
            force = (force <= 1e-2) ? force : force / 1.5;
        }
        if (IsKeyPressed(KEY_J)) {
            force =15;
        }
        //Inyection of dye
        int cy = g.ny / 2;
        int offset = 3;
        int half_width = 1; 
        if (IsKeyPressed(KEY_D)) inject_dye = !inject_dye;
        for (int j = 1; j <= g.ny; ++j) {
            bool top_line = inject_dye && (std::abs(j - (cy + offset)) <= half_width);
            bool bot_line = inject_dye && (std::abs(j - (cy - offset)) <= half_width);

            if (top_line) {
                dyeR(1, j) = 0.0;
                dyeG(1, j) = 0.85;
                dyeB(1, j) = 1.0;
            } else if (bot_line) {
                // Línea inferior: Naranja
                dyeR(1, j) = 1.0;
                dyeG(1, j) = 0.3;
                dyeB(1, j) = 0.1;
            } else {
                // Fluido limpio en el resto de la entrada (y cuando apagas con D)
                dyeR(1, j) = 0.0;
                dyeG(1, j) = 0.0;
                dyeB(1, j) = 0.0;
            }
        }


        // For the simulation loop we follow the stablished Chorin schematic
        // Add forces -> advect -> diffuse -> project -> advect dye 
        //Then we render the result

        //1. ---Mouse Forces---
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 pos = GetMousePosition();
            Vector2 delta = GetMouseDelta();

            //Transform screen pixels to physic coordinates [0, 1] doing while so a Y axis inversion

            double mx  =  pos.x / WIN_WIDTH  * DOMAIN_W;
            double my  =  DOMAIN_H * (1.0 - pos.y / WIN_HEIGHT);
            double mdx =  delta.x / WIN_WIDTH  * DOMAIN_W * 60.0;
            double mdy = -delta.y / WIN_HEIGHT * DOMAIN_H * 60.0;

            //Make the applied colour change with time almost randomly using trigonometric smooth functions
            double r_col = 0.5 + 0.5 * std::sin(sim_time * 1.5);
            double b_col = 0.5 + 0.5 * std::sin(sim_time * 1.5 + 2.094);
            double g_col = 0.5 + 0.5 * std::sin(sim_time * 1.5 + 4.188);

            applyMouseBrush(u, v, dyeR, dyeB, dyeG, mx, my, mdx, mdy, r_col, b_col, g_col, dt, g.h, force);
            applyChannelVelocityBC(u, v, U_in, true, nu);
        }

        //2. --- Advect ---
        if (use_advection) {
            advectU(un, u, v, dt, g.h);
            advectV(vn, u, v, dt, g.h);
            std::swap(u.a, un.a);
            std::swap(v.a, vn.a);
            applyChannelVelocityBC(u, v, U_in, true, nu);
            applyCylinderVelocityBC(u, v, cyl);
        }

        //3. --- Diffuse ---
        if (use_diffusion) {
            std::swap(u.a, un.a);
            std::swap(v.a, vn.a);
            diffuseChannelU(u, un, nu, dt, g.h, U_in, 30);
            diffuseChannelV(v, vn, nu, dt, g.h, 30);
            applyCylinderVelocityBC(u, v, cyl);
        }
        // 4. --- Project ---
        if (use_projection) {
            projectChannel(u, v, p, div, cyl, rhs, dt, g.h, U_in, nu, 80);
            applyCylinderVelocityBC(u, v, cyl);
        }

        //strohual meditions
        double v_probe = v(i_probe, j_probe);

        if (prev_v < 0.0 && v_probe >= 0.0) {
            if (last_cross_time > 0.0) {
                double period = sim_time - last_cross_time;
                double freq = 1.0 / period;
                strohual = (freq * cyl.diameter) / U_in;
            }
            last_cross_time = sim_time;
        }
        prev_v = v_probe;

        // 5. --- Advect Dye ---

        advectScalar(dyeR_new, dyeR, u, v, dt, g.h);
        advectScalar(dyeG_new, dyeG, u, v, dt, g.h);
        advectScalar(dyeB_new, dyeB, u, v, dt, g.h);
        std::swap(dyeR_new, dyeR);
        std::swap(dyeG_new, dyeG);
        std::swap(dyeB_new, dyeB);

        // Blank dye inside the cylinder AFTER the swap, on the field we just computed
        for (int i = 1; i <= g.nx; ++i) {
            for (int j = 1; j <= g.ny; ++j) {
                if (cyl.isSolid(i,j)) {
                    dyeR(i, j) = 0.0;
                    dyeG(i, j) = 0.0;
                    dyeB(i, j) = 0.0;
                }
            }
        }
        
        //Compute vorticity if the view mode is Vorticity
        if (view == ViewMode::VORTICITY) {
            computeVorticity(omega, u, v, g.h);
        }
        // ----------- RENDERING -------------
        Color* pixels = static_cast<Color*>(Img.data);
        #pragma omp parallel for schedule(static)
        for (int py = 0; py < NH; ++py) {
            int j = NH - py;
            for (int px = 0; px < NW; ++px) {
                int i = px + 1;

                if (view == ViewMode::DYE){
                    unsigned char r_c = static_cast<unsigned char>(std::clamp(dyeR(i, j), 0.0, 1.0) * 255.0);
                    unsigned char g_c = static_cast<unsigned char>(std::clamp(dyeG(i, j), 0.0, 1.0) * 255.0);
                    unsigned char b_c = static_cast<unsigned char>(std::clamp(dyeB(i, j), 0.0, 1.0) * 255.0);
                    pixels[py * NW + px] = Color{r_c, g_c, b_c, 255};
                    if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};} 
                }
                else if (view == ViewMode::VORTICITY){
                    pixels[py * NW + px] = vorticityToColor(omega(i,j), 15.0);
                    if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};}
                }
                else {
                    double uc = uAtCenter(u, i, j);
                    double vc = vAtCenter(v, i, j);
                    double speed = std::hypot(uc, vc);
                    pixels[py * NW + px] = speedToColor(speed, 4);
                    if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};}
                }
            }
        }

        UpdateTexture(tex, Img.data);

        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(tex, Vector2{0.0f, 0.0f}, 0.0f, static_cast<float>(SCALE), WHITE);
            DrawFPS(10, 10);
            if (nu > 0.0) {
                double Re = (U_in * L_ref) / nu;
                DrawText(TextFormat("nu =%.2e | Re = %.1f | St = %.3f", nu, Re, strohual), 10, 35, 15, YELLOW);
            }
            if (nu == 0.0) {
                DrawText("nu = 0.0 | Re = inf", 10, 35, 15, YELLOW);
            }
            DrawText(TextFormat("force =%.3e", force), 10, 60, 15, YELLOW);
            DrawText(TextFormat("[1] Adv:%s  [2] Diff:%s  [3] Pres:%s",
                                use_advection  ? "ON" : "OFF",
                                use_diffusion  ? "ON" : "OFF",
                                use_projection ? "ON" : "OFF"),
                                10, 85, 15, ORANGE);
            DrawText("Controls: [R] Reset | [I,O,P] Viscosity | [J,K,L] Brush | [1,2,3] EQ Terms | [D] Dye Injection | [TAB] Change View", 10, 110, 15, LIGHTGRAY);
        EndDrawing();

        #if defined(PLATFORM_WEB)
            emscripten_sleep(0); 
        #endif
    }
    UnloadTexture(tex);
    UnloadImage(Img);
    CloseWindow();

    return 0;
}