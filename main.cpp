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
#include "mouse.hpp"
#include "vorticity.hpp"
#include "BC.hpp"

int main(){

    // ---- SIM PARAMETERS ----

    const int NW         = 256;
    const int NH         = 128;
    const int SCALE      = 6;
    const int WIN_HEIGHT = NH * SCALE;
    const int WIN_WIDTH  = NW * SCALE;
    const double dt      = 1.0 / 60.0;
    double sim_time      = 0.0;
    double nu = 1e-4;
    double force = 15.0;
    const double U_ref = 1.0;
    const double L_ref = 1.0;
    Grid g{NW, NH, 1.0 / NH};

    // Physical size of the domain: NW*h x NH*h = 2 x 1 for a 2:1 channel
    const double DOMAIN_W = NW * g.h;
    const double DOMAIN_H = NH * g.h;
     
    // ---- Velocity Divergence and Pressure Fields ----

    FaceX u(g), un(g);
    FaceY v(g), vn(g);
    CellField div(g), p(g), omega(g);
    
    // ---- RGB Dye Fields ----

    CellField dyeR(g), dyeR_new(g);
    CellField dyeG(g), dyeG_new(g);
    CellField dyeB(g), dyeB_new(g);
    enum class ViewMode {DYE, VORTICITY, VELOCITY};
    ViewMode view = ViewMode::DYE;

    // ---- Inicialization of the program ----

    InitWindow(WIN_WIDTH, WIN_HEIGHT, "Navier-Stokes SIM (MAC 256x128)");
    SetTargetFPS(60);

    Image Img = GenImageColor(NW, NH, BLACK);
    Texture2D tex = LoadTextureFromImage(Img);
    SetTextureFilter(tex, TEXTURE_FILTER_BILINEAR);

    // ----- START OF THE SIMULATION LOOP -----

    while(!WindowShouldClose()) {
        sim_time += dt;

        //Reset the fluid with R
        if (IsKeyPressed(KEY_R)) {
        std::fill(u.a.begin(), u.a.end(), 0.0);
        std::fill(v.a.begin(), v.a.end(), 0.0);
        std::fill(p.a.begin(), p.a.end(), 0.0);
        std::fill(dyeR.a.begin(), dyeR.a.end(), 0.0);
        std::fill(dyeG.a.begin(), dyeG.a.end(), 0.0);
        std::fill(dyeB.a.begin(), dyeB.a.end(), 0.0);
        }

        //Switch Dye modes
        if (IsKeyPressed(KEY_TAB)) {
            if (view == ViewMode::DYE)              view = ViewMode::VORTICITY;
            else if (view == ViewMode::VORTICITY)   view = ViewMode::VELOCITY;
            else                                    view = ViewMode::DYE;
        }

        //change viscosity
        if (IsKeyPressed(KEY_EQUAL) || IsKeyPressed(KEY_KP_ADD)) {
            nu = (nu == 0) ? 1e-5 : nu *= 1.5;
        }
        if (IsKeyPressed(KEY_MINUS) || IsKeyPressed(KEY_KP_SUBTRACT)) {
            nu = (nu <= 1e-7) ? 0 : nu /= 1.5;
        }
        if (IsKeyPressed(KEY_ZERO) || IsKeyPressed(KEY_KP_0)) {
            nu = 0.0;
        }

        //Change force strengh
        if (IsKeyPressed(KEY_PERIOD)) {
            force = (force <= 1e-2) ? 5 : force *= 1.5;
        }
        if (IsKeyPressed(KEY_COMMA)) {
            force = (force <= 1e-2) ? force : force /= 1.5;
        }
        if (IsKeyPressed(KEY_M)) {
            force =15;
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
            applyVelocityBC(u, v);
        }

        //2. --- Advect ---
        advectU(un, u, v, dt, g.h);
        advectV(vn, u, v, dt, g.h);
        std::swap(u.a, un.a);
        std::swap(v.a, vn.a);
        applyVelocityBC(u, v);

        //3. --- Diffuse ---
        std::swap(u.a, un.a);
        std::swap(v.a, vn.a);
        DiffuseU(u, un, nu, dt, g.h, 25);
        DiffuseV(v, vn, nu, dt, g.h, 25);

        // 4. --- Project ---
        project(u, v, p, div, dt, g.h, 10);

        // 5. --- Advect Dye ---

        advectScalar(dyeR_new, dyeR, u, v, dt, g.h);
        advectScalar(dyeG_new, dyeG, u, v, dt, g.h);
        advectScalar(dyeB_new, dyeB, u, v, dt, g.h);
        std::swap(dyeR_new, dyeR);
        std::swap(dyeG_new, dyeG);
        std::swap(dyeB_new, dyeB);
        
        //Compute vorticity if the view mode is Vorticity
        if (view == ViewMode::VORTICITY) {
            computeVorticity(omega, u, v, g.h);
        }
        // ----------- RENDERING -------------
        Color* pixels = static_cast<Color*>(Img.data);
        for (int py = 0; py < NH; ++py) {
            int j = NH - py;
            for (int px = 0; px < NW; ++px) {
                int i = px + 1;

                if (view == ViewMode::DYE){
                    unsigned char r_c = static_cast<unsigned char>(std::clamp(dyeR(i, j), 0.0, 1.0) * 255.0);
                    unsigned char g_c = static_cast<unsigned char>(std::clamp(dyeG(i, j), 0.0, 1.0) * 255.0);
                    unsigned char b_c = static_cast<unsigned char>(std::clamp(dyeB(i, j), 0.0, 1.0) * 255.0);
                    pixels[py * NW + px] = Color{r_c, g_c, b_c, 255};   
                }
                else if (view == ViewMode::VORTICITY){
                    pixels[py * NW + px] = vorticityToColor(omega(i,j), 15.0);
                }
                else {
                    double uc = uAtCenter(u, i, j);
                    double vc = vAtCenter(v, i, j);
                    double speed = std::hypot(uc, vc);
                    pixels[py * NW + px] = speedToColor(speed, 0.8);
                }
            }
        }
        UpdateTexture(tex, Img.data);

        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(tex, Vector2{0.0f, 0.0f}, 0.0f, static_cast<float>(SCALE), WHITE);
            DrawFPS(10, 10);
            DrawText("Press [R] to reset", 10, 35, 18, LIGHTGRAY);
            if (nu > 0.0) {
                double Re = (U_ref * L_ref) / nu;
                DrawText(TextFormat("nu =%.2e | Re = %.1f", nu, Re), 10, 60, 20, YELLOW);
            }
            if (nu == 0.0) {
                DrawText("nu = 0.0 | Re = inf", 10, 60, 20, YELLOW);
            }
            DrawText(TextFormat("force =%.3e", force), 10, 85, 20, YELLOW);
        EndDrawing();
    }
    UnloadTexture(tex);
    UnloadImage(Img);
    CloseWindow();

    return 0;
}