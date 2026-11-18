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

int main(){

    // ---- SIM PARAMETERS ----

    const int N        = 128;
    const int SCALE    = 6;
    const int WIN_SIZE = N * SCALE;
    const double dt    = 1.0 / 60.0;
    double sim_time    = 0.0;
    Grid g{N, N, 1.0 / N};
     
    // ---- Velocity Divergence and Pressure Fields ----

    FaceX u(g), un(g);
    FaceY v(g), vn(g);
    CellField div(g), p(g);
    
    // ---- RGB Dye Fields ----

    CellField dyeR(g), dyeR_new(g);
    CellField dyeG(g), dyeG_new(g);
    CellField dyeB(g), dyeB_new(g);

    // ---- Inicialization of the program ----

    InitWindow(WIN_SIZE, WIN_SIZE, "Navier-Stokes Von Karman Simulation (MAC 128x128)");
    SetTargetFPS(60);

    Image Img = GenImageColor(N, N, BLACK);
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

        // For the simulation loop we follow the stablished Chorin schematic
        // Add forces -> advect -> diffuse -> project -> advect dye 
        //Then we render the result

        //1. ---Mouse Forces---
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 pos = GetMousePosition();
            Vector2 delta = GetMouseDelta();

            //Transform screen pixels to physic coordinates [0, 1] doing while so a Y axis inversion

            double mx = pos.x / WIN_SIZE;
            double my = 1.0 - (pos.y / WIN_SIZE);
            double mdx = delta.x / WIN_SIZE * 60.0;
            double mdy = - delta.y / WIN_SIZE * 60.0;

            //Make the applied colour change with time almost randomly using trigonometric smooth functions
            double r_col = 0.5 + 0.5 * std::sin(sim_time * 1.5);
            double b_col = 0.5 + 0.5 * std::sin(sim_time * 1.5 + 2.094);
            double g_col = 0.5 + 0.5 * std::sin(sim_time * 1.5 + 4.188);

            applyMouseBrush(u, v, dyeR, dyeB, dyeG, mx, my, mdx, mdy, r_col, b_col, g_col, dt, g.h);
            applyVelocityBC(u, v);
        }

        //2. --- Advect ---
        advectU(un, u, v, dt, g.h);
        advectV(vn, u, v, dt, g.h);
        std::swap(u.a, un.a);
        std::swap(v.a, vn.a);
        applyVelocityBC(u, v);

        // 3. --- Project ---
        project(u, v, p, div, dt, g.h, 60);

        // 4. --- Advect Dye ---

        advectScalar(dyeR_new, dyeR, u, v, dt, g.h);
        advectScalar(dyeG_new, dyeG, u, v, dt, g.h);
        advectScalar(dyeB_new, dyeB, u, v, dt, g.h);
        std::swap(dyeR_new, dyeR);
        std::swap(dyeG_new, dyeG);
        std::swap(dyeB_new, dyeB);

        // ----------- RENDERING -------------
        Color* pixels = static_cast<Color*>(Img.data);
        for (int py = 0; py < N; ++py) {
            int j = N - py;
            for (int px = 0; px < N; ++px) {
                int i = px + 1;

                unsigned char r_c = static_cast<unsigned char>(std::clamp(dyeR(i, j), 0.0, 1.0) * 255.0);
                unsigned char g_c = static_cast<unsigned char>(std::clamp(dyeG(i, j), 0.0, 1.0) * 255.0);
                unsigned char b_c = static_cast<unsigned char>(std::clamp(dyeB(i, j), 0.0, 1.0) * 255.0);

                pixels[py * N + px] = Color{r_c, g_c, b_c, 255};
            }
        }
        UpdateTexture(tex, Img.data);

        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(tex, Vector2{0.0f, 0.0f}, 0.0f, static_cast<float>(SCALE), WHITE);
            DrawFPS(10, 10);
            DrawText("Press [R] to reset", 10, 35, 18, LIGHTGRAY);
        EndDrawing();
    }
    UnloadTexture(tex);
    UnloadImage(Img);
    CloseWindow();

    return 0;
}