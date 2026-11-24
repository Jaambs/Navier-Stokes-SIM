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
#include <algorithm>
#include "core/grid.hpp"
#include "physics/BC.hpp"
#include "physics/cylinder.hpp"
#include "physics/mouse.hpp"
#include "physics/multigrid.hpp"
#include "physics/solver.hpp"
#include "physics/vorticity.hpp"
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
    const Real dt      = 1.0 / Real(60.0);
    Real sim_time      = 0.0f;
    SetTargetFPS(60);
    Real nu = Real(1.67e-3);
    Real force = Real(15.0);
    const Real U_in = 1.0f;
    const Real L_ref = 1.0f / 6.0f;
    Grid g{NW, NH, 1.0f / NH};
    CylinderMask cyl(g);
    MultigridSolver mg_poisson(g, cyl);
    bool inject_dye = true;
    bool use_advection = true;
    bool use_diffusion = true;
    bool use_projection = true;
    bool slip = false;

    // Physical size of the domain: NW*h x NH*h = 2 x 1 for a 2:1 channel
    const Real DOMAIN_W = NW * g.h;
    const Real DOMAIN_H = NH * g.h;
     
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

    // ---- Strohual Measurements ----

    int i_probe = static_cast<int>((cyl.xc_cyl + 3.0f * cyl.diameter) / g.h) + 1;
    int j_probe = g.ny / 2 + 1;
    //frecuency calculations and csv export variables
    Real prev_v = 0.0f;
    Real last_cross_time = 1.0f;
    Real strohual = 0.0f;

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
        std::fill(v.a.begin(), v.a.end(), 0.0f);
        std::fill(p.a.begin(), p.a.end(), 0.0f);
        std::fill(dyeR.a.begin(), dyeR.a.end(), 0.0f);
        std::fill(dyeG.a.begin(), dyeG.a.end(), 0.0f);
        std::fill(dyeB.a.begin(), dyeB.a.end(), 0.0f);
        prev_v = 0.0f;
        last_cross_time = -1.0f;
        strohual = 0.0f;
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
        if (IsKeyPressed(KEY_S)) slip = !slip;
        //change viscosity
        if (IsKeyPressed(KEY_O)) {
            nu = (nu <= 0.0f) ? 1e-5 : nu * 1.5f;
        }
        if (IsKeyPressed(KEY_P)) {
            nu = (nu <= Real(1e-7)) ? 0.0f : nu / 1.5f;
        }
        if (IsKeyPressed(KEY_I)) {
            nu = 0.0f;
        }

        //Change force strengh
        if (IsKeyPressed(KEY_K)) {
            force = (force <= Real(1e-2)) ? 5.0f : force * 1.5f;
        }
        if (IsKeyPressed(KEY_L)) {
            force = (force <= Real(1e-2)) ? force : force / 1.5f;
        }
        if (IsKeyPressed(KEY_J)) {
            force = 15;
        }
        if (IsKeyPressed(KEY_C)) {
            cyl.active = !cyl.active;
            mg_poisson.updateMask(cyl);
            if (cyl.active) {
                applyCylinderVelocityBC(u, v, cyl, slip);
                for (int i = 1; i <= g.nx; ++i) {
                    for (int j = 1; j <= g.ny; ++j) {
                        if (cyl.isSolid(i,j)) {
                            p(i,j) = Real(0.0);
                            dyeR(i, j) = dyeG(i, j) = dyeB(i ,j) = Real(0.0);
                        }
                    }
                }
            }
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
                dyeR(1, j) = 0.0f;
                dyeG(1, j) = 0.85f;
                dyeB(1, j) = 1.0f;
            } else if (bot_line) {
                // Línea inferior: Naranja
                dyeR(1, j) = 1.0f;
                dyeG(1, j) = 0.3f;
                dyeB(1, j) = 0.1f;
            } else {
                // Fluido limpio en el resto de la entrada (y cuando apagas con D)
                dyeR(1, j) = 0.0f;
                dyeG(1, j) = 0.0f;
                dyeB(1, j) = 0.0f;
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

            Real mx  =  pos.x / WIN_WIDTH  * DOMAIN_W;
            Real my  =  DOMAIN_H * (1.0f - pos.y / WIN_HEIGHT);
            Real mdx =  delta.x / WIN_WIDTH  * DOMAIN_W * 60.0f;
            Real mdy = -delta.y / WIN_HEIGHT * DOMAIN_H * 60.0f;

            //Make the applied colour change with time almost randomly using trigonometric smooth functions
            Real r_col = 0.5f + 0.5f * std::sin(Real(sim_time * 1.5f));
            Real b_col = 0.5f + 0.5f * std::sin(Real(sim_time * 1.5f + 2.094f));
            Real g_col = 0.5f + 0.5f * std::sin(Real(sim_time * 1.5f + 4.188f));

            applyMouseBrush(u, v, dyeR, dyeB, dyeG, mx, my, mdx, mdy, r_col, b_col, g_col, dt, g.h, force);
            applyChannelVelocityBC(u, v, U_in, true, slip);
        }

        //2. --- Advect ---
        if (use_advection) {
            advectU(un, u, v, dt, g.h);
            advectV(vn, u, v, dt, g.h);
            std::swap(u.a, un.a);
            std::swap(v.a, vn.a);
            applyChannelVelocityBC(u, v, U_in, true, slip);
            applyCylinderVelocityBC(u, v, cyl, slip);
        }

        //3. --- Diffuse ---
        if (use_diffusion) {
            std::swap(u.a, un.a);
            std::swap(v.a, vn.a);
            diffuseChannelU(u, un, nu, dt, g.h, U_in, 15);
            diffuseChannelV(v, vn, nu, dt, g.h, 15);
            applyCylinderVelocityBC(u, v, cyl, slip);
        }
        // 4. --- Project ---
        if (use_projection) {
            projectChannelMG(u, v, p, div, cyl, mg_poisson, dt, g.h, U_in, slip);
            applyCylinderVelocityBC(u, v, cyl, slip);
        }

        //strohual meditions
        Real v_probe = v(i_probe, j_probe);

        if (prev_v < 0.0f && v_probe >= 0.0f) {
            if (last_cross_time > 0.0f) {
                Real period = sim_time - last_cross_time;
                Real freq = 1.0f / period;
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
                if (cyl.active) {
                    if (cyl.isSolid(i,j)) {
                        dyeR(i, j) = 0.0f;
                        dyeG(i, j) = 0.0f;
                        dyeB(i, j) = 0.0f;
                    }
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
                    unsigned char r_c = static_cast<unsigned char>(std::clamp(dyeR(i, j), 0.0f, 1.0f) * 255.0f);
                    unsigned char g_c = static_cast<unsigned char>(std::clamp(dyeG(i, j), 0.0f, 1.0f) * 255.0f);
                    unsigned char b_c = static_cast<unsigned char>(std::clamp(dyeB(i, j), 0.0f, 1.0f) * 255.0f);
                    pixels[py * NW + px] = Color{r_c, g_c, b_c, 255};
                    if (cyl.active) {
                        if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};}
                    }
                }
                else if (view == ViewMode::VORTICITY){
                    pixels[py * NW + px] = vorticityToColor(omega(i,j), 15.0f);
                    if (cyl.active) {
                        if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};}
                    }
                }
                else {
                    Real uc = uAtCenter(u, i, j);
                    Real vc = vAtCenter(v, i, j);
                    Real speed = std::hypot(uc, vc);
                    pixels[py * NW + px] = speedToColor(speed, 4);
                    if (cyl.active) {
                        if (cyl.isSolid(i, j)) {pixels[py * NW + px] = Color{180, 180, 180, 255};}
                    }
                }
            }
        }

        UpdateTexture(tex, Img.data);

        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(tex, Vector2{0.0f, 0.0f}, 0.0f, static_cast<float>(SCALE), WHITE);
            DrawFPS(10, 10);
            if (nu > 0.0f) {
                Real Re = (U_in * L_ref) / nu;
                DrawText(TextFormat("nu =%.2e | Re = %.1f | St = %.3f", nu, Re, strohual), 10, 35, 15, YELLOW);
            }
            if (nu == 0.0f) {
                DrawText("nu = 0.0f | Re = inf", 10, 35, 15, YELLOW);
            }
            DrawText(TextFormat("force =%.3e", force), 10, 60, 15, YELLOW);
            DrawText(TextFormat("[1] Adv:%s  [2] Diff:%s  [3] Pres:%s  [S] Slip:%s",
                                use_advection  ? "ON" : "OFF",
                                use_diffusion  ? "ON" : "OFF",
                                use_projection ? "ON" : "OFF",
                                slip ? "ON" : "OFF"),
                                10, 85, 15, ORANGE);
            DrawText("Controls: [C] cylinder control | [R] Reset | [I,O,P] Viscosity | [J,K,L] Brush | [D] Dye Injection | [TAB] Change View", 10, 110, 15, LIGHTGRAY);
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