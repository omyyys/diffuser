#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <unistd.h>

#include "../raylib/include/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"

#include "../solver/solver.hpp"

#define FPS 60

#define DIFF_DEF 0.0004
#define VISC_DEF 20.0
#define PARTICLE_SIZE_DEF 6

#define STD_FONT_SIZE 20

#define TEST false

using namespace std;

void test_with_print() {
    Fluid fluid1(5, DIFF_DEF, 0.0, 1.0/FPS);
    Fluid fluid2(5, DIFF_DEF, VISC_DEF, 1.0/FPS);

    fluid1.add_density(2, 1, 1.0);
    fluid1.add_density(2, 2, 1.0);
    fluid1.add_density(2, 3, 1.0);
    fluid1.add_density(2, 4, 1.0);

    for (int i = 0; i < 60; i++) {
        fluid1.step();
    }

    fluid1.print_density();
    printf("fluid 1: diff=%lf, visc=%lf\n", fluid1.diff, fluid1.visc);

    printf("---------------------------\n");

    fluid2.add_density(2, 1, 1.0);
    fluid2.add_density(2, 2, 1.0);
    fluid2.add_density(2, 3, 1.0);
    fluid2.add_density(2, 4, 1.0);

    for (int i = 0; i < 60; i++) {
        fluid2.step();
    }

    fluid2.print_density();
    printf("fluid 2: diff=%lf, visc=%lf\n", fluid2.diff, fluid2.visc);
}

bool mouse_is_in_window(Vector2 mouse_pos, int nx, int scale) {
    if (mouse_pos.y <= 20 || mouse_pos.y >= nx*scale - 20) return false;

    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) nx*scale, (float) nx*scale };
    
    return CheckCollisionRecs(mouse_rect, screen_rect);    
}

int main() {
    // FILE *f = fopen("test.txt", "w");

    int scale = 10;
    int nx = 80;
    int gui_edge = nx*scale;

    float particle_size = PARTICLE_SIZE_DEF;
    float visc = VISC_DEF;
    float diff = DIFF_DEF;
    Color colour = ORANGE;
    bool blend = true;

    int timesteps = 0;

    Fluid fluid(nx, diff, visc, 1.0/FPS);

    if (TEST) {
        test_with_print();
        return 0;
    }

    InitWindow(gui_edge + 400, gui_edge, "test");
    SetTargetFPS(FPS);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        if (IsKeyDown(KEY_C)) fluid.reset();

        fluid.step();

        if (blend) BeginBlendMode(BLEND_ADDITIVE);
        fluid.render_density(scale, colour, particle_size);
        if (blend) EndBlendMode();

        DrawFPS(GetScreenWidth() - 100, 10);

        DrawText("Viscosity", 920, 70, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 100, 200, 25 }, "0", "200", &visc, 0.0, 200.0);
        fluid.visc = visc;

        DrawText("Diffusion", 920, 160, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 190, 200, 25 }, "0", "0.001", &diff, 0.0, 0.001);
        fluid.diff = diff;

        DrawText("Scale", 920, 250, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 280, 200, 25 }, "0", "10.0", &particle_size, 0.0, 10.0);

        DrawText("Colour", 920, 340, STD_FONT_SIZE, WHITE);
        GuiColorPicker((Rectangle) { 920, 370, 200, 200 }, NULL, &colour);

        DrawText("Additive Blending", 920, 630, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 630, 20, 20 }, NULL, &blend);

        DrawText(TextFormat("timesteps: %d", timesteps), 20, 20, STD_FONT_SIZE, WHITE);

        DrawText("drag mouse to draw fluid", 920, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
        DrawText("c to clear screen", 920, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_pos = GetMousePosition();
            
            if (mouse_is_in_window(mouse_pos, nx, scale)) {
                fluid.add_density(round(mouse_pos.x / scale), 
                                  round(mouse_pos.y / scale), 1.0);

                fluid.add_velocity(round(mouse_pos.x / scale), 
                                   round(mouse_pos.y / scale), 1.0, 1.0);
            }
        } 
 
        EndDrawing();

        timesteps++;
    }

    CloseWindow();

    return 0;
}
