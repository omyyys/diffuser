#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "raylib/include/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui/src/raygui.h"

#include "solver/solver.hpp"

#define FPS 240

#define STD_FONT_SIZE 20

bool mouse_is_in_window(Vector2 mouse_pos, int nx, int scale) {
    if (mouse_pos.y < 20 || mouse_pos.y > nx*scale - 20) return false;

    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) nx*scale, (float) nx*scale };
    
    return CheckCollisionRecs(mouse_rect, screen_rect);    
}
    
int main() {
    int scale = 10;
    int nx = 80;
    int gui_edge = nx*scale;

    float particle_size = 5;
    float visc = 10.0;
    float diff = 0.0;
    Color colour = DARKBLUE;
    
    FluidCube fluid = FluidCube(nx, diff, visc, 1.0/FPS);

    InitWindow(gui_edge + 400, gui_edge, "test");
    SetTargetFPS(FPS);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);

        DrawText("Viscosity", 920, 70, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 100, 200, 25 }, "0", "100", &visc, 0, 100);
        fluid.visc = visc;

        DrawText("Diffusion", 920, 160, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 190, 200, 25 }, "0", "0.0005", &diff, 0, 0.0005);
        fluid.diff = diff;

        DrawText("Scale", 920, 250, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 280, 200, 25 }, "0", "10", &particle_size, 0, 10);

        DrawText("Colour", 920, 340, STD_FONT_SIZE, WHITE);
        GuiColorPicker((Rectangle) { 920, 370, 200, 200 }, NULL, &colour);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_pos = GetMousePosition();
            
            if (mouse_is_in_window(mouse_pos, nx, scale)) {
                fluid.add_density(mouse_pos.x / scale, mouse_pos.y / scale, 1.0);
                fluid.add_velocity(mouse_pos.x / scale, mouse_pos.y / scale, 1.0, 1.0);
            }
        }

        fluid.step();
        fluid.render_density(scale, colour, particle_size);

        DrawFPS(GetScreenWidth() - 100, 10);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
