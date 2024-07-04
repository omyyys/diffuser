#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "raylib/include/raylib.h"
#include "solver/solver.hpp"

#define FPS 240

bool mouse_is_in_window(Vector2 mouse_pos, int scale) {
    if (mouse_pos.y < 20 || mouse_pos.y > GetScreenHeight() - 20) return false;

    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) GetScreenWidth(), (float) GetScreenHeight() };
    
    return CheckCollisionRecs(mouse_rect, screen_rect);    
}

    
int main() {
    int scale = 10;
    int nx = 80;
    
    FluidCube fluid = FluidCube(nx, 0, 2, 1.0/FPS);

    InitWindow(nx*scale, nx*scale, "test");
    SetTargetFPS(FPS);

    while (!WindowShouldClose()) {
        ClearBackground(BLACK);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_pos = GetMousePosition();
            
            if (mouse_is_in_window(mouse_pos, scale)) {
                fluid.add_density(mouse_pos.x / scale, mouse_pos.y / scale, 1.0);
                fluid.add_velocity(mouse_pos.x / scale, mouse_pos.y / scale, 1.0, 1.0);

            }
        }

        fluid.step();
        fluid.render_density(scale);

        DrawFPS(GetScreenWidth() - 100, 10);
    }

    CloseWindow();

    return 0;
}
