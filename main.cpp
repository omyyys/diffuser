#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "raylib/include/raylib.h"
#include "solver/solver.hpp"
    
int main() {
    const int scale = 10;
    const int nx = 50;
    const int fps = 240;

    FluidCube fluid = FluidCube(nx, 0, 0, 1.0/fps);

    InitWindow(nx*scale, nx*scale, "test");
    SetTargetFPS(fps);

    while (!WindowShouldClose()) {
        ClearBackground(BLACK);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            fluid.add_density(GetMouseX() / scale, GetMouseY() / scale, 1.0);
            fluid.add_velocity(GetMouseX() / scale, GetMouseY() / scale, 1.0, 1.0);
        }

        fluid.step();
        fluid.render_density(scale);

        DrawFPS(GetScreenWidth() - 100, 10);
    }

    CloseWindow();

    return 0;
}
