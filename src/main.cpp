#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <unistd.h>

#include "../raylib/include/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"

#include "../solver/solver.hpp"

#define FPS 120

#define DIFF_DEF 0.05
#define VISC_DEF 10.0
#define PARTICLE_SIZE_DEF 5
#define RANDOM_SPLAT_COUNT 200

#define STD_FONT_SIZE 20

#define TEST false
#define BENCHMARK false
#define BENCHMARK_TIMESTEPS 1000

using namespace std;

double rand_double(double start, double end) { 
    return lin_interpolation(rand(), 0, RAND_MAX, start, end);
}

int rand_int(int start, int end) {
    return (rand() % (end - start + 1)) + start; 
} 

void random_spawn(Fluid &fluid, int splats) {
    for (int i = 0; i < splats; i++) {
        int index_x = rand_int(3, fluid.nx - 3); 
        int index_y = rand_int(3, fluid.nx - 3); 

        fluid.add_density(index_x, index_y, rand_double(0.5, 1.0));
        fluid.add_velocity(index_x, index_y, rand_double(0.5, 1.0), rand_double(0.5, 1.0));
    }
}

void test_with_print() {
    Fluid fluid1(7, DIFF_DEF, 0.0, 1.0/FPS);
    Fluid fluid2(7, DIFF_DEF, 10.0, 1.0/FPS);

    fluid1.add_density(2, 1, 1.0);
    fluid1.add_density(2, 2, 1.0);
    fluid1.add_density(2, 3, 1.0);
    fluid1.add_density(2, 4, 1.0);

    for (int i = 0; i < 60; i++) { 
        fluid1.step();
    }

    fluid1.print_density();
    printf("fluid 1: diff=%lf, visc=%lf, general test\n", fluid1.diff, fluid1.visc);

    fluid2.add_density(2, 1, 1.0);
    fluid2.add_density(2, 2, 1.0);
    fluid2.add_density(2, 3, 1.0);
    fluid2.add_density(2, 4, 1.0);

    for (int i = 0; i < 60; i++) { 
        fluid2.step();
    }

    fluid2.print_density();
    printf("fluid 2: diff=%lf, visc=%lf, general test\n", fluid2.diff, fluid2.visc);
}

bool mouse_is_in_window(Vector2 mouse_pos, int nx, int scale) {
    if (mouse_pos.y <= 20 || mouse_pos.y >= nx*scale - 20) return false;

    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) nx*scale, (float) nx*scale };
    
    return CheckCollisionRecs(mouse_rect, screen_rect);    
}

int main() {
    srand(time(NULL));
    
    int scale = 20;
    int nx = 40;
    int gui_edge = nx*scale;

    float particle_size = PARTICLE_SIZE_DEF;
    float visc = VISC_DEF;
    float diff = DIFF_DEF;
    Color colour = ORANGE;
    bool blend = true; 

    int timesteps = 0;
    int fps_sum = 0;

    Fluid fluid(nx, diff, visc, 1.0/FPS);

    // start_splat(fluid);

    if (TEST) {
        test_with_print();
        return 0;
    }

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(gui_edge + 400, gui_edge, "test");
    GuiSetStyle(DEFAULT, TEXT_SIZE, STD_FONT_SIZE); 
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
        GuiSlider((Rectangle) { 920, 100, 200, 25 }, "0", "50", &visc, 0.0, 50.0);
        fluid.visc = visc;

        DrawText("Diffusion", 920, 160, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 190, 200, 25 }, "0", "1.0", &diff, 0.0, 0.1);
        fluid.diff = diff;

        DrawText("Scale", 920, 250, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 280, 200, 25 }, "0", "10.0", &particle_size, 0.0, 10.0);

        DrawText("Colour", 920, 340, STD_FONT_SIZE, WHITE);
        GuiColorPicker((Rectangle) { 920, 370, 200, 200 }, NULL, &colour);

        DrawText("Additive Blending", 920, 630, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 630, 20, 20 }, NULL, &blend);

        if (GuiButton((Rectangle) { 920, 680, 150, 40 }, "Random splat")) {
            for (int i = 0; i < 5; i++)
                random_spawn(fluid, RANDOM_SPLAT_COUNT);
        }
        
        DrawText(TextFormat("timesteps: %d", timesteps), 20, 20, STD_FONT_SIZE, WHITE);

        DrawText("drag mouse to draw fluid", 20, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
        DrawText("c to clear screen", 20, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);

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
        if (BENCHMARK) {
            fps_sum += GetFPS();
            if (timesteps >= BENCHMARK_TIMESTEPS) break;
        }
    }

    CloseWindow();

    if (BENCHMARK) printf("avg fps over %d timesteps: %0.1lf\n", timesteps, double(fps_sum) / timesteps);

    return 0;
}
