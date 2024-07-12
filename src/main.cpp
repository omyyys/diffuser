#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"
 
#include "../solver/solver.hpp"

#define FPS 120

#define DIFF_DEF 0.005
#define VISC_DEF 0.0

#define STD_FONT_SIZE 20

#define TEST false
#define BENCHMARK false
#define BENCHMARK_TIMESTEPS 1000

#define WIDTH 800
#define HEIGHT WIDTH
#define GUI_WIDTH 400

typedef enum {
    DENSITY = 0,
    VELOCITY,
    VORTICITY,
    FIELDS_COUNT,
    
} FieldType;

int resolutions[] = {
    80, 100, 160 
};

using namespace std;

double rand_double(double start, double end) { 
    return lin_interpolation(rand(), 0, RAND_MAX, start, end);
}

void random_spawn(Fluid &fluid, int splats) {
    for (int i = 0; i < splats; i++) {
        int index_x = GetRandomValue(3, fluid.nx - 3); 
        int index_y = GetRandomValue(3, fluid.nx - 3); 

        fluid.add_density(index_x, index_y, rand_double(0.5, 1.0));
        fluid.add_velocity(index_x, index_y, rand_double(0.5, 1.0), rand_double(0.5, 1.0));
    }
}

Color get_colour(int timesteps) {
    return ColorFromHSV (
        timesteps % 360,
        100, 
        60
    );
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

bool mouse_is_in_window(Vector2 mouse_pos, int scale) {
    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) WIDTH, (float) HEIGHT };

    Rectangle collision = GetCollisionRec(mouse_rect, screen_rect);

    return (collision.width * collision.height) >= (mouse_rect.width * mouse_rect.height);    
}

int main() {
    srand(time(NULL));

    int current_resolution_index = 1; 
    int nx = resolutions[current_resolution_index];
    int scale = WIDTH / nx; 

    float particle_size = nx / 20.0;
    float visc = VISC_DEF;
    float diff = DIFF_DEF;
    Color colour = {.r = 0, .g = 0, .b = 92, .a = 255};
    bool blend = true;
    bool colour_change = false; 
    FieldType rendering_field = DENSITY;

    int timesteps = 0;
    int fps_sum = 0;
    
    Fluid fluid(nx, diff, visc, 1.0/FPS);
    
    if (TEST) {
        test_with_print();
        return 0;
    }

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(WIDTH + GUI_WIDTH, HEIGHT, "test");
    GuiSetStyle(DEFAULT, TEXT_SIZE, STD_FONT_SIZE); 
    SetTargetFPS(FPS);

    Vector2 pmouse = GetMousePosition(); 

    while (!WindowShouldClose()) {
        if (GetScreenWidth() != (WIDTH + GUI_WIDTH) || GetScreenHeight() != HEIGHT) 
            SetWindowSize(WIDTH + GUI_WIDTH, HEIGHT);

        BeginDrawing();
        ClearBackground(BLACK);

        if (IsKeyPressed(KEY_C)) fluid.reset();

        if (IsKeyPressed(KEY_R)) {
            current_resolution_index = (current_resolution_index + 1) % (sizeof(resolutions) / sizeof(int));
            nx = resolutions[current_resolution_index];
            fluid = Fluid(nx, diff, visc, 1.0/FPS);
            scale = WIDTH / nx;
            particle_size = nx / 20.0;
        }

        if (IsKeyPressed(KEY_RIGHT)) {
            rendering_field = (FieldType) ((rendering_field + 1) % FIELDS_COUNT); 

        } else if (IsKeyPressed(KEY_LEFT)) {
            rendering_field = (FieldType) ((rendering_field - 1) % FIELDS_COUNT);
        }

        if (rendering_field < 0) rendering_field = (FieldType) (FIELDS_COUNT-1);

        fluid.step();

        if (blend) BeginBlendMode(BLEND_ADDITIVE);

        if (rendering_field == DENSITY) fluid.render_density(scale, colour, particle_size);
        else if (rendering_field == VELOCITY) fluid.render_velocity(scale, colour, particle_size);
        else if (rendering_field == VORTICITY) fluid.render_vorticity(scale, colour, particle_size);     

        if (blend) EndBlendMode();

        // Change colour every 5 frames
        if (colour_change && timesteps % 10 == 0) colour = get_colour(timesteps);

        fluid.fade_density(5e-4 + visc * 0.001);

        DrawFPS(GetScreenWidth() - 100, 10);

        DrawText("Viscosity", 920, 50, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 80, 200, 25 }, "0", "10", &visc, 0.0, 10.0);
        fluid.visc = visc;

        DrawText("Diffusion", 920, 140, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 170, 200, 25 }, "0", "0.1", &diff, 0.0, 0.1);
        fluid.diff = diff;

        /*
        DrawText("Size", 920, 230, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 260, 200, 25 }, "0", "10", &particle_size, 0.0, 10.0);
        */

        DrawText("Colour", 920, 320, STD_FONT_SIZE, WHITE);
        GuiColorPicker((Rectangle) { 920, 350, 200, 200 }, NULL, &colour);

        DrawText("Additive Blending", 920, 610, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 610, 20, 20 }, NULL, &blend);

        DrawText("Auto Colours", 920, 660, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 660, 20, 20 }, NULL, &colour_change);

        if (GuiButton((Rectangle) { 920, 710, 200, 40 }, "Random splat")) {
            for (int i = 0; i < 5; i++)
                random_spawn(fluid, nx);
        }
   
        DrawText(TextFormat("timesteps: %d", timesteps), 20, 20, STD_FONT_SIZE, WHITE);
        DrawText(TextFormat("grid resolution: %d", nx), 20, 50, STD_FONT_SIZE, WHITE);

        if (rendering_field == DENSITY)
            DrawText("rendering density", 20, 80, STD_FONT_SIZE, WHITE);

        else if (rendering_field == VELOCITY)
            DrawText("rendering velocity", 20, 80, STD_FONT_SIZE, WHITE);

        else if (rendering_field == VORTICITY)
            DrawText("rendering vorticity", 20, 80, STD_FONT_SIZE, WHITE);

        else 
            DrawText(TextFormat("error: field enum not exhausted: %d", rendering_field), 20, 80, STD_FONT_SIZE, WHITE);
        
        DrawText("drag mouse to draw fluid", 20, GetScreenHeight() - 140, STD_FONT_SIZE, WHITE);
        DrawText("arrow keys to switch rendered field", 20, GetScreenHeight() - 110, STD_FONT_SIZE, WHITE);
        DrawText("c to clear screen", 20, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
        DrawText("r to change grid resolution", 20, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);
        
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mouse_pos = GetMousePosition();
            Vector2 mouse_travelled = Vector2Subtract(mouse_pos, pmouse);
            (void) mouse_travelled;
            
            if (mouse_is_in_window(mouse_pos, scale)) {
                fluid.add_density(round(mouse_pos.x / scale), 
                                  round(mouse_pos.y / scale), 1.0);

                fluid.add_velocity(round(mouse_pos.x / scale), 
                                   round(mouse_pos.y / scale), 3.0, 3.0); 
            }

        }
 
        EndDrawing();

        timesteps++;
        if (BENCHMARK) {
            fps_sum += GetFPS();
            if (timesteps >= BENCHMARK_TIMESTEPS) break;
        }

        pmouse = GetMousePosition();
    }

    CloseWindow();

    if (BENCHMARK) printf("avg fps over %d timesteps: %0.1lf\n", timesteps, double(fps_sum) / timesteps);

    return 0;
}
