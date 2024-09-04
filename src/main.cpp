#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"
 
#include "../solver/solver.hpp"

constexpr int FPS = 120;

constexpr double DIFF_DEF = 0.0;
constexpr double VISC_DEF = 0.0;
constexpr double PARTICLE_SIZE_DEF = 5.0; 

constexpr int STD_FONT_SIZE = 20;

constexpr bool BENCHMARK = false;
constexpr int BENCHMARK_TIMESTEPS = 1000;

constexpr int WIDTH = 800;
constexpr int HEIGHT = 800;
constexpr int GUI_WIDTH = 400;

constexpr int resolutions[] = { 80, 100, 160, 200 };

typedef enum {
    DENSITY = 0,
    VELOCITY,
    VORTICITY,
    FIELDS_COUNT,
    
} FieldType;

using namespace std;

double rand_double(double start, double end) { 
    return lin_interpolation(rand(), 0, RAND_MAX, start, end)    ;
}

void random_spawn(Fluid &fluid, int splats) {
    for (int i = 0; i < splats; i++) {
        int index_x = GetRandomValue(3, fluid.nx - 3); 
        int index_y = GetRandomValue(3, fluid.nx - 3); 

        fluid.add_density(index_x, index_y, rand_double(-100, 100));
        fluid.add_velocity(index_x, index_y, rand_double(-100, 100), rand_double(-100, 100));
    }
}

Color get_colour(int timesteps) {
    return ColorFromHSV (
        timesteps % 360,
        1, 
        0.6
    );
}

bool mouse_is_in_window(Vector2 mouse_pos, int scale) {
    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) WIDTH, (float) HEIGHT };

    Rectangle collision = GetCollisionRec(mouse_rect, screen_rect);

    return (collision.width * collision.height) >= (mouse_rect.width * mouse_rect.height);    
}

int main() {
    srand(time(NULL));

    int current_resolution_index = 2; 
    int nx = resolutions[current_resolution_index];
    int scale = WIDTH / nx; 
   
    float visc = VISC_DEF;
    float diff = DIFF_DEF;
    Color colour = {.r = 235, .g = 128, .b = 52, .a = 255};
    bool blend = true;
    bool colour_change = false; 
    bool paused = false;
    FieldType rendering_field = DENSITY;

    int timesteps = 0;
    int fps_sum = 0;
    
    Fluid fluid(nx, diff, visc, 1.0/FPS);

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(WIDTH + GUI_WIDTH, HEIGHT, "test");
    GuiSetStyle(DEFAULT, TEXT_SIZE, STD_FONT_SIZE); 
    SetTargetFPS(FPS);
    SetExitKey(KEY_NULL);

    Vector2 mouse = Vector2Zero();
    Vector2 prev_mouse = Vector2Zero(); 

    while (!WindowShouldClose()) {
        if (GetScreenWidth() != (WIDTH + GUI_WIDTH) || GetScreenHeight() != HEIGHT) 
            SetWindowSize(WIDTH + GUI_WIDTH, HEIGHT);

        prev_mouse = mouse;
        mouse = GetMousePosition();

        BeginDrawing();
        ClearBackground(BLACK);

        if (IsKeyPressed(KEY_C)) fluid.reset();

        if (IsKeyPressed(KEY_R)) {
            current_resolution_index = (current_resolution_index + 1) % (sizeof(resolutions) / sizeof(int));
            nx = resolutions[current_resolution_index];
            fluid = Fluid(nx, diff, visc, 1.0/FPS);
            scale = WIDTH / nx;
        }

        if (IsKeyPressed(KEY_RIGHT)) {
            rendering_field = (FieldType) ((rendering_field + 1) % FIELDS_COUNT); 

        } else if (IsKeyPressed(KEY_LEFT)) {
            rendering_field = (FieldType) ((rendering_field - 1) % FIELDS_COUNT);
        }

        if (IsKeyPressed(KEY_ESCAPE)) paused = !paused;

        if (rendering_field < 0) rendering_field = (FieldType) (FIELDS_COUNT-1);

        if (!paused) {
            fluid.step();
            fluid.fade_density(5e-4 + visc * 0.001);
            timesteps++;
        }

        if (blend) BeginBlendMode(BLEND_ADDITIVE);

        if (rendering_field == DENSITY) { 
            fluid.render_density(scale, colour, PARTICLE_SIZE_DEF);

        } else if (rendering_field == VELOCITY) {
            fluid.render_velocity(scale, colour, PARTICLE_SIZE_DEF);

        } else if (rendering_field == VORTICITY) {
            fluid.render_vorticity(scale, colour, PARTICLE_SIZE_DEF); 
        }

        if (blend) EndBlendMode();

        // Change colour every 10 frames
        if (colour_change && timesteps % 4 == 0) colour = get_colour(timesteps); 

        DrawFPS(GetScreenWidth() - 100, 10);

        DrawText("Viscosity", 920, 50, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 80, 200, 25 }, "0", "5", &visc, 0.0, 5.0);
        fluid.visc = visc;

        DrawText("Diffusion", 920, 140, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 170, 200, 25 }, "0", "0.1", &diff, 0.0, 0.1);
        fluid.diff = diff;

        DrawText("Colour", 920, 320, STD_FONT_SIZE, WHITE);
        GuiColorPicker((Rectangle) { 920, 350, 200, 200 }, NULL, &colour);

        DrawText("Additive Blending", 920, 610, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 610, 20, 20 }, NULL, &blend);

        DrawText("Auto Colours", 920, 660, STD_FONT_SIZE, WHITE);
        GuiCheckBox((Rectangle) { 1120, 660, 20, 20 }, NULL, &colour_change);

        if (GuiButton((Rectangle) { 920, 710, 200, 40 }, "Random splat")) {
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

        if (paused) DrawText("paused", 20, 110, STD_FONT_SIZE, WHITE);
        ;
        DrawText("arrow keys to switch field", 20, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
        DrawText("c to clear screen", 20, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);
        DrawText("r to change grid resolution", 320, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
        DrawText("esc to pause", 320, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);
        
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && mouse_is_in_window(mouse, scale)) { 
            Vector2 mouse_travelled = Vector2Subtract(mouse, prev_mouse);

            fluid.add_density(round(mouse.x / scale), 
                              round(mouse.y / scale), 1.0);

            fluid.add_velocity(round(mouse.x / scale), 
                               round(mouse.y / scale), mouse_travelled.x, mouse_travelled.y); 
        }
 
        EndDrawing();

        if (BENCHMARK) {
            fps_sum += GetFPS();
            if (timesteps >= BENCHMARK_TIMESTEPS) break;
        }
    }

    CloseWindow();

    if (BENCHMARK) printf("avg fps over %d timesteps: %0.1lf\n", timesteps, double(fps_sum) / timesteps);

    return 0;
}
