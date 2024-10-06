// TODO: Perlin noise for fade and for random spawn
// TODO: Help icons for info about the simulation

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"
 
#include "../solver/solver.hpp"

using namespace std;

constexpr int FPS = 120;

constexpr double DIFF_DEF = 0.0;
constexpr double VISC_DEF = 0.0;
constexpr double PARTICLE_SIZE_DEF = 5.0; 
constexpr int COLOUR_CHANGE_PERIOD = 4;
constexpr double BASE_FADE = 5e-4;

constexpr int STD_FONT_SIZE = 20;

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

typedef struct {
    int current_resolution_index; 
    int nx;
    int scale;
    float visc;
    float diff;
    Color colour;
    bool blend;
    bool colour_change;
    bool random_spawn;
    int random_splat_period;
    bool paused;
    FieldType rendering_field;

    Vector2 mouse;
    Vector2 mouse_prev;

    int timesteps;

} UI_State;

void set_default_state(UI_State& state) {
    state.current_resolution_index = 2; 
    state.nx = resolutions[state.current_resolution_index];
    state.scale = WIDTH / state.nx; 
    state.visc = VISC_DEF;
    state.diff = DIFF_DEF;
    state.colour = {.r = 65, .g = 70, .b = 186, .a = 255};
    state.blend = true;
    state.colour_change = false; 
    state.random_spawn = false; 
    state.random_splat_period = FPS;
    state.paused = false;
    state.rendering_field = DENSITY;

    state.mouse = Vector2Zero();
    state.mouse_prev = Vector2Zero(); 

    state.timesteps = 0;
}

double rand_double(double start, double end) {
    return translate_value(rand(), 0, RAND_MAX, start, end);
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

void normalise_window_size() {
    if (GetScreenWidth() != (WIDTH + GUI_WIDTH) || GetScreenHeight() != HEIGHT) {
        SetWindowSize(WIDTH + GUI_WIDTH, HEIGHT);
    }
}

bool mouse_is_in_window(Vector2 mouse_pos, int scale) {
    Rectangle mouse_rect = (Rectangle) { mouse_pos.x, mouse_pos.y, (float) scale, (float) scale };
    Rectangle screen_rect = (Rectangle) { 0, 0, (float) WIDTH, (float) HEIGHT };

    Rectangle collision = GetCollisionRec(mouse_rect, screen_rect);

    return (collision.width * collision.height) >= (mouse_rect.width * mouse_rect.height);    
}

void render(UI_State &state, Fluid &fluid) {
    if (state.rendering_field < 0) {
        state.rendering_field = (FieldType) (FIELDS_COUNT-1);
    }

    if (state.blend) BeginBlendMode(BLEND_ADDITIVE);

    if (state.rendering_field == DENSITY) { 
        fluid.render_density(state.scale, state.colour, PARTICLE_SIZE_DEF);

    } else if (state.rendering_field == VELOCITY) {
        fluid.render_velocity(state.scale, state.colour, PARTICLE_SIZE_DEF);

    } else if (state.rendering_field == VORTICITY) {
        fluid.render_vorticity(state.scale, state.colour, PARTICLE_SIZE_DEF); 
    }

    if (state.blend) EndBlendMode();
}

void get_ui_input(UI_State &state, Fluid &fluid) {
    DrawText("Viscosity", 920, 50, STD_FONT_SIZE, WHITE);
    GuiSlider((Rectangle) { 920, 80, 200, 25 }, "0", "1", &state.visc, 0.0, 0.01);
    fluid.visc = state.visc;

    DrawText("Diffusion", 920, 140, STD_FONT_SIZE, WHITE);
    GuiSlider((Rectangle) { 920, 170, 200, 25 }, "0", "0.01", &state.diff, 0.0, 0.01);
    fluid.diff = state.diff;

    DrawText("Colour", 920, 250, STD_FONT_SIZE, WHITE);
    GuiColorPicker((Rectangle) { 920, 280, 200, 200 }, NULL, &state.colour);

    DrawText("Additive Blending", 920, 550, STD_FONT_SIZE, WHITE);
    GuiCheckBox((Rectangle) { 1120, 550, 20, 20 }, NULL, &state.blend);

    DrawText("Auto Colours", 920, 600, STD_FONT_SIZE, WHITE);
    GuiCheckBox((Rectangle) { 1120, 600, 20, 20 }, NULL, &state.colour_change);

    DrawText("Auto Spawn", 920, 650, STD_FONT_SIZE, WHITE);
    GuiCheckBox((Rectangle) { 1120, 650, 20, 20 }, NULL, &state.random_spawn);

    if (state.random_spawn) {
        float random_spawn_period_float = state.random_splat_period;
        DrawText("Spawn Interval (Frames)", 920, 700, STD_FONT_SIZE, WHITE);
        GuiSlider((Rectangle) { 920, 730, 200, 25 }, "30", "240", &random_spawn_period_float, 30, 240);
        state.random_splat_period = random_spawn_period_float;
    }
}

void draw_info(UI_State state) {
    DrawText(TextFormat("grid resolution: %d", state.nx), 20, 20, STD_FONT_SIZE, WHITE);

    if (state.rendering_field == DENSITY) {
        DrawText("rendering density", 20, 50, STD_FONT_SIZE, WHITE);

    } else if (state.rendering_field == VELOCITY) {
        DrawText("rendering velocity", 20, 50, STD_FONT_SIZE, WHITE);

    } else if (state.rendering_field == VORTICITY) {
        DrawText("rendering vorticity", 20, 50, STD_FONT_SIZE, WHITE);

    } else {
        DrawText(TextFormat("error: field enum not exhausted: %d", state.rendering_field), 20, 50, STD_FONT_SIZE, WHITE);
    }

    if (state.paused) DrawText("paused", 20, 80, STD_FONT_SIZE, WHITE);
        
    DrawText("arrow keys to switch field", 20, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
    DrawText("c to clear screen", 20, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);
    DrawText("r to change grid resolution", 320, GetScreenHeight() - 80, STD_FONT_SIZE, WHITE);
    DrawText("esc to pause", 320, GetScreenHeight() - 50, STD_FONT_SIZE, WHITE);

    DrawFPS(GetScreenWidth() - 100, 10);
}

void handle_input(UI_State &state, Fluid &fluid) {
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && mouse_is_in_window(state.mouse, state.scale)) { 
        Vector2 mouse_travelled = Vector2Subtract(state.mouse, state.mouse_prev);

        fluid.add_density(state.mouse.x / state.scale, 
                            state.mouse.y / state.scale, 1.0);

        fluid.add_velocity(state.mouse.x / state.scale, 
                            state.mouse.y / state.scale, mouse_travelled.x, mouse_travelled.y); 
    }

    if (IsKeyPressed(KEY_C)) fluid.reset();

    if (IsKeyPressed(KEY_R)) {
        state.current_resolution_index = (state.current_resolution_index + 1) % (sizeof(resolutions) / sizeof(int));
        state.nx = resolutions[state.current_resolution_index];
        fluid = Fluid(state.nx, state.diff, state.visc, 1.0/FPS);
        state.scale = WIDTH / state.nx;
    }

    if (IsKeyPressed(KEY_RIGHT)) {
        state.rendering_field = (FieldType) ((state.rendering_field + 1) % FIELDS_COUNT); 

    } else if (IsKeyPressed(KEY_LEFT)) {
        state.rendering_field = (FieldType) ((state.rendering_field  - 1) % FIELDS_COUNT);
    }

    if (IsKeyPressed(KEY_ESCAPE)) {
        state.paused = !state.paused;
    }
}

int main() {
    srand(time(NULL));

    UI_State state;
    set_default_state(state);

    Fluid fluid(state.nx, state.diff, state.visc, 1.0/FPS);

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(WIDTH + GUI_WIDTH, HEIGHT, "test");
    GuiSetStyle(DEFAULT, TEXT_SIZE, STD_FONT_SIZE); 
    SetTargetFPS(FPS);
    SetExitKey(KEY_NULL);

    while (!WindowShouldClose()) {
        normalise_window_size();

        state.mouse_prev = state.mouse;
        state.mouse = GetMousePosition();

        BeginDrawing();
        ClearBackground(BLACK);

        if (!state.paused) {
            fluid.step();
            fluid.fade_density(BASE_FADE + state.visc * 0.001);
            state.timesteps++;
        }

        if (state.colour_change && state.timesteps % COLOUR_CHANGE_PERIOD == 0) {
            state.colour = get_colour(state.timesteps); 
        }

        if (state.random_spawn && state.timesteps % state.random_splat_period == 0) {
            random_spawn(fluid, state.nx);
        }

        render(state, fluid);
        get_ui_input(state, fluid);
        draw_info(state);
        handle_input(state, fluid);
 
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
