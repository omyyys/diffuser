// Uses Raylib 5.5

// TODO: Help icons for info about the simulation

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "../raygui/src/raygui.h"
 
#include "../include/solver.hpp"

using namespace std;

constexpr int FPS = 240;

constexpr double DIFF_DEF = 0.0;
constexpr double VISC_DEF = 0.0;
constexpr double PARTICLE_SIZE_DEF = 5.0; 
constexpr double VELOCITY_ARROW_LENGTH_DEF = 20.0; 

constexpr int COLOUR_CHANGE_PERIOD = 4;
constexpr double BASE_FADE = 3e-4;

constexpr int STD_FONT_SIZE = 20;

constexpr int WIDTH = 800;
constexpr int HEIGHT = 800;
constexpr int GUI_WIDTH = 400;

constexpr float HUE_SHIFT_FACTOR = 0.5;
constexpr int MAX_HUE = 360;
constexpr float SATURATION = 1;
constexpr float BRIGHTNESS = 0.6;

constexpr Color COLOUR_DEF = {.r = 65, .g = 70, .b = 186, .a = 255};

constexpr int resolutions[] = { 80, 100, 160, 200 };

int resolution_next_index(int current) {
    return (current + 1) % (sizeof(resolutions) / sizeof(int));
}

typedef enum {
    DENSITY = 0,
    VELOCITY,
    VORTICITY,
    FIELDS_COUNT,
    
} FieldType;

FieldType field_type_next(FieldType current) {
    return (FieldType) ((current + 1) % FIELDS_COUNT);
}

FieldType field_type_prev(FieldType current) {
    return (FieldType) ((current - 1) % FIELDS_COUNT);
}

typedef struct {
    int resolution_index; 
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

void ui_state_init(UI_State& state) {
    state.resolution_index = 2; 
    state.nx = resolutions[state.resolution_index];
    state.scale = WIDTH / state.nx; 
    state.visc = VISC_DEF;
    state.diff = DIFF_DEF;
    state.colour = COLOUR_DEF;
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

void toggle_pause(UI_State &state) {
    state.paused = !state.paused;
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
        fmod(timesteps * HUE_SHIFT_FACTOR, MAX_HUE),
        SATURATION, 
        BRIGHTNESS
    );
}

void normalise_window_size() {
    if (GetScreenWidth() != (WIDTH + GUI_WIDTH) || GetScreenHeight() != HEIGHT) {
        SetWindowSize(WIDTH + GUI_WIDTH, HEIGHT);
    }
}

bool mouse_is_in_window(Vector2 mouse_pos, int scale) {
    return (mouse_pos.x >= 0 && mouse_pos.x <= WIDTH - scale &&
            mouse_pos.y >= 0 && mouse_pos.y <= HEIGHT - scale);
}

void render(UI_State &state, Fluid &fluid) {
    if (state.rendering_field < 0) {
        state.rendering_field = (FieldType) (FIELDS_COUNT-1);
    }

    if (state.blend) BeginBlendMode(BLEND_ADDITIVE);

    if (state.rendering_field == DENSITY) { 
        fluid.render_density(state.scale, state.colour, PARTICLE_SIZE_DEF);

    } else if (state.rendering_field == VELOCITY) {
        fluid.render_velocity(state.scale, state.colour, VELOCITY_ARROW_LENGTH_DEF);

    } else if (state.rendering_field == VORTICITY) {
        fluid.render_vorticity(state.scale, state.colour, PARTICLE_SIZE_DEF); 
    }

    if (state.blend) EndBlendMode();
}

void get_ui_input(UI_State &state, Fluid &fluid) {
    DrawText("Viscosity", 920, 50, STD_FONT_SIZE, WHITE);
    GuiSlider((Rectangle) { 920, 80, 200, 25 }, "0", "0.1", &state.visc, 0.0, 0.1);
    fluid.viscosity = state.visc;

    DrawText("Diffusion", 920, 140, STD_FONT_SIZE, WHITE);
    GuiSlider((Rectangle) { 920, 170, 200, 25 }, "0", "0.01", &state.diff, 0.0, 0.01);
    fluid.diffusion = state.diff;

    DrawText("Colour", 920, 250, STD_FONT_SIZE, WHITE);
    GuiColorPicker((Rectangle) { 920, 280, 200, 200 }, NULL, &state.colour);

    DrawText("Blending", 920, 550, STD_FONT_SIZE, WHITE);
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
        Vector2 mouse_travelled = state.mouse - state.mouse_prev;

        fluid.add_density(state.mouse.x / state.scale, 
                            state.mouse.y / state.scale, Vector2Length(mouse_travelled));

        fluid.add_velocity(state.mouse.x / state.scale, 
                            state.mouse.y / state.scale, mouse_travelled.x, mouse_travelled.y); 
    }

    if (IsKeyPressed(KEY_C)) fluid.reset();

    if (IsKeyPressed(KEY_R)) {
        state.resolution_index = resolution_next_index(state.resolution_index);
        state.nx = resolutions[state.resolution_index];
        state.scale = WIDTH / state.nx;
        fluid = Fluid(state.nx, state.diff, state.visc, 1.0/FPS);
    }

    if (IsKeyPressed(KEY_RIGHT)) {
        state.rendering_field = field_type_next(state.rendering_field);

    } else if (IsKeyPressed(KEY_LEFT)) {
        state.rendering_field = field_type_prev(state.rendering_field);
    }

    if (IsKeyPressed(KEY_ESCAPE)) {
        toggle_pause(state);
    }
}

int main() {
    srand(time(NULL));

    UI_State state;
    ui_state_init(state);

    Fluid fluid(state.nx, state.diff, state.visc, 1.0/FPS);

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(WIDTH + GUI_WIDTH, HEIGHT, "Diffuser");
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
        handle_input(state, fluid);
        draw_info(state);
 
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
