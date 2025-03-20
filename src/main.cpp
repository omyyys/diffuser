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
constexpr int HELP_ICON_FONT_SIZE = 15;
constexpr int MAIN_MENU_TITLE_SIZE = 100;
constexpr int MAIN_MENU_STD_FONT_SIZE = 40;

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

std::string resolution_to_description(int resolution) {
    switch (resolution) {
        case 80:
            return "very low";
        case 100:
            return "low";
        case 160:
            return "medium";
        case 200:
            return "high";
        default:
            return "unknown resolution";
    }
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

typedef enum {
    APP_STATE_MAIN_MENU,
    APP_STATE_SIMULATION

} App_State;

typedef struct {
    bool diffusion_help;
    bool viscosity_help;
    bool blending_help;
    bool spawn_interval_help;
    
} Help_Icons_State;

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

void gui_set_font_size(int size) {
    GuiSetStyle(DEFAULT, TEXT_SIZE, size);
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

void show_help_message(Help_Icons_State &help_state) {
    std::string message;

    if (help_state.viscosity_help) {
        message = "Viscosity is the resistance of a fluid to flow, \n\n and can roughly be characterised as its thickness. \n\n A higher viscosity will make the fluid flow slower. \n\n For example, honey has a higher viscosity than water.";

    } else if (help_state.diffusion_help) {
        message = "Diffusion is the process which causes the fluid to \n\n spread out from areas of high concentration to areas of \n\n low concentration. A higher diffusion will cause the fluid \n\n to spread out more quickly.";

    } else if (help_state.blending_help) {
        message = "Additive blending adds the RGB components of any \n\n overlapping colours. This has the effect of making parts \n\n of the fluid with a high concentration appear brighter.";

    } else if (help_state.spawn_interval_help) {
        message = "The spawn interval is the time in frames between each \n\n splat of fluid, if auto mode is selected. For example, setting \n\n A lower interval will cause the fluid to spawn more \n\n frequently.";

    } else {
        return;
    }
    
    int result = GuiMessageBox((Rectangle){ 100, 250, 600, 300 }, "Help", message.c_str(), "OK");
    if (result != -1) help_state = {0};
}

void draw_help_icon(float x_pos, float y_pos, float size, Help_Icons_State &help_state, bool &help_state_field) {
    gui_set_font_size(HELP_ICON_FONT_SIZE);
    if (GuiButton((Rectangle){ x_pos, y_pos, size, size }, "?")) {
        help_state = {0};
        help_state_field = true;
    }
    gui_set_font_size(STD_FONT_SIZE);
}

void get_ui_input(UI_State &state, Help_Icons_State &help_state, Fluid &fluid) {
    DrawText("Viscosity", 920, 50, STD_FONT_SIZE, WHITE);

    draw_help_icon(1020, 50, STD_FONT_SIZE, help_state, help_state.viscosity_help);

    GuiSlider((Rectangle) { 920, 80, 200, 25 }, "0", "0.1", &state.visc, 0.0, 0.1);
    fluid.viscosity = state.visc;

    DrawText("Diffusion", 920, 140, STD_FONT_SIZE, WHITE);

    draw_help_icon(1020, 140, STD_FONT_SIZE, help_state, help_state.diffusion_help);

    GuiSlider((Rectangle) { 920, 170, 200, 25 }, "0", "0.01", &state.diff, 0.0, 0.01);
    fluid.diffusion = state.diff;

    DrawText("Colour", 920, 250, STD_FONT_SIZE, WHITE);
    GuiColorPicker((Rectangle) { 920, 280, 200, 200 }, NULL, &state.colour);

    DrawText("Blending", 920, 550, STD_FONT_SIZE, WHITE);

    draw_help_icon(1015, 550, STD_FONT_SIZE, help_state, help_state.blending_help);

    GuiCheckBox((Rectangle) { 1120, 550, 20, 20 }, NULL, &state.blend);

    DrawText("Auto Colours", 920, 600, STD_FONT_SIZE, WHITE);
    GuiCheckBox((Rectangle) { 1120, 600, 20, 20 }, NULL, &state.colour_change);

    DrawText("Auto Spawn", 920, 650, STD_FONT_SIZE, WHITE);
    GuiCheckBox((Rectangle) { 1120, 650, 20, 20 }, NULL, &state.random_spawn);

    if (state.random_spawn) {
        float random_spawn_period_float = state.random_splat_period;
        DrawText("Spawn Interval", 920, 700, STD_FONT_SIZE, WHITE);

        draw_help_icon(1085, 700, STD_FONT_SIZE, help_state, help_state.spawn_interval_help);

        GuiSlider((Rectangle) { 920, 730, 200, 25 }, "30", "240", &random_spawn_period_float, 30, 240);
        state.random_splat_period = random_spawn_period_float;
    }

    show_help_message(help_state);
}

void draw_info(UI_State state) {
    DrawText(TextFormat("grid resolution: %s", resolution_to_description(state.nx).c_str()), 20, 20, STD_FONT_SIZE, WHITE);

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

void simulation_update(UI_State &ui_state, Help_Icons_State &help_state, Fluid &fluid) {
    ui_state.mouse_prev = ui_state.mouse;
    ui_state.mouse = GetMousePosition();

    if (!ui_state.paused) {
        fluid.step();
        fluid.fade_density(BASE_FADE + ui_state.visc * 0.001);
        ui_state.timesteps++;
    }

    if (ui_state.colour_change && ui_state.timesteps % COLOUR_CHANGE_PERIOD == 0) {
        ui_state.colour = get_colour(ui_state.timesteps); 
    }

    if (ui_state.random_spawn && ui_state.timesteps % ui_state.random_splat_period == 0) {
        random_spawn(fluid, ui_state.nx);
    }

    render(ui_state, fluid);
    get_ui_input(ui_state, help_state, fluid);
    handle_input(ui_state, fluid);
    draw_info(ui_state);
}

void main_menu_update(App_State &app_state, Texture bg_texture) {
    DrawTexture(bg_texture, GetScreenWidth() - bg_texture.width, GetScreenHeight() - bg_texture.height, WHITE);
    DrawText("diffuser", 100, 100, MAIN_MENU_TITLE_SIZE, WHITE);

    if (GuiButton((Rectangle){ 500, 300, 500, 100 }, "Start")) {
        gui_set_font_size(STD_FONT_SIZE);
        app_state = APP_STATE_SIMULATION;
    }

    if (GuiButton((Rectangle){ 500, 450, 500, 80 }, "Exit")) exit(0);
}

int main() {
    srand(time(NULL));

    UI_State ui_state;
    ui_state_init(ui_state);

    App_State app_state = APP_STATE_MAIN_MENU;

    Help_Icons_State help_icons_state = {0};

    Fluid fluid(ui_state.nx, ui_state.diff, ui_state.visc, 1.0/FPS);

    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(WIDTH + GUI_WIDTH, HEIGHT, "Diffuser");
    gui_set_font_size(MAIN_MENU_STD_FONT_SIZE);
    SetTargetFPS(FPS);
    SetExitKey(KEY_NULL);

    Texture2D main_menu_bg = LoadTexture("menu_bg.png");

    while (!WindowShouldClose()) {
        normalise_window_size();

        BeginDrawing();
        ClearBackground(BLACK);

        if (app_state == APP_STATE_MAIN_MENU) {
            main_menu_update(app_state, main_menu_bg);

        } else if (app_state == APP_STATE_SIMULATION) {
            simulation_update(ui_state, help_icons_state, fluid);
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}
