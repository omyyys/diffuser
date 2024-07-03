#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <unistd.h>

#include "raylib/include/raylib.h"
#include "solver/solver.hpp"

#define IX(x, y) (x + y*N)

double translate_range(double value, double cur_min, double cur_max, double res_min, double res_max) {
    value = (value > cur_max) ? cur_max : value;
    value = (value < cur_min) ? cur_min : value;

    double cur_range = cur_max - cur_min;
    double res_range = res_max - res_min;

    double normalised_value = (value - cur_min) / cur_range;

    return res_min + normalised_value * res_range;
}

struct FluidCube {
    int size;
    double dt;
    double diff;
    double visc;
    
    double *s;
    double *density;
    
    double *Vx;
    double *Vy;

    double *Vx0;
    double *Vy0;
};


FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, double dt) {
    FluidCube *cube = (FluidCube*) malloc(sizeof(*cube));
    int N = size;
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = (double*) calloc(N * N, sizeof(double));
    cube->density = (double*) calloc(N * N, sizeof(double));
    
    cube->Vx = (double*) calloc(N * N, sizeof(double));
    cube->Vy = (double*) calloc(N * N, sizeof(double));
    
    cube->Vx0 = (double*) calloc(N * N, sizeof(double));
    cube->Vy0 = (double*) calloc(N * N, sizeof(double));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube) {
    free(cube->s);
    free(cube->density);
    
    free(cube->Vx);
    free(cube->Vy);
    
    free(cube->Vx0);
    free(cube->Vy0);
    
    free(cube);
}

static void set_bnd(int b, double *x, int N)
{
    for (int k = 1; k < N - 1; k++) {
        for (int i = 1; i < N - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
        }
    }
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
            x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
        }
    }
    
    x[IX(0, 0)]       = 0.5f * (x[IX(1, 0)]
                                  + x[IX(0, 1)]);

    x[IX(0, N-1)]     = 0.5f * (x[IX(1, N-1)]
                                  + x[IX(0, N-2)]);

    x[IX(N-1, 0)]     = 0.5f * (x[IX(N-2, 0)]
                                  + x[IX(N-1, 1)]);

    x[IX(N-1, N-1)]   = 0.5f * (x[IX(N-2, N-1)]
                                  + x[IX(N-1, N-2)]);
}

static void lin_solve(int b, double *x, double *x0, double a, double c, int iter, int N)
{
    double cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a*(    x[IX(i+1, j)]
                                +x[IX(i-1, j)]
                                +x[IX(i  , j+1)]
                                +x[IX(i  , j-1)]
                        )) * cRecip;
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, double *x, double *x0, double diff, double dt, int iter, int N)
{
    double a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, double *d, double *d0, double *velocX, double *velocY, double dt, int N)
{
    double i0, i1, j0, j1;

    double dtx = dt * (N - 2);
    double dty = dt * (N - 2);
    
    double s0, s1, t0, t1;
    double tmp1, tmp2, x, y;
    
    double Ndouble = N;
    double idouble, jdouble;

    int i, j;
    
    for (j = 1, jdouble = 1; j < N - 1; j++, jdouble++) { 
        for (i = 1, idouble = 1; i < N - 1; i++, idouble++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];

            x    = idouble - tmp1; 
            y    = jdouble - tmp2;
                
            if(x < 0.5f) x = 0.5f; 
            if(x > Ndouble + 0.5f) x = Ndouble + 0.5f; 
            i0 = floorf(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > Ndouble + 0.5f) y = Ndouble + 0.5f; 
            j0 = floorf(y);
            j1 = j0 + 1.0f; 
                
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
                
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
                
            d[IX(i, j)] = 
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) + 
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }

    set_bnd(b, d, N);
}

static void project(double *velocX, double *velocY, double *p, double *div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                     velocX[IX(i+1, j)]
                    -velocX[IX(i-1, j)]
                    +velocY[IX(i  , j+1)]
                    -velocY[IX(i  , j-1)]
                )/N;
            p[IX(i, j)] = 0;
        }
    }

    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                            -p[IX(i-1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                            -p[IX(i, j-1)]) * N;
        }
    }

    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

void FluidCubeStep(FluidCube *cube)
{
    int N          = cube->size;
    double visc     = cube->visc;
    double diff     = cube->diff;
    double dt       = cube->dt;
    double *Vx      = cube->Vx;
    double *Vy      = cube->Vy;
    double *Vx0     = cube->Vx0;
    double *Vy0     = cube->Vy0;
    double *s       = cube->s;
    double *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    
    project(Vx0, Vy0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);
    
    project(Vx, Vy, Vx0, Vy0, 4, N);
    
    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, dt, N);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, double amount) {
    int N = cube->size;
    cube->density[IX(x, y)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, double amountX, double amountY) {
    int N = cube->size;
    int index = IX(x, y);
    
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
}

void fluid_field_print(double *field, int N) {
    for (int i = 0; i < N; i++) {
        printf("{ ");
        for (int j = 0; j < N; j++) {
            printf("%f ", field[IX(i, j)]);
        }
        printf("}\n");
    } 
}

bool double_eq(double a, double b, double eps = 1e-5) {
    return fabs(a - b) < eps;
}

Color get_colour(double dens) {
    unsigned char dens_alpha = translate_range(dens, 0.0, 1.0, 0, 255);
    return (Color) {
        .r = dens_alpha,
        .g = dens_alpha,
        .b = 255,
        .a = dens_alpha
    };
}

void render_density(FluidCube *fluid, int N, int scale) {
    BeginDrawing();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dens = fluid->density[IX(i, j)];
            DrawRectangle(i*scale, j*scale, 3*scale, 3*scale, get_colour(dens));
        }
    }

    EndDrawing();
}

int main() {
    const int scale = 10;
    const int N = 60;
    const int fps = 240;

    FluidCube *fluid = FluidCubeCreate(N, 0, 0, 1.0/fps);

    InitWindow(N*scale, N*scale, "test");
    SetTargetFPS(fps);

    while (!WindowShouldClose()) {
        ClearBackground(BLACK);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            FluidCubeAddDensity(fluid, GetMouseX() / scale, GetMouseY() / scale, 1.0);
            FluidCubeAddVelocity(fluid, GetMouseX() / scale, GetMouseY() / scale, 1.0, 1.0);
        }

        FluidCubeStep(fluid);
        render_density(fluid, N, scale);

        DrawFPS(GetScreenWidth() - 100, 10);
    }

    CloseWindow();

    return 0;
}
