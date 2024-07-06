#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "../raylib/include/raylib.h"
#include "solver.hpp"

double avg(int count, ...) {
    double result = 0.0;

    va_list args;
    va_start(args, count);

    for (int i = 0; i < count; i++) {
        result += va_arg(args, double);
    }

    return result / count;
}

double constrain(double value, double range_start, double range_end) {
    value = (value > range_end) ? range_end : value;
    value = (value < range_start) ? range_start : value;

    return value;
}

double translate_range(double value, double cur_min, double cur_max, double res_min, double res_max) {
    value = constrain(value, cur_min, cur_max); 

    double cur_range = cur_max - cur_min;
    double res_range = res_max - res_min;

    double normalised_value = (value - cur_min) / cur_range;

    return res_min + normalised_value * res_range;
}

Color get_colour(Color base, double dens) { 
    base.a = translate_range(dens, 0.0, 1.0, 0, 255);
    return base;
}

bool double_eq(double a, double b, double eps = 1e-5) {
    return fabs(a - b) < eps;
}

Fluid::Fluid(int _nx, double _diff, double _visc, double _dt) {
    nx = _nx;
    diff = _diff;
    visc = _visc;
    dt = _dt;
    iter = 5;

    s = std::vector<double>(nx*nx); 
    density = std::vector<double>(nx*nx);

    Vx = std::vector<double>(nx*nx);
    Vy = std::vector<double>(nx*nx);
    Vx0 = std::vector<double>(nx*nx);
    Vy0 = std::vector<double>(nx*nx);
}

void Fluid::set_bnd(int b, Grid &x) {
    for (int k = 1; k < nx - 1; k++) {
        for (int i = 1; i < nx - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, nx-1)] = b == 2 ? -x[IX(i, nx-2)] : x[IX(i, nx-2)];
        }
    }
    for (int k = 1; k < nx - 1; k++) {
        for (int j = 1; j < nx - 1; j++) {
            x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
            x[IX(nx-1, j)] = b == 1 ? -x[IX(nx-2, j)] : x[IX(nx-2, j)];
        }
    }
    
    x[IX(0, 0)] = avg(2, x[IX(1, 0)], x[IX(0, 1)]);
    x[IX(0, nx-1)] = avg(2, x[IX(1, nx-1)], x[IX(0, nx-2)]);
    x[IX(nx-1, 0)] = avg(2, x[IX(nx-2, 0)], x[IX(nx-1, 1)]);
    x[IX(nx-1, nx-1)] = avg(2, x[IX(nx-2, nx-1)], x[IX(nx-1, nx-2)]);
}

void Fluid::lin_solve(int b, Grid &x, Grid &x0, double a, double c) {
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < nx - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                x[IX(i, j)] =  (x0[IX(i, j)] + a *
                               (x[IX(i+1, j)] +
                                x[IX(i-1, j)] +
                                x[IX(i, j+1)] +
                                x[IX(i, j-1)] )) * 
                                (1.0/c);
            }
        }
        set_bnd(b, x);
    }
}

void Fluid::diffuse(int b, double visc_or_diff, Grid &x, Grid &x0)
{
    double a = dt * visc_or_diff * nx * nx;
    lin_solve(b, x, x0, a, 4*a + 1);
}

void Fluid::advect(int b, Grid &d, Grid &d0, Grid &velocX, Grid &velocY) {
    double i0, i1, j0, j1;

    double dtx = dt * (nx - 2);
    double dty = dt * (nx - 2);
    
    double s0, s1, t0, t1;
    double tmp1, tmp2, x, y;
    
    double nxdouble = nx;
    double idouble, jdouble;

    int i, j;
    
    for (j = 1, jdouble = 1; j < nx - 1; j++, jdouble++) { 
        for (i = 1, idouble = 1; i < nx - 1; i++, idouble++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];

            x    = idouble - tmp1; 
            y    = jdouble - tmp2;
                
            if(x < 0.5f) x = 0.5f; 
            if(x > nxdouble + 0.5f) x = nxdouble + 0.5f; 
            i0 = floorf(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > nxdouble + 0.5f) y = nxdouble + 0.5f; 
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

    set_bnd(b, d);
}

void Fluid::project(Grid &velocX, Grid &velocY, Grid &p, Grid &div)
{
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            div[IX(i, j)] = -0.5f * (
                    velocX[IX(i+1, j)] -
                    velocX[IX(i-1, j)] +
                    velocY[IX(i, j+1)] -
                    velocY[IX(i, j-1)]
                ) / nx;

            p[IX(i, j)] = 0;
        }
    }

    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * nx;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * nx;
        }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void Fluid::step() {
    diffuse(1, visc, Vx0, Vx);
    diffuse(2, visc, Vy0, Vy);
    
    project(Vx0, Vy0, Vx, Vy);
    
    advect(1, Vx, Vx0, Vx0, Vy0);
    advect(2, Vy, Vy0, Vx0, Vy0);
    
    project(Vx, Vy, Vx0, Vy0);
    
    diffuse(0, diff, s, density);
    advect(0, density, s, Vx, Vy);
}

void Fluid::add_density(int x, int y, double amount) {
    density[IX(x, y)] = amount + density[IX(x, y)]; 
}

void Fluid::add_velocity(int x, int y, double amountX, double amountY) {
    int index = IX(x, y);
    
    Vx[index] = amountX + Vx[index];
    Vy[index] = amountY + Vy[index];
}

void Fluid::render_density(int scale, Color colour, float particle_size) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            double dens = density[IX(i, j)];
            DrawRectangle(i*scale, j*scale, particle_size*scale, particle_size*scale, 
                          get_colour(colour, dens));
        }
    }
}

void Fluid::print_density(FILE *f) {
    for (int i = 0; i < nx; i++) {
        fprintf(f, "{ ");
        for (int j = 0; j < nx; j++) {
            fprintf(f, "%lf ", density[IX(i, j)]);
        }
        fprintf(f, "}\n");
    } 
}

void Fluid::reset() {
    std::fill(s.begin(), s.end(), 0);
    std::fill(density.begin(), density.end(), 0);

    std::fill(Vx.begin(), Vy.end(), 0);
    std::fill(Vy.begin(), Vy.end(), 0);
    std::fill(Vx0.begin(), Vx0.end(), 0);
    std::fill(Vy0.begin(), Vy0.end(), 0);
}
