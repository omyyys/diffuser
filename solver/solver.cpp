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

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max) {
    value = constrain(value, cur_min, cur_max); 

    double cur_range = cur_max - cur_min;
    double res_range = res_max - res_min;

    double normalised_value = (value - cur_min) / cur_range;

    return res_min + normalised_value * res_range;
}

Color get_colour(Color base, double dens) { 
    base.a = lin_interpolation(dens, 0.0, 1.0, 0, 255);
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

    v = {
        .x = std::vector<double>(nx*nx),
        .y = std::vector<double>(nx*nx) 
    };

    v0 = {
        .x = std::vector<double>(nx*nx),
        .y = std::vector<double>(nx*nx) 
    };
}

void Fluid::set_bnd(Grid &x) {
    for (int k = 1; k < nx - 1; k++) {
        for (int i = 1; i < nx - 1; i++) {
            x[IX(i, 0)] = 0;
            x[IX(i, nx-1)] = 0;
        }
    }

    for (int k = 1; k < nx - 1; k++) {
        for (int j = 1; j < nx - 1; j++) {
            x[IX(0, j)] = 0;
            x[IX(nx-1, j)] = 0;
        }
    }
    
    x[IX(0, 0)] = 0;
    x[IX(0, nx-1)] = 0;
    x[IX(nx-1, 0)] = 0;
    x[IX(nx-1, nx-1)] = 0;
}

void Fluid::gauss_seidel(Grid &x, Grid &x0, double a, double c) {
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
        set_bnd(x);
    }
}

void Fluid::advect(Grid &d, Grid &d0, Grid &velocX, Grid &velocY) {
    double i0, i1, j0, j1;

    double dtx = dt * (nx - 2);
    double dty = dt * (nx - 2);
    
    double s0, s1, t0, t1;
    double x, y;
    
    double nxdouble = nx;
    double idouble, jdouble;

    int i, j;
    
    for (j = 1, jdouble = 1; j < nx - 1; j++, jdouble++) { 
        for (i = 1, idouble = 1; i < nx - 1; i++, idouble++) {
            x = idouble - dtx * velocX[IX(i, j)]; 
            y = jdouble - dty * velocY[IX(i, j)];
                
            if (x < 0.5) x = 0.5;
            x = constrain(x, 0, nxdouble + 0.5);
            i0 = floorf(x); 
            i1 = i0 + 1.0;

            if (y < 0.5) y = 0.5; 
            y = constrain(y, 0, nxdouble + 0.5);
            j0 = floorf(y);
            j1 = j0 + 1.0; 
                
            s1 = x - i0; 
            s0 = 1.0 - s1; 
            t1 = y - j0; 
            t0 = 1.0 - t1;
                
            d[IX(i, j)] = 
                s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + 
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }

    set_bnd(d);
}

void Fluid::project(Grid &velocX, Grid &velocY, Grid &p, Grid &div)
{
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            div[IX(i, j)] = -0.5 * (
                    velocX[IX(i+1, j)] -
                    velocX[IX(i-1, j)] +
                    velocY[IX(i, j+1)] -
                    velocY[IX(i, j-1)]
                ) / nx;

            p[IX(i, j)] = 0;
        }
    }

    set_bnd(div); 
    set_bnd(p);
    gauss_seidel(p, div, 1, 6);
    
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]) * nx;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]) * nx;
        }
    }

    set_bnd(velocX);
    set_bnd(velocY);
}

void Fluid::diffuse_velocity() {
    double a = dt * visc * nx * nx;

    gauss_seidel(v0.x, v.x, a, 4*a + 1);
    gauss_seidel(v0.y, v.y, a, 4*a + 1);
}

void Fluid::diffuse_density() {
    double a = dt * diff * nx * nx;

    gauss_seidel(s, density, a, 4*a + 1);
}

void Fluid::step() {
    diffuse_velocity();

    project(v0.x, v0.y, v.x, v.y);

    advect(v.x, v0.x, v0.x, v0.y);
    advect(v.y, v0.y, v0.x, v0.y);

    diffuse_density();
    advect(density, s, v.x, v.y);
}

void Fluid::add_density(int x, int y, double amount) {
    density[IX(x, y)] = amount + density[IX(x, y)]; 
}

void Fluid::add_velocity(int x, int y, double amountX, double amountY) {
    int index = IX(x, y);
    
    v.x[index] = amountX + v.x[index];
    v.y[index] = amountY + v.y[index];
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

    std::fill(v.x.begin(), v.y.end(), 0);
    std::fill(v.y.begin(), v.y.end(), 0);
    std::fill(v0.x.begin(), v0.x.end(), 0);
    std::fill(v0.y.begin(), v0.y.end(), 0);
}
