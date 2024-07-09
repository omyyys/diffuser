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
    iter = 20;

    s = std::vector(nx, std::vector<double>(nx, 0));
    density = std::vector(nx, std::vector<double>(nx, 0));

    v = {
        .x = std::vector(nx, std::vector<double>(nx, 0)),
        .y = std::vector(nx, std::vector<double>(nx, 0)) 
    };

    v0 = {
        .x = std::vector(nx, std::vector<double>(nx, 0)),
        .y = std::vector(nx, std::vector<double>(nx, 0)) 
    };
}

void Fluid::set_bnd(Grid &x) {
    for (int _ = 0; _ < nx; _++) {
        for (int i = 0; i < nx; i++) {
            x[i][0] = 0;
            x[i][nx-1] = 0;
        }
    }

    for (int _ = 0; _ < nx; _++) {
        for (int j = 0; j < nx; j++) {
            x[0][j] = 0;
            x[nx-1][j] = 0;
        }
    }
}

void Fluid::gauss_seidel(Grid &x, Grid &x0, double a) {
    double c = 4*a + 1; 

    for (int k = 0; k < iter; k++) {
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < nx - 1; j++) {
                x[i][j] =  ((x0[i][j]) + a *
                           (x[i+1][j] +
                            x[i-1][j] +
                            x[i][j+1] +
                            x[i][j-1] )) * (1.0/c);
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

    for (double j = 1; j < nx - 1; j++) { 
        for (double i = 1; i < nx - 1; i++) {
            x = i - dtx * velocX[i][j]; 
            y = j - dty * velocY[i][j];

            x = constrain(x, 0.5, nxdouble + 0.5);
            i0 = floorf(x); 
            i1 = i0 + 1.0;

            y = constrain(y, 0.5, nxdouble + 0.5);
            j0 = floorf(y);
            j1 = j0 + 1.0; 
                
            s1 = x - i0; 
            s0 = 1.0 - s1; 
            t1 = y - j0; 
            t0 = 1.0 - t1;
                
            d[i][j] = 
                s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) + 
                s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
        }
    }

    set_bnd(d);
}

void Fluid::project_velocity() {
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            v.y[i][j] = -0.5 * (
                    v0.x[i+1][j] -
                    v0.x[i-1][j] +
                    v0.y[i][j+1] -
                    v0.y[i][j-1]
                ) / nx;

            v.x[i][j] = 0;
        }
    }

    set_bnd(v.x); 
    set_bnd(v.y);
    gauss_seidel(v.x, v.y, 1);
    
    for (int j = 1; j < nx - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            v0.x[i][j] -= 0.5 * (v.x[i+1][j] - v.x[i-1][j]) * nx;
            v0.y[i][j] -= 0.5 * (v.x[i][j+1] - v.x[i][j-1]) * nx;
        }
    }

    set_bnd(v0.x);
    set_bnd(v0.y);
}

void Fluid::diffuse_velocity() {
    double a = dt * visc * (nx-2) * (nx-2);

    gauss_seidel(v0.x, v.x, a);
    gauss_seidel(v0.y, v.y, a);
}

void Fluid::diffuse_density() {
    double a = dt * diff * (nx-2) * (nx-2);

    gauss_seidel(s, density, a);
}

void Fluid::step() {
    diffuse_velocity();

    project_velocity();

    advect(v.x, v0.x, v0.x, v0.y);
    advect(v.y, v0.y, v0.x, v0.y);

    diffuse_density();
    advect(density, s, v.x, v.y);
}

void Fluid::add_density(int x, int y, double amount) {
    density[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amount; 
}

void Fluid::add_velocity(int x, int y, double amountX, double amountY) {    
    v.x[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amountX;
    v.y[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amountY;
}

void Fluid::render_density(int scale, Color colour, float particle_size) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            double dens = density[i][j];
            DrawRectangle(i*scale, j*scale, particle_size*scale, particle_size*scale, 
                          get_colour(colour, dens));
        }
    }
}

void Fluid::print_density(FILE *f) {
    for (int i = 0; i < nx; i++) {
        fprintf(f, "{ ");
        for (int j = 0; j < nx; j++) {
            fprintf(f, "%lf ", density[i][j]);
        }
        fprintf(f, "}\n");
    } 
}

void Fluid::reset() {
    std::fill(s.begin(), s.end(), std::vector<double>(nx));
    std::fill(density.begin(), density.end(), std::vector<double>(nx));

    std::fill(v.x.begin(), v.y.end(), std::vector<double>(nx));
    std::fill(v.y.begin(), v.y.end(), std::vector<double>(nx));
    std::fill(v0.x.begin(), v0.x.end(), std::vector<double>(nx));
    std::fill(v0.y.begin(), v0.y.end(), std::vector<double>(nx));
}
