#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#include "solver.hpp"

inline bool double_eq(double a, double b, double eps = 1e-6) {
    return fabs(a - b) < eps;
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

Color get_colour_w_alpha(Color base, double val, double min, double max) { 
    base.a = lin_interpolation(val, min, max, 0, 255);
    return base;
}

Color get_colour_temp(Color color, double val, double min, double max) {
    Vector3 complement_hsv = ColorToHSV(color);
    complement_hsv.x = ((int) complement_hsv.x + 180) % 360;
    Color complement = ColorFromHSV(complement_hsv.x, complement_hsv.y, complement_hsv.z); 

    return (Color){
        (unsigned char) lin_interpolation(val, min, max, complement.r, color.r),
        (unsigned char) lin_interpolation(val, min, max, complement.g, color.g),
        (unsigned char) lin_interpolation(val, min, max, complement.b, color.b),
        (unsigned char) lin_interpolation(abs(val), 0, max, 0, 255), 
    };
}

Scalar_Field::Scalar_Field(int _rows, int _cols) {
    rows = _rows;
    cols = _cols;
    elements = std::vector(_rows, std::vector<double>(_cols, 0));
}

std::vector<double> Scalar_Field::operator[] (int index) const {
    return elements[index];
}

std::vector<double>& Scalar_Field::operator[] (int index) {
    return elements[index];
}

Fluid::Fluid(int _nx, double _diff, double _visc, double _dt) {
    nx = _nx;
    diff = _diff;
    visc = _visc;
    dt = _dt;

    s = Scalar_Field(nx, nx);
    density = Scalar_Field(nx, nx);

    v = {
        .x = Scalar_Field(nx, nx),
        .y = Scalar_Field(nx, nx)
    };

    v0 = {
        .x = Scalar_Field(nx, nx),
        .y = Scalar_Field(nx, nx)
    };
}

void Fluid::set_bnd(Scalar_Field &x, double bnd_val) {
    for (int _ = 0; _ < nx; _++) {
        for (int i = 0; i < nx; i++) {
            x[i][0] = bnd_val;
            x[i][nx-1] = bnd_val;
        }
    }

    for (int _ = 0; _ < nx; _++) {
        for (int j = 0; j < nx; j++) {
            x[0][j] = bnd_val;
            x[nx-1][j] = bnd_val;
        }
    }
}

void Fluid::gauss_seidel(Scalar_Field &x, Scalar_Field &x0, double a, int iter) {
    double c = 4*a + 1; 

    for (int k = 0; k < iter; k++) {
        for (int i = 1; i < x.rows - 1; i++) {
            for (int j = 1; j < x.cols - 1; j++) {
                x[i][j] =  (x0[i][j] + a *
                           (x[i+1][j] +
                            x[i-1][j] +
                            x[i][j+1] +
                            x[i][j-1] )) * (1.0/c);
            }
        }
        set_bnd(x);
    }
}

Scalar_Field Fluid::div(Vector_Field field, double h) {
    Scalar_Field result(field.x.rows, field.x.cols);

    for (int i = 1; i < field.x.rows - 1; i++) {
        for (int j = 1; j < field.y.cols - 1; j++) {
            result[i][j] = -0.5 * (
                    field.x[i+1][j] -
                    field.x[i-1][j] +
                    field.y[i][j+1] -
                    field.y[i][j-1]

                ) / h;
        }
    }

    return result;
}

// https://scicomp.stackexchange.com/questions/21915/discrete-definitions-of-curl-nabla-times-f
Scalar_Field Fluid::curl(Vector_Field field, double h) {
    Scalar_Field result(field.x.rows, field.x.cols);

    for (int i = 1; i < field.x.rows - 1; i++) {
        for (int j = 1; j < field.y.rows - 1; j++) {
            double dv_dx = (field.y[i+1][j] - field.y[i-1][j]) / (2*h);
            double du_dy = (field.x[i][j+1] - field.x[i][j]) / (2*h);

            result[i][j] = dv_dx - du_dy;
        }
    }

    return result;
}

void Fluid::advect(Scalar_Field &d, Scalar_Field &d0, Scalar_Field &v_x, Scalar_Field &v_y) {
    int ix, iy;

    double dtx = dt * (nx - 2);
    double dty = dt * (nx - 2);
    
    double s0, s1, t0, t1;
    double x, y;

    for (double j = 1; j < nx - 1; j++) { 
        for (double i = 1; i < nx - 1; i++) {
            x = constrain(i - dtx * v_x[i][j], 0.5, nx + 0.5);
            ix = constrain(floorf(x), 0, nx-2);

            y = constrain(j - dty * v_y[i][j], 0.5, nx + 0.5);
            iy = constrain(floorf(y), 0, nx-2); 
                
            s1 = x - ix; 
            s0 = 1.0 - s1; 
            t1 = y - iy; 
            t0 = 1.0 - t1;
                
            d[i][j] = 
                s0 * (t0 * d0[ix][iy] + t1 * d0[ix][iy+1]) + 
                s1 * (t0 * d0[ix+1][iy] + t1 * d0[ix+1][iy+1]);
        }
    }

    set_bnd(d);
}

void Fluid::project_velocity() {
    for (auto& vec : v.x.elements) std::fill(vec.begin(), vec.end(), 0);
    v.y = div(v0, nx);

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
    double a = dt * visc * nx * nx;

    gauss_seidel(v0.x, v.x, a);
    gauss_seidel(v0.y, v.y, a);
}

void Fluid::diffuse_density() {
    double a = dt * diff * nx * nx;

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

void Fluid::add_density(int x, int y, double amt) {
    density[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt; 
}

void Fluid::add_velocity(int x, int y, double amt_x, double amt_y) {    
    v.x[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_x;
    v.y[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_y;
}

void Fluid::render_density(int scale, Color colour, float particle_size) {
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            double dens = density[i][j];
            DrawRectangle(i*scale, j*scale, particle_size*scale, particle_size*scale, 
                          get_colour_w_alpha(colour, dens, 0, 1));
        }
    }
}

void Fluid::render_velocity(int scale, Color colour, float line_length) {
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            Vector2 position = {.x = (float) (i*scale), .y = (float) (j*scale) };
            Vector2 velocity = {.x = (float) v.x[i][j], .y = (float) v.y[i][j]};
            Vector2 velocity_normal = Vector2Scale(Vector2Normalize(velocity), 4*line_length);

            if (Vector2Equals(velocity_normal, Vector2Zero())) DrawPixelV(position, colour);
            DrawLineEx(position, Vector2Add(velocity_normal, position), 3, colour);
        }
    }
}

void Fluid::render_vorticity(int scale, Color colour, float particle_size) {
    Scalar_Field vorticity = curl(v);
    set_bnd(vorticity);

    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < nx - 2; j++) {
            double vort = vorticity[i][j];
            double colour_rng = visc == 0 ? 0.01 : 0.0005;
            DrawRectangle(i*scale, j*scale, (particle_size/4)*scale, (particle_size/4)*scale, 
                          get_colour_temp(colour, vort, -colour_rng, colour_rng));
        }
    }
}

void Fluid::fade_density(double fade) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            density[i][j] = constrain(density[i][j] - fade, 0, 1);
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
    for (auto& vec : s.elements) std::fill(vec.begin(), vec.end(), 0);
    for (auto& vec : density.elements) std::fill(vec.begin(), vec.end(), 0);

    for (auto& vec : v.x.elements) std::fill(vec.begin(), vec.end(), 0);
    for (auto& vec : v.y.elements) std::fill(vec.begin(), vec.end(), 0);
    for (auto& vec : v0.x.elements) std::fill(vec.begin(), vec.end(), 0);
    for (auto& vec : v0.y.elements) std::fill(vec.begin(), vec.end(), 0);
}
