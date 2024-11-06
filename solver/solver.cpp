#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#include "solver.hpp"

// Checks if two doubles are equal, with tolerance = eps
inline bool double_eq(double a, double b, double eps = 1e-6) {
    return fabs(a - b) < eps;
}

// Clamps value to range [range_start, range_end], making sure it can never be out of that range
inline double constrain(double value, double range_start, double range_end) {
    return (value > range_end) ? range_end : (value < range_start) ? range_start : value;
}

// Translate value from range [cur_min, cur_max] to [res_min, res_max]. For example, 0.5 tranlated from [0, 1] to [0, 10] is 5.
double translate_value(double value, double cur_min, double cur_max, double res_min, double res_max) {
    value = constrain(value, cur_min, cur_max); 

    double cur_range = cur_max - cur_min;
    double res_range = res_max - res_min;

    double normalised_value = (value - cur_min) / cur_range;

    return res_min + normalised_value * res_range;
}

// Perform binlinear intepolation on point (x, y) with respect to cell (cell_x, cell_y)
double bilinear_interpolation(Scalar_Field &field, double x, double y, int cell_x, int cell_y) {
    double weight_right_x = x - cell_x;
    double weight_left_x = 1.0 - weight_right_x;
    double weight_top_y = y - cell_y;
    double weight_bottom_y = 1.0 - weight_top_y;

    return weight_left_x * (weight_bottom_y * field[cell_x][cell_y] + weight_top_y * field[cell_x][cell_y+1]) + 
           weight_right_x * (weight_bottom_y * field[cell_x+1][cell_y] + weight_top_y * field[cell_x+1][cell_y+1]);
}

// Set transparency (alpha) of colour based on value in the range [min, max]
Color get_colour_with_alpha(Color base, double val, double min, double max) { 
    base.a = translate_value(val, min, max, 0, 255);
    return base;
}

// Get a colour with a "temperature" (brightness) corresponding to the value in the range [min, max]
Color get_colour_temp(Color color, double val, double min, double max) {
    Vector3 complement_hsv = ColorToHSV(color);
    complement_hsv.x = ((int) complement_hsv.x + 180) % 360;
    Color complement = ColorFromHSV(complement_hsv.x, complement_hsv.y, complement_hsv.z); 

    return (Color){
        (unsigned char) translate_value(val, min, max, complement.r, color.r),
        (unsigned char) translate_value(val, min, max, complement.g, color.g),
        (unsigned char) translate_value(val, min, max, complement.b, color.b),
        (unsigned char) translate_value(abs(val), 0, max, 0, 255), 
    };
}

Scalar_Field::Scalar_Field(int _rows, int _cols) {
    rows = _rows;
    cols = _cols;
    elements = std::vector(_rows, std::vector<double>(_cols, 0));
}

std::vector<double> Scalar_Field::operator[](int index) const {
    return elements[index];
}

std::vector<double>& Scalar_Field::operator[](int index) {
    return elements[index];
}

Scalar_Field Scalar_Field::operator+(Scalar_Field b) {
    Scalar_Field result(rows, cols);

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            result[i][j] = elements[i][j] + b[i][j];
        }
    }

    return result;
}

Scalar_Field Scalar_Field::operator-(Scalar_Field b) {
    Scalar_Field result(rows, cols);

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            result[i][j] = elements[i][j] - b[i][j];
        }
    }

    return result;
}

void Scalar_Field::reset() {
    for (auto& vec : elements) std::fill(vec.begin(), vec.end(), 0);
}

Fluid::Fluid(int _nx, double _diff, double _visc, double _dt) {
    nx = _nx;
    diff = _diff;
    visc = _visc;
    dt = _dt;

    density_temp = Scalar_Field(nx, nx);
    density = Scalar_Field(nx, nx);

    vel = Vector_Field (
        Scalar_Field(nx, nx),
        Scalar_Field(nx, nx)
    );

    vel0 = Vector_Field (
        Scalar_Field(nx, nx),
        Scalar_Field(nx, nx)
    );
}

Vector_Field::Vector_Field(Scalar_Field _x, Scalar_Field _y) {
    x = _x;
    y = _y;
}

Vector_Field Vector_Field::operator+(Vector_Field b) {
    return Vector_Field(x + b.x, y + b.y);
}

Vector_Field Vector_Field::operator-(Vector_Field b) {
    return Vector_Field(x - b.x, y - b.y);
}

// https://en.wikipedia.org/wiki/Gradient
Vector_Field Vector_Field::grad(double h) {
    Vector_Field result(x, y);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.rows - 1; j++) {
            result.x[i][j] = 0.5 * (x[i+1][j] - x[i-1][j]) * h;
            result.y[i][j] = 0.5 * (x[i][j+1] - x[i][j-1]) * h;
        }
    }

    return result;
}

// https://en.wikipedia.org/wiki/Divergence
Scalar_Field Vector_Field::div(double h) {
    Scalar_Field result(x.rows, x.cols);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.cols - 1; j++) {
            result[i][j] = -0.5 * (
                    x[i+1][j] -
                    x[i-1][j] +
                    y[i][j+1] -
                    y[i][j-1]

                ) / h;
        }
    }

    return result;
}

// https://scicomp.stackexchange.com/questions/21915/discrete-definitions-of-curl-nabla-times-f
Scalar_Field Vector_Field::curl(double h) {
    Scalar_Field result(x.rows, x.cols);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.rows - 1; j++) {
            double dv_dx = (y[i+1][j] - y[i-1][j]) / (2*h);
            double du_dy = (x[i][j+1] - x[i][j-1]) / (2*h);

            result[i][j] = dv_dx - du_dy;
        }
    }

    return result;
}

// Set all boundary cells in a field x to bnd_val
void Fluid::set_bounds(Scalar_Field &x, double bnd_val) {
    for (int i = 0; i < nx; i++) {
        x[i][0] = x[i][nx-1] = bnd_val;
        x[0][1] = x[nx-1][i] = bnd_val;
    }
}

// Use the Gauss-Seidel method to solve a linear system of equations, whose coefficients are in the scalar fields x and x0.
void Fluid::gauss_seidel_solve(Scalar_Field &x, Scalar_Field &x0, double a, int iter) {
    for (int k = 0; k < iter; k++) {
        for (int i = 1; i < x.rows - 1; i++) {
            for (int j = 1; j < x.cols - 1; j++) {
                x[i][j] =  (x0[i][j] + a *
                               (x[i+1][j] +
                                x[i-1][j] +
                                x[i][j+1] +
                                x[i][j-1] )) * 1.0/(4*a + 1);
            }
        }

        set_bounds(x);
    }
}

// Apply advection onto the density using dens0, vel_x, and vel_y
void Fluid::advect(Scalar_Field &dens, Scalar_Field &dens0, Scalar_Field &vel_x, Scalar_Field &vel_y) {
    int cell_x, cell_y;
    double dtx = dt * (nx - 2);
    double x, y;

    for (int i = 1; i < nx - 1; i++) { 
        for (int j = 1; j < nx - 1; j++) {
            x = constrain(i - dtx * vel_x[i][j], 0.5, nx + 0.5);
            cell_x = constrain(floorf(x), 0, nx-2);

            y = constrain(j - dtx * vel_y[i][j], 0.5, nx + 0.5);
            cell_y = constrain(floorf(y), 0, nx-2); 
                
            dens[i][j] = bilinear_interpolation(dens0, x, y, cell_x, cell_y);
        }
    }

    set_bounds(dens);
}

// Project the velocity by applying u = u0 − ∇p
void Fluid::project_velocity() {
    vel.x.reset();
    vel.y = vel0.div(nx);

    set_bounds(vel.y);

    gauss_seidel_solve(vel.x, vel.y, 1);
    vel0 = vel0 - vel.grad(nx);

    set_bounds(vel0.x);
    set_bounds(vel0.y);
}

// DIffuse the velocity by solving the diffusion equation d(vel)/dt = (visc)(∇^2)(vel)
void Fluid::diffuse_velocity() {
    double a = dt * visc * nx * nx;

    gauss_seidel_solve(vel0.x, vel.x, a);
    gauss_seidel_solve(vel0.y, vel.y, a);
}

// DIffuse the density by solving the diffusion equation d(dens)/dt = (diff)(∇^2)(dens)
void Fluid::diffuse_density() {
    double a = dt * diff * nx * nx;
    gauss_seidel_solve(density_temp, density, a);
}

// Sequentically apply all the steps needed to complete one timestep of the similation
void Fluid::step() {
    diffuse_velocity();

    project_velocity();

    advect(vel.x, vel0.x, vel0.x, vel0.y);
    advect(vel.y, vel0.y, vel0.x, vel0.y);

    diffuse_density();

    advect(density, density_temp, vel.x, vel.y);
}

void Fluid::add_density(int x, int y, double amt) {
    density[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt; 
}

void Fluid::add_velocity(int x, int y, double amt_x, double amt_y) {    
    vel.x[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_x;
    vel.y[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_y;
}

void Fluid::render_density(int scale, Color colour, float particle_size) {
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            double dens = density[i][j];
            
            if (double_eq(dens, 0.0)) continue;
            DrawRectangle(i*scale, j*scale, particle_size*scale, particle_size*scale, 
                    get_colour_with_alpha(colour, dens, 0, 1)); 
        }
    }
}

void Fluid::render_velocity(int scale, Color colour, float line_length) {
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            Vector2 position = {.x = (float) (i*scale), .y = (float) (j*scale) };
            Vector2 velocity = {.x = (float) vel.x[i][j], .y = (float) vel.y[i][j]};
            Vector2 velocity_normal = Vector2Scale(Vector2Normalize(velocity), 4*line_length);

            if (Vector2Equals(velocity_normal, Vector2Zero())) DrawPixelV(position, colour);
            DrawLineEx(position, Vector2Add(velocity_normal, position), 3, colour);
        }
    }
}

void Fluid::render_vorticity(int scale, Color colour, float particle_size) {
    Scalar_Field vorticity = vel.curl();
    set_bounds(vorticity);

    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < nx - 2; j++) {
            double vort = vorticity[i][j];
            if (double_eq(vort, 0.0)) continue;

            double colour_rng = visc == 0 ? 0.02 : 0.05;
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

void Fluid::reset() {
    density_temp.reset();
    density.reset();
    vel.x.reset();
    vel.y.reset();
    vel0.x.reset();
    vel0.y.reset();
}

void Fluid::debug_print_density(FILE *f) {
    for (int i = 0; i < nx; i++) {
        fprintf(f, "{ ");
        for (int j = 0; j < nx; j++) {
            fprintf(f, "%lf ", density[i][j]);
        }
        fprintf(f, "}\n");
    } 
}
