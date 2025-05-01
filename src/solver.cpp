#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "../raylib/include/raylib.h"
#include "../raylib/include/raymath.h"

#include "../include/solver.hpp"

constexpr int VELOCITY_LINE_THICKNESS = 3;

inline bool double_eq(double a, double b, double eps = 1e-6) {
    return fabs(a - b) <= eps * fmax(fabs(a), fabs(b));
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

double bilinear_interpolation(ScalarField &field, double x, double y, int cell_x, int cell_y) {
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

ScalarField::ScalarField(int _rows, int _cols) {
    rows = _rows;
    cols = _cols;
    elements = std::vector(_rows, std::vector<double>(_cols, 0));
}

std::vector<double> ScalarField::operator[](int index) const {
    return elements[index];
}

std::vector<double>& ScalarField::operator[](int index) {
    return elements[index];
}

ScalarField ScalarField::operator+(ScalarField b) {
    ScalarField result(rows, cols);

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            result[i][j] = elements[i][j] + b[i][j];
        }
    }

    return result;
}

ScalarField ScalarField::operator-(ScalarField b) {
    ScalarField result(rows, cols);

    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            result[i][j] = elements[i][j] - b[i][j];
        }
    }

    return result;
}

void ScalarField::set_bounds(double bnd_val) {
    if (rows == 0 || cols == 0) return;

    for (int i = 0; i < rows; i++) {
        elements[i][0] = elements[i][cols - 1] = bnd_val;
    }

    for (int j = 0; j < cols; j++) {
        elements[0][j] = elements[rows - 1][j] = bnd_val;
    }
}

void ScalarField::reset() {
    for (auto& vec : elements) std::fill(vec.begin(), vec.end(), 0);
}

Fluid::Fluid(int _nx, double _diff, double _visc, double _dt) {
    nx = _nx;
    diffusion = _diff;
    viscosity = _visc;
    dt = _dt;

    density_temp = ScalarField(nx, nx);
    density = ScalarField(nx, nx);

    velocity = VectorField (
        ScalarField(nx, nx),
        ScalarField(nx, nx)
    );

    velocity_prev = VectorField (
        ScalarField(nx, nx),
        ScalarField(nx, nx)
    );
}

VectorField::VectorField(ScalarField _x, ScalarField _y) {
    x = _x;
    y = _y;
}

VectorField VectorField::operator+(VectorField b) {
    return VectorField(x + b.x, y + b.y);
}

VectorField VectorField::operator-(VectorField b) {
    return VectorField(x - b.x, y - b.y);
}

// https://en.wikipedia.org/wiki/Gradient
VectorField VectorField::grad(double h) {
    VectorField result(x, y);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.rows - 1; j++) {
            result.x[i][j] = 0.5 * (x[i+1][j] - x[i-1][j]) * h;
            result.y[i][j] = 0.5 * (x[i][j+1] - x[i][j-1]) * h;
        }
    }

    return result;
}

// https://en.wikipedia.org/wiki/Divergence
ScalarField VectorField::div(double h) {
    ScalarField result(x.rows, x.cols);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.cols - 1; j++) {
            result[i][j] = -0.5 * (
                    x[i+1][j] - x[i-1][j] +
                    y[i][j+1] - y[i][j-1]

                ) / h;
        }
    }

    return result;
}

// https://scicomp.stackexchange.com/questions/21915/discrete-definitions-of-curl-nabla-times-f
ScalarField VectorField::curl(double h) {
    ScalarField result(x.rows, x.cols);

    for (int i = 1; i < x.rows - 1; i++) {
        for (int j = 1; j < y.rows - 1; j++) {
            double dv_dx = (y[i+1][j] - y[i-1][j]) / (2*h);
            double du_dy = (x[i][j+1] - x[i][j-1]) / (2*h);

            result[i][j] = dv_dx - du_dy;
        }
    }

    return result;
}

// Use the Gauss-Seidel method to solve a linear system of equations, whose coefficients are in the scalar fields x and x0.
void Fluid::gauss_seidel_solve(ScalarField &x, ScalarField &x0, double field_coefficient, int iter) {
    for (int k = 0; k < iter; k++) {
        for (int i = 1; i < x.rows - 1; i++) {
            for (int j = 1; j < x.cols - 1; j++) {
                x[i][j] =  (x0[i][j] + field_coefficient *
                               (x[i+1][j] + x[i-1][j] +
                                x[i][j+1] + x[i][j-1] )) * 1.0/(4*field_coefficient + 1);
            }
        }

        x.set_bounds();
    }
}

void Fluid::advect(ScalarField &field, ScalarField &field_prev, ScalarField &vel_x, ScalarField &vel_y) {
    int cell_x, cell_y;
    double dtx = dt * (nx - 2);
    double x, y;

    for (int i = 1; i < nx - 1; i++) { 
        for (int j = 1; j < nx - 1; j++) {
            x = constrain(i - dtx * vel_x[i][j], 0.5, nx + 0.5);
            cell_x = constrain(floorf(x), 0, nx-2);

            y = constrain(j - dtx * vel_y[i][j], 0.5, nx + 0.5);
            cell_y = constrain(floorf(y), 0, nx-2); 
                
            field[i][j] = bilinear_interpolation(field_prev, x, y, cell_x, cell_y);
        }
    }

    field.set_bounds();
}

void Fluid::project_velocity() {
    velocity.x.reset();
    velocity.y = velocity_prev.div(nx);

    velocity.y.set_bounds();

    gauss_seidel_solve(velocity.x, velocity.y, 1);
    velocity_prev = velocity_prev - velocity.grad(nx);

    velocity_prev.x.set_bounds();
    velocity_prev.y.set_bounds();
}

void Fluid::diffuse_velocity() {
    double velocity_diffusion = dt * viscosity * nx * nx;

    gauss_seidel_solve(velocity_prev.x, velocity.x, velocity_diffusion);
    gauss_seidel_solve(velocity_prev.y, velocity.y, velocity_diffusion);
}

void Fluid::diffuse_density() {
    double density_diffusion = dt * diffusion * nx * nx;
    gauss_seidel_solve(density_temp, density, density_diffusion);
}

void Fluid::step() {
    diffuse_velocity();
    project_velocity();

    advect(velocity.x, velocity_prev.x, velocity_prev.x, velocity_prev.y);
    advect(velocity.y, velocity_prev.y, velocity_prev.x, velocity_prev.y);

    diffuse_density();

    advect(density, density_temp, velocity.x, velocity.y);
}

void Fluid::add_density(int x, int y, double amt) {
    density[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt; 
}

void Fluid::add_velocity(int x, int y, double amt_x, double amt_y) {    
    velocity.x[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_x;
    velocity.y[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)] += amt_y;
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
            Vector2 cell_velocity = {.x = (float) velocity.x[i][j], .y = (float) velocity.y[i][j]};
            Vector2 cell_velocity_normal = Vector2Normalize(cell_velocity) * line_length;

            if (cell_velocity_normal == Vector2Zero()) DrawPixelV(position, colour);
            DrawLineEx(position, cell_velocity_normal + position, VELOCITY_LINE_THICKNESS, colour);
        }
    }
}

void Fluid::render_vorticity(int scale, Color colour, float particle_size) {
    ScalarField vorticity = velocity.curl();
    vorticity.set_bounds();

    for (int i = 2; i < nx - 2; i++) {
        for (int j = 2; j < nx - 2; j++) {
            double vort = vorticity[i][j];
            if (double_eq(vort, 0.0)) continue;

            double colour_rng = viscosity == 0 ? 0.02 : 0.05;
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
    velocity.x.reset();
    velocity.y.reset();
    velocity_prev.x.reset();
    velocity_prev.y.reset();
}