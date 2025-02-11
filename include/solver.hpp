#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

double translate_value(double value, double cur_min, double cur_max, double res_min, double res_max);
double constrain(double value, double range_start, double range_end);

class ScalarField {
    public:
    int rows, cols;
    std::vector<std::vector<double>> elements;

    ScalarField() = default;
    ScalarField(int rows, int cols);
    std::vector<double> operator[] (int index) const; 
    std::vector<double>& operator[] (int index); 
    ScalarField operator+(ScalarField b); 
    ScalarField operator-(ScalarField b);

    void set_bounds(double bnd_val = 0.0);
    void reset();
};

class VectorField {
    public:
    ScalarField x, y;

    VectorField() = default;
    VectorField(ScalarField _x, ScalarField _y);
    VectorField operator+(VectorField b); 
    VectorField operator-(VectorField b);

    VectorField grad(double h = 1.0);
    ScalarField curl(double h = 1.0);
    ScalarField div(double h = 1.0);
};

class Fluid {
    public:

    int nx;
    double dt;
    double diffusion;
    double viscosity;
    
    ScalarField density_temp;
    ScalarField density;
    
    VectorField velocity;
    VectorField velocity_prev;

    Fluid() = default;
    Fluid(int _nx, double _diff, double _visc, double _dt);

    void add_density(int x, int y, double amt); 
    void add_velocity(int x, int y, double amt_x, double amt_y);

    void gauss_seidel_solve(ScalarField &x, ScalarField &x0, double field_coefficient, int iter = 20);
    void advect(ScalarField &field, ScalarField &field_prev, ScalarField &vel_x, ScalarField &vel_y);
    void project_velocity();
    void diffuse_velocity();
    void diffuse_density();
    void step(); 
    void render_density(int scale, Color colour, float particle_size);
    void render_velocity(int scale, Color colour, float line_length);
    void render_vorticity(int scale, Color colour, float particle_size);
    void fade_density(double fade);
    void reset();
};
