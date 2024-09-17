#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max);
double constrain(double value, double range_start, double range_end);

class Scalar_Field {
    public:
    int rows, cols;
    std::vector<std::vector<double>> elements;

    Scalar_Field() = default;
    Scalar_Field(int rows, int cols);
    std::vector<double> operator[] (int index) const; 
    std::vector<double>& operator[] (int index); 
};

class Vector_Field {
    public:
    Scalar_Field x, y;

    Scalar_Field curl(double h = 1.0);
    Scalar_Field div(double h = 1.0);
};

class Fluid {
    public:

    int nx;
    double dt;
    double diff;
    double visc;
    
    Scalar_Field s;
    Scalar_Field density;
    
    Vector_Field v;
    Vector_Field v0;

    Fluid() = default;
    Fluid(int _nx, double _diff, double _visc, double _dt);

    void add_density(int x, int y, double amt); 
    void add_velocity(int x, int y, double amt_x, double amt_y);

    void set_bnd(Scalar_Field &x, double bnd_val = 0.0);
    void gauss_seidel(Scalar_Field &x, Scalar_Field &x0, double a, int iter = 20);
    void advect(Scalar_Field &d, Scalar_Field &d0, Scalar_Field &vx, Scalar_Field &vy);
    void project_velocity();
    void diffuse_velocity();
    void diffuse_density();
    void step(); 
    void render_density(int scale, Color colour, float particle_size);
    void render_velocity(int scale, Color colour, float line_length);
    void render_vorticity(int scale, Color colour, float particle_size);
    void fade_density(double fade);
    void print_density(FILE *out = stdout);
    void reset();
};
