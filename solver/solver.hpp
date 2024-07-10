#include <cmath>
#include <cstdio>
#include <vector>

#define Matrix Raylib_Matrix
#include "../raylib/include/raylib.h"
#undef Matrix

class Matrix {
    public:
    
    int rows, cols;
    std::vector<std::vector<double>> elements;

    Matrix() = default;
    Matrix(int rows, int cols);
    std::vector<double> operator[] (int index) const; 
    std::vector<double>& operator[] (int index); 
};

struct Matrix_Vector {
    Matrix x, y;
};

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max);
double constrain(double value, double range_start, double range_end);

class Fluid {
    public:

    int nx;
    int iter;
    double dt;
    double diff;
    double visc;
    
    Matrix s;
    Matrix density;
    
    Matrix_Vector v;
    Matrix_Vector v0;

    Fluid() = default;
    Fluid(int _nx, double _diff, double _visc, double _dt);

    void set_bnd(Matrix &x);
    void gauss_seidel(Matrix &x, Matrix &x0, double a);
    void advect(Matrix &d, Matrix &d0, Matrix &velocX, Matrix &velocY);
    void project_velocity();
    void diffuse_velocity();
    void diffuse_density();
    void step();
    void add_density(int x, int y, double amount); 
    void add_velocity(int x, int y, double amountX, double amountY);
    void render_density(int scale, Color colour, float particle_size);
    void render_velocity(int scale, Color colour, float line_length);
    void fade_density(double fade);
    void print_density(FILE *out = stdout);
    void reset();
};
