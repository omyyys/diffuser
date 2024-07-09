#include <cmath>
#include <cstdio>
#include <vector>

#define Matrix Raylib_Matrix
#include "../raylib/include/raylib.h"
#undef Matrix

typedef std::vector<std::vector<double>> Matrix;

struct Matrix_Vector {
    Matrix x, y;
};

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max);

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
    void fade_density(double fade);
    void print_density(FILE *out = stdout);
    void reset();
};
