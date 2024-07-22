#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max);
double constrain(double value, double range_start, double range_end);

class Scalar_Field3D {
    public:
    
    int rows, cols, height;
    std::vector<std::vector<std::vector<double>>> elements;

    Scalar_Field3D() = default;
    Scalar_Field3D(int rows, int cols, int height);
    std::vector<std::vector<double>> operator[] (int index) const; 
    std::vector<std::vector<double>>& operator[] (int index); 
};

struct Vector_Field3D {
    Scalar_Field3D x, y, z;
};

class Fluid3D {
    public:

    int nx;
    double dt;
    double diff;
    double visc;
    
    Scalar_Field3D s;
    Scalar_Field3D density;
    
    Vector_Field3D v;
    Vector_Field3D v0;

    Fluid3D() = default;
    Fluid3D(int _nx, double _diff, double _visc, double _dt);
 
    static Scalar_Field3D curl(Vector_Field3D field, double h = 1.0);
    static Scalar_Field3D div(Vector_Field3D field, double h = 1.0);

    void add_density(int x, int y, int z, double amt); 
    void add_velocity(int x, int y, int z, double amt_x, double amt_y, double amt_z);

    void set_bnd(Scalar_Field3D &x, double bnd_val = 0.0);
    void gauss_seidel(Scalar_Field3D &x, Scalar_Field3D &x0, double a, int iter = 20);
    void advect(Scalar_Field3D &d, Scalar_Field3D &d0, Scalar_Field3D &vx, Scalar_Field3D &vy, Scalar_Field3D &vz);
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
