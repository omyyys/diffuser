#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

typedef std::vector<std::vector<double>> Grid;

struct Grid_Vector {
    Grid x, y;
};

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max);

class Fluid {
    public:

    int nx;
    int iter;
    double dt;
    double diff;
    double visc;
    
    Grid s;
    Grid density;
    
    Grid_Vector v;
    Grid_Vector v0;

    Fluid() = default;
    Fluid(int _nx, double _diff, double _visc, double _dt);

    void set_bnd(Grid &x);
    void gauss_seidel(Grid &x, Grid &x0, double a);
    void advect(Grid &d, Grid &d0, Grid &velocX, Grid &velocY);
    void project_velocity();
    void diffuse_velocity();
    void diffuse_density();
    void step();
    void add_density(int x, int y, double amount); 
    void add_velocity(int x, int y, double amountX, double amountY);
    void render_density(int scale, Color colour, float particle_size);
    void print_density(FILE *out = stdout);
    void reset();
};
