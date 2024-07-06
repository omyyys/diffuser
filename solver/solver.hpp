#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

typedef std::vector<double> Grid;

struct Grid_Vector {
    Grid x, y;
};

#define IX(x, y) (x + y*nx)

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

    Fluid(int _nx, double _diff, double _visc, double _dt);

    void set_bnd(int b, Grid &x);
    void gauss_geidel(int b, Grid &x, Grid &x0, double a, double c);
    void diffuse (int b, double diff_or_visc, Grid &x, Grid &x0);
    void advect(int b, Grid &d, Grid &d0, Grid &velocX, Grid &velocY);
    void project(Grid &velocX, Grid &velocY, Grid &p, Grid &div);
    void diffuse_velocity();
    void diffuse_density();
    void velocity_step();
    void density_step();
    void step();
    void add_density(int x, int y, double amount); 
    void add_velocity(int x, int y, double amountX, double amountY);
    void render_density(int scale, Color colour, float particle_size);
    void print_density(FILE *out = stdout);
    void reset();
};
