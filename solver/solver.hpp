#include <cmath>
#include <cstdio>
#include <vector>

#include "../raylib/include/raylib.h"

typedef std::vector<double> Grid;

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
    
    Grid Vx;
    Grid Vy;

    Grid Vx0;
    Grid Vy0;

    Fluid(int _nx, double _diff, double _visc, double _dt);

    void set_bnd(int b, Grid &x);
    void lin_solve(int b, Grid &x, Grid &x0, double a, double c);
    void diffuse (int b, double diff_or_visc, Grid &x, Grid &x0);
    void advect(int b, Grid &d, Grid &d0, Grid &velocX, Grid &velocY);
    void project(Grid &velocX, Grid &velocY, Grid &p, Grid &div);
    void step();
    void add_density(int x, int y, double amount); 
    void add_velocity(int x, int y, double amountX, double amountY);
    void render_density(int scale, Color colour, float particle_size);
    void print_density(FILE *out = stdout);
    void reset();
};
