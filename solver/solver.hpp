#include <cstdio>

#define IX(x, y) (x + y*nx)

class FluidCube {
    public:

    int nx;
    int iter;
    double dt;
    double diff;
    double visc;
    
    double *s;
    double *density;
    
    double *Vx;
    double *Vy;

    double *Vx0;
    double *Vy0;

    FluidCube(int _nx, double _diff, double _visc, double _dt); 

    void set_bnd(int b, double *x);
    void lin_solve(int b, double *x, double *x0, double a, double c);
    void diffuse (int b, double *x, double *x0);
    void advect(int b, double *d, double *d0, double *velocX, double *velocY);
    void project(double *velocX, double *velocY, double *p, double *div);
    void step();
    void add_density(int x, int y, double amount); 
    void add_velocity(int x, int y, double amountX, double amountY);
    void render_density(int scale);
    void print_density(FILE *out = stdout);
};
