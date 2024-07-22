#include <cmath>

#include "solver-3d.hpp"

inline bool double_eq(double a, double b, double eps = 1e-6) {
    return fabs(a - b) < eps;
}

double constrain(double value, double range_start, double range_end) {
    value = (value > range_end) ? range_end : value;
    value = (value < range_start) ? range_start : value;

    return value;
}

double lin_interpolation(double value, double cur_min, double cur_max, double res_min, double res_max) {
    value = constrain(value, cur_min, cur_max); 

    double cur_range = cur_max - cur_min;
    double res_range = res_max - res_min;

    double normalised_value = (value - cur_min) / cur_range;

    return res_min + normalised_value * res_range;
}


Scalar_Field3D::Scalar_Field3D(int _rows, int _cols, int _height) {
    rows = _rows;
    cols = _cols;
    height = _height;
    
    elements = std::vector<std::vector<std::vector<double>>>(rows, std::vector<std::vector<double>>(cols, std::vector<double>(height)));
}

std::vector<std::vector<double>> Scalar_Field3D::operator[] (int index) const {
    return elements[index];
}

std::vector<std::vector<double>>& Scalar_Field3D::operator[] (int index) {
    return elements[index];
}

Fluid3D::Fluid3D(int _nx, double _diff, double _visc, double _dt) {
    nx = _nx;
    diff = _diff;
    visc = _visc;
    dt = _dt;

    s = Scalar_Field3D(nx, nx, nx);
    density = Scalar_Field3D(nx, nx, nx);

    v = {
        .x = Scalar_Field3D(nx, nx, nx),
        .y = Scalar_Field3D(nx, nx, nx)
    };

    v0 = {
        .x = Scalar_Field3D(nx, nx, nx),
        .y = Scalar_Field3D(nx, nx, nx)
    };
}

Scalar_Field3D Fluid3D::div(Vector_Field3D field, double h) {
    Scalar_Field3D result(field.x.rows, field.x.cols, field.x.height);

    for (int i = 1; i < field.x.rows - 1; i++) {
        for (int j = 1; j < field.y.cols - 1; j++) {
            for (int k = 1; k < field.z.height - 1; k++) {
                result[i][j][k] = -0.5 * (
                        field.x[i+1][j][k] -
                        field.x[i-1][j][k] +
                        field.y[i][j+1][k] -
                        field.y[i][j-1][k] +
                        field.z[i][j][k+1] -
                        field.z[i][j][k-1]
                    ) / h;
            }
        }
    }

    return result;
}

void Fluid3D::set_bnd(Scalar_Field3D &x, double bnd_val) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            for (int k = 0; k < nx; k++) {
                if (i == 0 || i == nx - 1 || j == 0 || j == nx - 1 || k == 0 || k == nx - 1) {
                    x[i][j][k] = 0;
                }
            }
        }
    } 
}

void Fluid3D::gauss_seidel(Scalar_Field3D &x, Scalar_Field3D &x0, double a, int iter) {
    double c = 4*a + 1; 

    for (int _ = 0; _ < iter; _++) {
        for (int i = 1; i < x.rows - 1; i++) {
            for (int j = 1; j < x.cols - 1; j++) {
                for (int k = 1; k < x.height - 1; k++) {
                    x[i][j][k] = (x0[i][j][k] + a *
                                    (x[i+1][j][k] +
                                     x[i-1][j][k] +
                                     x[i][j+1][k] +
                                     x[i][j-1][k] +
                                     x[i][j][k+1] +
                                     x[i][j][k-1] )) * (1.0/c);
                }
            }
        }
        set_bnd(x);
    }
}

void Fluid3D::advect(Scalar_Field3D &d, Scalar_Field3D &d0, Scalar_Field3D &v_x, Scalar_Field3D &v_y, Scalar_Field3D &v_z) {
    int ix, iy, iz;

    double dtx = dt * (nx - 2);
    
    double s0, s1, t0, t1, u0, u1;
    double x, y, z;

    for (double i = 1; i < nx - 1; i++) { 
        for (double j = 1; j < nx - 1; j++) {
            for (double k = 1; k < nx - 1; k++) {
                x = constrain(i - dtx * v_x[i][j][k], 0.5, nx + 0.5);
                ix = constrain(floorf(x), 0, nx-2);

                y = constrain(j - dtx * v_y[i][j][k], 0.5, nx + 0.5);
                iy = constrain(floorf(y), 0, nx-2);

                z = constrain(j - dtx * v_z[i][j][k], 0.5, nx + 0.5);
                iz = constrain(floorf(z), 0, nx-2);
                
                s1 = x - ix; 
                s0 = 1.0 - s1; 
                t1 = y - iy; 
                t0 = 1.0 - t1;
                u1 = z - iz;
                u0 = 1.0f - u1;
                
                d[i][j][k] = s0 * ( t0 * (u0 * d0[ix][iy][iz]
                                +u1 * d0[ix][iy][iz])
                                +( t1 * (u0 * d0[ix][iy][iz]
                                +u1 * d0[ix][iy][iz])))
                                +s1 * ( t0 * (u0 * d0[ix][iy][iz]
                                +u1 * d0[ix][iy][iz])
                                +( t1 * (u0 * d0[ix][iy][iz]
                                +u1 * d0[ix][iy][iz])));
            }
        }
    }

    set_bnd(d);
}

void Fluid3D::project_velocity() {
    for (auto& vec_outer : v.x.elements) {
        for (auto& vec : vec_outer) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    v.y = div(v0, nx);

    set_bnd(v.y);
    gauss_seidel(v.x, v.y, 1);
    
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            for (int k = 1; k < nx - 1; k++) {
                v0.x[i][j][k] -= 0.5 * (v.x[i+1][j][k] - v.x[i-1][j][k]) * nx;
                v0.y[i][j][k] -= 0.5 * (v.x[i][j+1][k] - v.x[i][j-1][k]) * nx;
                v0.z[i][j][k] -= 0.5 * (v.x[i][j+1][k+1] - v.x[i][j-1][k-1]) * nx;
            }
        }
    }

    set_bnd(v0.x);
    set_bnd(v0.y);
    set_bnd(v0.z);
}

void Fluid3D::diffuse_velocity() {
    double a = dt * visc * nx * nx;

    gauss_seidel(v0.x, v.x, a);
    gauss_seidel(v0.y, v.y, a);
}

void Fluid3D::diffuse_density() {
    double a = dt * diff * nx * nx;

    gauss_seidel(s, density, a);
}

void Fluid3D::step() {
    diffuse_velocity();

    project_velocity();

    advect(v.x, v0.x, v0.x, v0.y, v0.z);
    advect(v.y, v0.y, v0.x, v0.y, v0.z);

    diffuse_density();
    advect(density, s, v.x, v.y, v.z);
}

void Fluid3D::add_density(int x, int y, int z, double amt) {
    density[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)][constrain(z, 0, nx-1)] += amt; 
}

void Fluid3D::add_velocity(int x, int y, int z, double amt_x, double amt_y, double amt_z) {    
    v.x[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)][constrain(z, 0, nx-1)] += amt_x;
    v.y[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)][constrain(z, 0, nx-1)] += amt_y;
    v.z[constrain(x, 0, nx-1)][constrain(y, 0, nx-1)][constrain(z, 0, nx-1)] += amt_z;
}

void Fluid3D::fade_density(double fade) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            for (int k = 0; k < nx; k++) {
                density[i][j][k] = constrain(density[i][j][k] - fade, 0, 1);
            }
        }
    }
}

void Fluid3D::reset() {
    for (auto& outer_vec : s.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }
   
    for (auto& outer_vec : density.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v.x.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v.y.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v.z.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v0.x.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v0.y.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }

    for (auto& outer_vec : v0.z.elements) {
        for (auto& vec : outer_vec) {
            std::fill(vec.begin(), vec.end(), 0);
        }
    }
}
