#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<random>

using namespace std;

typedef vector<double> Vector;

class DIGITAL_FILTER {

    private:

    int Ny, Nz, n_cells;
    int Nx_max, Ny_max, Nz_max;
    double rand1, rand2;
    double d_i, d_v;  

    Vector r_xs, r_ys, r_zs; 
    Vector b_u, b_v, b_w, b_tilde;  
    Vector y, dy;
    double dx, dz;
    double Iz_inn, Iz_out;
    Vector Ix, Iy, Iz, N_ys, N_xs, N_zs;

    public:

    DIGITAL_FILTER(int Ny, int Nz, double d_i);

    double Box_Muller_rand(mt19937& gen, uniform_real_distribution<double>& dist);
    void white_noise();
    void calculate_Ns();
    void calculate_bs();
    void fake_grid();
    void test(); 
    void find_mean_variance(Vector& v);
};
