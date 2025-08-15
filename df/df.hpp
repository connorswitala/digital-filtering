#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include <string>
#include <fstream>

using namespace std;

typedef vector<double> Vector;

class DIGITAL_FILTER {

    private:

    int Ny, Nz, n_cells;
    int Ny_max, Nz_max;
    double rand1, rand2;
    double d_i, d_v;  

    Vector ru_ys, ru_zs, rv_ys, rv_zs, rw_ys, rw_zs; 
    vector<int> Nu_ys, Nu_zs, Nv_ys, Nv_zs, Nw_ys, Nw_zs;
    Vector bu_y, bu_z, bv_y, bv_z, bw_y, bw_z, b_tilde;  
    Vector holder, u_fluc, v_fluc, w_fluc;

    vector<int> N_holder; // size of 6. First two are Nz then Ny max of u'. Then v' and w'.

    Vector y, yc, z, dy, dz;
    double Iz_inn, Iz_out;
    Vector Iy, Iz;

    public:

    DIGITAL_FILTER(int Ny, int Nz, double d_i);

    void generate_white_noise();
    void calculate_filter_properties(); 
    void apply_RST(); 
    void fake_grid();
    void filter(); 
    void test(); 
    void display_data(Vector& v);
    void find_mean_variance(Vector& v);

    void write_tecplot(const string &filename);
};
