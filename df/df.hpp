#pragma once

#include <iostream>
// #include <mpi.h>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono> 
#include "../pcg-cpp/include/pcg_random.hpp"

#define NOW chrono::high_resolution_clock::now();

constexpr double pi_c = -2 * M_PI; 

using namespace std;

typedef vector<double> Vector;

class DIGITAL_FILTER {

    private:

    int size, rank;         // Storage for size and rank of entire CFD simulation.
    int n_ranks;            // Storage for number of ranks you want to use for digital filtering.
    int Ny, Nz, n_cells;    // Storage for CFD domain without ghost cells.
    int Ny_max, Nz_max;     // Intermediate storage to find filter widths.
    double rand1, rand2;    // Intermediate storage for white-noise random number generator.
    double d_i, d_v;        // Inlet boundary layer height and viscous length scale. 

    Vector ru_ys, rv_ys, rw_ys;     // Data structures for random numbers in y_sweep filtering. 
    Vector ru_zs, rv_zs, rw_zs;     // Data structures for random numbers in z_sweep filtering.

    vector<int> Nu_ys, Nv_ys, Nw_ys;    // Data structures to hold filter half-widths for each CFD cell in y-sweep.
    vector<int> Nu_zs, Nv_zs, Nw_zs;    // Data structures to hold filter half-widths for each CFD cell in z-sweep.

    Vector bu_y, bv_y, bw_y;    // Data structures for filter coefficients in y-sweep.
    Vector bu_z, bv_z, bw_z;    // Data structures for filter coefficients in y-sweep.

    Vector u_fluc, v_fluc, w_fluc;              // Data structures for final filtered velocity fluctuations.
    Vector u_fluc_old, v_fluc_old, w_fluc_old;  // Data structures to hold old fluctuations.

    Vector R11, R21, R22, R33;

    int timestep;
    double dt;

    /**
     *  The following vector has a size of 6. N_holder[0] holds Nz_max for u', 
     *  N_holder[1] holds Ny_max for u'. The next 4 are for v' and w' Nz_max
     *  is even number indices [0, 2, 4] and Ny_max  is odd [1, 3, 5].
     */
    vector<int> N_holder; 

    Vector y, yc, z, dy, dz;

    double Iz_inn, Iz_out;
    Vector Iy, Iz;

    public:

    DIGITAL_FILTER(int Ny, int Nz, double d_i, double d_v, int n_ranks);

    void generate_white_noise();
    void calculate_filter_properties(); 
    void correlate_fields(); 
    void fake_grid();
    void filtering_sweeps(); 
    void filter(double dt_input); 
    void test(); 
    void display_data(Vector& v);
    void find_mean_variance(Vector& v);
    void get_RST(const string filename, int N_values, double rho_e, double U_e,  double mu_e, Vector rho_y);
    
    void read_grid(int N_x_vertices, int N_y_vertices);

    void write_tecplot(const string &filename);
    void write_tecplot_line(const string &filename);
};
