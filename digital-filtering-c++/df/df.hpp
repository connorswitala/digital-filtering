#pragma once

#include <iostream>
// #include <mpi.h>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <chrono> 
#include "../pcg-cpp/include/pcg_random.hpp"

#define NOW chrono::high_resolution_clock::now(); // This is a macro to get the current time.
constexpr double pi_c = -2.0 * 3.14159265358979323846; // This is a constant used in finding the filter coefficients.

using namespace std;
typedef vector<double> Vector;

//  FilterField is a struct that contains all necessary 
//  data structures for filring the velocity fluctuations.
//  There is one for each velocity component, u, v, and w.
struct FilterField {
    Vector  by, bz, 
            r_ys, r_zs,
            rms_added, rms, 
            filt_old, filt, fluc;

    vector<int> N_ys, N_zs,
                by_offsets, bz_offsets;
    double Iz_inn, Iz_out, Lt;
    int Nz_max, Ny_max; 
};

// DFConfig is a structure that contains all the inputs necessary
// to use digital filtering.
struct DFConfig {
    double  d_i,    // Inlet delta (boundary layer height)
            rho_e,  // Freestream density
            U_e,    // Freestream velocity 
            mu_e;   // Freestream viscosity

    int vel_file_offset,    // Number of header lines in fluctuation file
        vel_file_N_values;  // Number of fluctuation data points

    string  grid_file,      // Name of grid file (might be deprecated)
            vel_fluc_file;  // Name of fluctuation file
};

// Digital_Filter class
class DIGITAL_FILTER {

    private:

    int Ny, Nz, n_cells;        // Storage for CFD domain.
    double N_in;                // Duan file number of values to rrea din.

    Vector rho_fluc, T_fluc;    // Thermodynamic fluctuations
    Vector R11, R21, R22, R33;  // Reynolds stress terms
    Vector R11_in, R21_in, R22_in, R33_in;  // Reynolds stress terms

    int rms_counter;            // Counter to find mean of fluctuations
    double dt;                  // Input dt for filtering

    // Geometry vectors
    Vector y,       // y vertices for entire inflow size (Ny * Nz)
        yc,         // y cell-center for entire inflow
        dy,         // cell heights for entire inflow
        dz,         // cell widths for entire inflow
        yc_d,       // y cell centers normalized by d_i
        yin_d,      // y locations normalized by delta (from Duan data)
        y_in,       // y locations (from Duan data)
        ydline,     // Assuming y stretching is constant across z, this is all y locations normalized by our delta for DF
        yline;      // Same as above but not normalized;

    double d_i, d_v;                        // Inlet boundary layer height and viscous length scale. 
    double rho_e, U_e, T_e;                 // Freestream flow parameters
    double U_w, rho_w, T_w, P, mu, gcon;    // Wall values and global constants
    string line_file;                       // line.dat filename
    Vector Us, Ts, rhos, Ps, Ms;            // Mean wall-normal profiles
    double u_tau, tau_w;                    // Importans wall quantities

    Vector T_rms, rho_rms, T_rms_added, rho_rms_added;  // RMS data structures
    public:

    FilterField u, v, w; 

    DIGITAL_FILTER(DFConfig config);

    // =====: Filtering functions :=====
    void read_grid();                                   // Fills geometry data structures and find vector sizes
    void allocate_data_structures(FilterField& F);      // Allocates data structures for the filter.
    void calculate_filter_properties(FilterField& F);   // Calculates filter properties
    void get_RST_in();                                  // Reads in fluctuation file and calculates RST
    void generate_white_noise();                        // Generates white noise with mean = 0 and variance of 1
    void filtering_sweeps(FilterField& F);              // Filters the velocity fluctuations in y and z sweeps.
    void correlate_fields(FilterField& F);              // Correlates new fluctuations from old ones
    void apply_RST_scaling();                           // Scales fluctuations by RST
    void filter(double dt_input);                       // Runs all the filtering processes and updates the fluctuations.
    void get_rho_T_fluc();                              // Calculates the fluctuations for temperature and density.
    void read_line_file();                              // Reads in line.dat file to get wall normal profile.

    // ====: Debugging functions :=====
    //  (to be removed when finished)      
    void display_data(Vector& v);                       // Displays the data in the vector.

    // =====: RMS functions :=====
    //  (to be removed when finished)
    void allocate_rms_structures(FilterField& F);       // Allocates data structures for the RMS values.
    void rms_add();                                     // Adds square of fluctuations to rms sum every timestep
    void get_rms();                                     // Function designed to be called on its own
    void plot_rms();                                    // Plots RMS data for u', v', and w'

    // =====: Plotting functions : =====
    //  (to be removed when finished)
    void write_tecplot(const string &filename);         // Plots fluctuations
    void plot_RST_lerp();                               // Plots Reynolds stress terms
    void write_csv(const std::string &filename);

    Vector linear_interpolate(
    const vector<double>& y_data,
    const vector<double>& f_data,
    const vector<double>& y_new);
};


