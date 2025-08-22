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

#define NOW chrono::high_resolution_clock::now(); // This is a macro to get the current time.
constexpr double pi_c = -2 * 3.14159265358979323846; // This is a constant used in finding the filter coefficients.

using namespace std;
typedef vector<double> Vector;

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
struct df_config {
    double d_i, rho_e, U_e, mu_e;
    int vel_file_offset, vel_file_N_values;
    string grid_file, vel_fluc_file;
};

class DIGITAL_FILTER {

    private:

    string grid_file;           // File containing the grid data.
    string vel_fluc_file;       // File containing the velocity fluctuations.
    int vel_file_offset;        // Number of lines to skip in fluctuation header file
    int vel_file_N_values;      // Number of lines to read in fluctuation file

    int Ny, Nz, n_cells;        // Storage for CFD domain without ghost cells.
    double rand1, rand2;        // Intermediate storage for white-noise random number generator.
    double d_i, d_v;            // Inlet boundary layer height and viscous length scale. 
    double rho_e, U_e, mu_e;    // Freestream flow parameters
    Vector rhoy, Uy, My, Ty;    // Mean distribution of flow variable in wall normal direction.
    double mu_wall;             // Wall viscosity

    Vector rho_fluc, T_fluc;
    Vector R11, R21, R22, R33;                  // Reynolds stress terms

    int rms_counter;
    double dt;

    /**
     *  The following vector has a size of 6. N_holder[0] holds Nz_max for u', 
     *  N_holder[1] holds Ny_max for u'. The next 4 are for v' and w' Nz_max
     *  is even number indices [0, 2, 4] and Ny_max  is odd [1, 3, 5].
     */
    vector<int> N_holder; 

    Vector y, yc, z, dy, dz;    // Geometry data

    // Coefficient of integral length scales
    double Iz_inn, Iz_out; 
    Vector Iy, Iz;

    public:

    FilterField u, v, w; 

    DIGITAL_FILTER(df_config config);

    void allocate_data_structures(FilterField& F);          // Allocates data structures for the filter.

    void generate_white_noise();                // Generates white noise for the filtering.
    void calculate_filter_properties(FilterField& F);         // Calculates filter properties for the filtering.
    void correlate_fields(FilterField& F);                    // Correlates new fluctuations from old ones and then scales using RST.
    void apply_RST_scaling();                // Used for t = 0 and only scales with RST.
    void filtering_sweeps(FilterField& F);                    // Filters the velocity fluctuations in y and z sweeps.
    void filter(double dt_input);               // Runs all the filtering processes and updates the fluctuations.
    void get_rho_T_fluc();                      // Calculates the fluctuations for temperature and density.
    void set_old();                             // Sets the old fluctuations to the new ones. 
    void test();                                // Used for testing purposes.        
    void display_data(Vector& v);               // Displays the data in the vector.
    void find_mean_variance(Vector& v);         // Finds the mean and variance of the vector.
    


    void allocate_rms_structures(FilterField F); // Allocates data structures for the RMS values.
    void rms_add();
    void get_rms();
    void plot_rms(); 

    /** 
     *  The following function opens a file that contains rms values for u', w', and w' that are normalized by the Morokovin
     *  velocity scale u* = u_tau * sqrt(rho_w / rho(y)). It then multiplies these values by a turbulent boundary layer estimate
     *  for wall shear stress depending on an estimate turbulent skin friction coefficient and for Re_x using the inlet boundary 
     *  layer height, d_i. Once the fluctuations are calculated, it find the reynolds stresses R11, R21, R22, and R33 for each 
     *  cell to be used when scaling the fluctuations.
     */
    void get_RST();   

    void rho_y_test(); // This function creates a test distribution of density in the wall-normal direction to be used with Morkovin scaling.
             
    void read_grid(); // This function reads the grid from a file and populates the y, yc, z, dy, dz vectors.

    // The following functions write the data to a file in Tecplot format to visualze the distubution of the fluctuations.
    void write_tecplot(const string &filename);
    void write_tecplot_line(const string &filename);
};
