#include "../df/df.hpp"

int main() {

    // Create configuration struct
    df_config config;

    // Configuration variables. 
    config.grid_file = "grid_file_here";
    config.vel_fluc_file = "../files/M6Tw025_Stat.dat";
    config.rho_e = 0.044;
    config.U_e = 869.1;
    config.mu_e = 1.8e-5;
    config.d_i = 0.0013;
    config.vel_file_offset = 142;
    config.vel_file_N_values = 330;
  

    // Constructor
    DIGITAL_FILTER df(config);

    // Call filter procudure with timestep

    double dt = 1e-8;
    df.filter(dt);

    return 0;
}