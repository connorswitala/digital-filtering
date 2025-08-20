#include "../df/df.hpp"

int main() {

    string grid_file = "grid_file_here";
    string vel_fluc_file = "../files/M6Tw025_Stat.dat";
    double U_e = 869.1, rho_e = 0.044, mu_e = 1.8e-5;
    double d_i = 0.0013;

    DIGITAL_FILTER df(d_i, rho_e, U_e, mu_e, grid_file, vel_fluc_file, 142, 330);
    df.test();

    return 0;
}