#include "../df/df.hpp"

using namespace std;

int main() {

    string file1 = "hi";
    string file2 = "../files/M6Tw025_Stat.dat";
    double d_i = 0.0013;
    double rho_e = 0.01, U_e = 800, mu_e = 1.8e-5;
    int offset = 142; 

    DIGITAL_FILTER df(d_i, rho_e, U_e, mu_e, file1, file2, offset, 330);
    df.test();

    return 0;
}