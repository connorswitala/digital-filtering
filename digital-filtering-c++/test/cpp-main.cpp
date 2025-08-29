#include "../df/df.hpp"

int main() {

    // Create configuration struct
    DFConfig config;

    // Configuration variables. 


    // Constructor
    DIGITAL_FILTER df(config);

    // Call filter procudure with timestep
    double dt = 1e-5;
    // df.filter(dt);
    df.get_rms();

    return 0;
}