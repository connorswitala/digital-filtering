#include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(int Ny, int Nz, double d_i) : Ny(Ny), Nz(Nz), d_i(d_i) {


    d_v = 0.002 * d_i; 
    dx = 0.0001;
    dz = 0.0001;

    n_cells = Ny * Nz;  

    b_tilde = Vector(n_cells);
    b_u = Vector(n_cells);
    b_v = Vector(n_cells);
    b_w = Vector(n_cells);

    Ix = Vector(n_cells); 
    Iy = Vector(n_cells);
    Iz = Vector(n_cells);
    N_ys = Vector(n_cells);
    N_xs = Vector(n_cells);
    N_zs = Vector(n_cells);
}

double DIGITAL_FILTER::Box_Muller_rand(mt19937& gen, uniform_real_distribution<double>& dist) {

    double u1 = dist(gen);
    double u2 = dist(gen);
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1) * cos(2 * M_PI * u2));
}

void DIGITAL_FILTER::white_noise() {

    static mt19937 gen(random_device{}());
    static normal_distribution<> dist(0.0, 1.0);

    for (int i = 0; i < n_cells; i++) {

        // Generating white noise for u'
        u_rand[i] = dist(gen);
        v_rand[i] = dist(gen);
        w_rand[i] = dist(gen);
    }
}

void DIGITAL_FILTER::calculate_Ns() {

    // Find N for filtering in y-direction and create size of r (white noise)

    Iz_out = 0.4 * d_i;
    Iz_inn = 150 * d_v;
    Nx_max = 0; Ny_max = 0; Nz_max = 0; 

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {

            double y_m = ( y[(j + 1) * (Nz + 1) + i] - y[j * (Nz + 1) + i] ) / 2;  
            Iz[j * Nz + i] = Iz_inn + (Iz_out - Iz_inn) * ( (1 + tanh(y_m / d_i - 0.2) / 0.03) / 2 ); 

            double n_int = Iz[j * Nz + i] / dz;
            N_zs[j * Nz + i] = static_cast<int>(n_int);           // Need to be integer or a double

            if (n_int > Nz_max) Nz_max = static_cast<int>(n_int); 
        }
    }

    r_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);  

    // Find N for filtering in z-direction and create size of r

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {
            Iy[j * Nz + i] = 0.67 * Iz[j * Nz + i];

            double n_int = Iy[j * Nz + i] / dy[j * Nz + i]; 
            N_zs[j * Nz + i] = static_cast<int>(n_int);

            if (n_int > Ny_max) Ny_max = static_cast<int>(n_int);
        }
    }

    r_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);


    // Find N for filtering in x-direction and create size of r

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {
            Ix[j * Nz + i] = 
        }
    }




}


void DIGITAL_FILTER::calculate_bs() {

    // Find N for filtering in y-direction

    Iy_out = 0.4 * d_i;
    Iy_inn = 150 * d_v;

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {

            double y_m = ( y[(j + 1) * (Nz + 1) + i] - y[j * (Nz + 1) + i] ) / 2;  
            Iy[j * Nz + i] = Iy_inn + (Iy_out - Iy_inn) * ( (1 + tanh(y_m / d_i - 0.2) / 0.03) / 2 ); 
            Iz[j * Nz + i] = 0.67 * Iy[j * Nz + i];  
            Ix[j * Nz + i] = 0.8 * d_i; 

            double n_int = Iy[j * Nz + i] / dy[j * Nz + i];
            ny[i] = static_cast<int>(n_int);         
        }
    }


    for (int i = 0; i < n_cells; ++i) {
        b_tilde[i]  exp()
    }



}




void DIGITAL_FILTER::test() {
    white_noise();
    find_mean_variance(u_rand); 
}
void DIGITAL_FILTER::fake_grid() {
    double y_max = 0.01;
    double eta = 0.0;  
    double a = 2.0; 

    y = Vector((Ny + 1) * (Nz + 1), 0.0); 
    dy = Vector(Ny * Nz, 0.0); 

    for (int i = 0; i < Nz + 1; ++i) {
        for (int j = 0; j < Ny + 1; ++j) {
            eta = ( (j + 1) * y_max / (Ny + 1)) / y_max;   
            y[j * (Nz + 1) + i] = y_max * (1 - tanh(a * eta) / tanh(a)); 
        }
    }

    for (int i = 0; i < Nz; ++i) {
        for (int j = 0; j < Ny; ++j) {
            dy[j * Nz + i] = y[(j + 1) * (Nz + 1) + i] - y[j * (Nz + 1) + i];
        }
    }
}
void DIGITAL_FILTER::find_mean_variance(Vector& v) {
    double mean_sum = 0.0, var_sum = 0.0;

    for (int i = 0; i < n_cells; ++i) {
        mean_sum += v[i];
    }

    double mean = mean_sum / n_cells;

    for (int i = 0; i < n_cells; ++i) {
        var_sum += (v[i] - mean) * (v[i] - mean);
    }

    double variance = var_sum / n_cells;

    cout << "Mean: " << mean << "\t Variance: " << variance << endl;
}
