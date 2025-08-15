#include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(int Ny, int Nz, double d_i) : Ny(Ny), Nz(Nz), d_i(d_i) {


    d_v = 0.002 * d_i; 

    n_cells = Ny * Nz;  

    Iy = Vector(n_cells);
    Iz = Vector(n_cells);
    Nu_ys = vector<int>(n_cells);
    Nu_zs = vector<int>(n_cells);
    Nv_ys = vector<int>(n_cells);
    Nv_zs = vector<int>(n_cells);
    Nw_ys = vector<int>(n_cells);
    Nw_zs = vector<int>(n_cells);
    N_holder = vector<int>(6); 
    holder = Vector(n_cells);
    u_fluc = Vector(n_cells);
    v_fluc = Vector(n_cells);
    w_fluc = Vector(n_cells);

    
}

void DIGITAL_FILTER::generate_white_noise() {

    static mt19937 gen(random_device{}());
    static normal_distribution<> dist(0.0, 1.0);

    for (auto &val : ru_ys) {
        val = dist(gen);
    }

    for (auto &val : ru_zs) {
        val = dist(gen);
    }

    for (auto &val : rv_ys) {
        val = dist(gen);
    }

    for (auto &val : rv_zs) {
        val = dist(gen);
    }

    for (auto &val : rw_ys) {
        val = dist(gen);
    }

    for (auto &val : rw_zs) {
        val = dist(gen);
    }
}

/**
 *  This function uses prescribed integral length scales Iz and Iy to find the filter width
 *  for each cell as well as calculating the convolution coefficients and create the storage for 
 *  white noise vectors.
 */
void DIGITAL_FILTER::calculate_filter_properties() {

    /**
    *       Find N for u' filtering in the z-direction.
    */

    Iz_out = 0.4 * d_i;
    Iz_inn = 150 * d_v; 

    Ny_max = 0; Nz_max = 0; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            Iz[j * Nz + k] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[j * Nz + k] / d_i - 0.2) / 0.03));       // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[j * Nz + k] / dz[j * Nz + k]);                                       // This finds the intermediate N value with double precision
            Nu_zs[j * Nz + k] = 2 * static_cast<int>(n_int);                                                // This casts the intermediate variable as an integer to be used in convolution coefficient calculation
            
            if (n_int > Nz_max) Nz_max =  2 * static_cast<int>(n_int);                                      // This line checks if the most recently calculated N is bigger than the last one to create boundaries for r_k
        }
    }

    N_holder[0] = Nz_max;       // For ease of storage, I pust Ny and Nz max for u, v, and w in a vector.
    ru_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);  // Create size of random data.

    // Calculate filter coefficients for u' in z-direction
    bu_z = Vector(2 * Nz_max + 1, 0.0); 

    double sum = 0.0;
    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Nz_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 
    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bu_z[i] = exp(-2 * M_PI * abs(kk) / Nz_max) / sum;
    }

    /**
     * Find N for u' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            Iy[j * Nz + k] = 0.67 * Iz[j * Nz + k];
            double n_int = max(1.0, Iy[j * Nz + k] / dy[j * Nz + k]); 

            Nu_ys[j * Nz + k] = static_cast<int>(n_int);
            if (n_int > Ny_max) Ny_max = static_cast<int>(n_int);
        }
    }

    N_holder[1] = Ny_max;
    ru_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);

    // Calculate filter coefficients for u' in y direction
    bu_y = Vector(2 * Ny_max + 1, 0.0);
    sum = 0.0;
    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Ny_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bu_y[i] = exp(-2 * M_PI * abs(kk) / Ny_max) / sum;
    }

    

    /** /////////////////////////////////////////////////////////////////////////////////////////////
    *       Find N for v' filtering in the z-direction.
    */

    Iz_out = 0.3 * d_i;
    Iz_inn = 75 * d_v; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            Iz[j * Nz + k] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[j * Nz + k] / d_i - 0.2) / 0.03));       // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[j * Nz + k] / dz[j * Nz + k]);                                                 // This finds the intermediate N value with double precision
            Nv_zs[j * Nz + k] = static_cast<int>(n_int);                                                     // This casts the intermediate variable as an integer to be used in convolution coefficient calculation
            
            if (n_int > Nz_max) Nz_max = static_cast<int>(n_int);                                           // This line checks if the most recently calculated N is bigger than the last one to create boundaries for r_k
        }
    }

    N_holder[2] = Nz_max; 
    rv_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);  

    // Calculate filter coefficients for v' in z-direction
    bv_z = Vector(2 * Nz_max + 1, 0.0);

    sum = 0.0;
    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Nz_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bv_z[i] = exp(-2 * M_PI * abs(kk) / Nz_max) / sum;
    }

    /**
     * Find N for v' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            Iy[j * Nz + k] = 0.67 * Iz[j * Nz + k];

            double n_int = max(1.0, Iy[j * Nz + k] / dy[j * Nz + k]);  
            Nv_ys[j * Nz + k] = static_cast<int>(n_int);

            if (n_int > Ny_max) Ny_max = static_cast<int>(n_int);
        }
    }

    N_holder[3] = Ny_max; 
    rv_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);

    // Calculate filter coefficients for v' in y-direction
    bv_y = Vector(2 * Ny_max + 1, 0.0);

    sum = 0.0;
    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Ny_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bv_y[i] = exp(-2 * M_PI * abs(kk) / Ny_max) / sum;
    }

    /** /////////////////////////////////////////////////////////////////////////////////////////////////
     * Find N for w' filtering in the z-direction.
     */

    Iz_out = 0.4 * d_i;
    Iz_inn = 150 * d_v; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            Iz[j * Nz + k] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[j * Nz + k]/ d_i - 0.2) / 0.03));       // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[j * Nz + k] / dz[j * Nz + k]);                                                 // This finds the intermediate N value with double precision
            Nw_zs[j * Nz + k] = static_cast<int>(n_int);                                                     // This casts the intermediate variable as an integer to be used in convolution coefficient calculation

            if (n_int > Nz_max) Nz_max = static_cast<int>(n_int);                                           // This line checks if the most recently calculated N is bigger than the last one to create boundaries for r_k
        }
    }
    
    N_holder[4] = Nz_max;
    rw_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0); 

    // Calculate filter coefficients for w' in z-direction
    bw_z = Vector(2 * Nz_max + 1, 0.0);

    sum = 0.0;
    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Nz_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bw_z[i] = exp(-2 * M_PI * abs(kk) / Nz_max) / sum;
    } 

    /**
     * Find N for w' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            Iy[j * Nz + k] = 0.67 * Iz[j * Nz + k];

            double n_int = max(1.0, Iy[j * Nz + k] / dy[j * Nz + k]);  
            Nw_ys[j * Nz + k] = static_cast<int>(n_int);

            if (n_int > Ny_max) Ny_max = static_cast<int>(n_int);
        }
    }

    N_holder[5] = Ny_max;
    rw_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);

    // Calculate filter coefficients for w' in y-direction
    bw_y = Vector(2 * Ny_max + 1, 0.0);

    sum = 0.0;
    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max; 
        double val = exp(-2 * M_PI *  abs(kk) / Ny_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bw_y[i] = exp(-2 * M_PI * abs(kk) / Ny_max) / sum;
    }

}
void DIGITAL_FILTER::filter() {

    // Filter u' in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_holder[1]) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_holder[0]) + N_holder[0] + k);

            double sum = 0.0;
            for (int i = 0; i < 2 * Nu_ys[j * Nz + k] + 1; ++i) {
                int kk = i - Nu_ys[j * Nz + k];
                sum += bu_y[N_holder[1] + kk] * ru_ys[r_idy + kk * Nz];
            }

            ru_zs[r_idz] = sum; 
        }
    }

    // Filter u' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_holder[0]) + N_holder[0] + k);

            double sum = 0.0;
            for (int i = 0; i < 2 * Nu_zs[j * Nz + k] + 1; ++i) {
                int kk = i - Nu_zs[j * Nz + k];
                sum += bu_z[N_holder[0] + kk] * ru_zs[r_idz + kk]; 
            }

            u_fluc[j * Nz + k] = sum;
        }
    }

    // Filter v' in y-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_holder[3]) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_holder[2]) + N_holder[2] + k);

            double sum = 0.0;
            for (int i = 0; i < 2 * Nv_ys[j * Nz + k] + 1; ++i) {
                int kk = i - Nv_ys[j * Nz + k];
                sum += bv_y[N_holder[3] + kk] * rv_ys[r_idy + kk * Nz];
            }

            rv_zs[r_idz] = sum; 
        }
    }

    // Filter v' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_holder[2]) + N_holder[2] + k);

            double sum = 0.0;
            for (int i = 0; i < 2 * Nv_zs[j * Nz + k] + 1; ++i) {
                int kk = i - Nv_zs[j * Nz + k];
                sum += bv_z[N_holder[2] + kk] * rv_zs[r_idz + kk]; 
            }

            v_fluc[j * Nz + k] = sum;
        }
    }

    // Filter w' in y-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_holder[5]) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_holder[4]) + N_holder[4] + k);
            double sum = 0.0;

            for (int i = 0; i < 2 * Nw_ys[j * Nz + k] + 1; ++i) {
                int kk = i - Nw_ys[j * Nz + k];
                sum += bw_y[N_holder[5] + kk] * rw_ys[r_idy + kk * Nz];
            }
     
            rw_zs[r_idz] = sum; 
        }
    }

    // Filter w' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_holder[4]) + N_holder[4] + k);

            double sum = 0.0;
            for (int i = 0; i < 2 * Nw_zs[j * Nz + k] + 1; ++i) {
                int kk = i - Nw_zs[j * Nz + k];
                sum += bw_z[N_holder[4] + kk] * rw_zs[r_idz + kk]; 
            }

            w_fluc[j * Nz + k] = sum;
        }
    }

}

void DIGITAL_FILTER::apply_RST() {
    double u_tau = 0.05;

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            double eta = yc[idx] / d_i;
            double K = 1.5 * u_tau * u_tau * (1.0 - tanh((eta - 0.2) / 0.05));
            double sigma = sqrt(2.0/3.0 * K);

            u_fluc[idx] *= sigma;
            v_fluc[idx] *= sigma;
            w_fluc[idx] *= sigma;

        }
    }
}

void DIGITAL_FILTER::test() {
    fake_grid();
    calculate_filter_properties();
    generate_white_noise();
    filter();
    apply_RST();

    string filename = "velocity_fluc.dat";
    write_tecplot(filename);
}
void DIGITAL_FILTER::fake_grid() {
    double y_max = 0.01;
    double eta = 0.0;  
    double a = 3.0; 

    y = Vector((Ny + 1) * (Nz + 1), 0.0); 
    yc = Vector(Ny * Nz, 0.0);
    z = Vector((Ny + 1) * (Nz + 1), 0.0); 
    dy = Vector(Ny * Nz, 0.0); 
    dz = Vector(Ny * Nz, 0.0);


    for (int j = Ny + 1; j >= 0; --j) {
        for (int k = 0; k < Nz + 1; ++k) {
            eta = ( (j) * y_max / (Ny + 1)) / y_max;   
            y[abs(j - Ny) * (Nz + 1) + k] = y_max * (1 - tanh(a * eta) / tanh(a)); 
            z[j * (Nz + 1) + k] = k * 0.000133;
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            dy[j * Nz + k] = y[(j + 1) * (Nz + 1) + k] - y[j * (Nz + 1) + k];
            dz[j * Nz + k] = 0.000133;        // This is just a holder for dz
            yc[j * Nz + k] = 0.25 * (y[j * (Nz + 1) + k] + y[(j + 1) * (Nz + 1) + k] + y[j * (Nz + 1) + k + 1] + y[(j + 1) * (Nz + 1) + k + 1]);
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
void DIGITAL_FILTER::display_data(Vector& v) {

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            cout << j << ", " << k << ": " << v[j * Nz + k] << endl;
        }
    }

}

void DIGITAL_FILTER::write_tecplot(const string &filename) {

    ofstream file(filename);
    file << "VARIABLES = \"z\", \"y\", \"u_fluc\", \"v_fluc\", \"w_fluc\" \n";
    file << "ZONE T=\"Flow Field\", I=" << Nz + 1 << ", J=" << Ny + 1 << ", F=BLOCK\n";
    file << "VARLOCATION=([3-5]=CELLCENTERED)\n";


    // Loop over y (rows) and z (columns)
    for (int j = 0; j < Ny + 1; ++j) {
        for (int k = 0; k < Nz + 1; ++k) {
            int idx = j * (Nz + 1) + k;  // row-major: y-major
            file << z[idx] << endl;
        }
    }


    for (int j = 0; j < Ny + 1; ++j) {
        for (int k = 0; k < Nz + 1; ++k) {
            int idx = j * (Nz + 1) + k;  // row-major: y-major
            file << y[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << u_fluc[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << v_fluc[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << w_fluc[idx] << endl;
        }
    }

    file.close();
    cout << "Finished plotting." << endl;
}
