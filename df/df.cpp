#include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(int Ny, int Nz, double d_i, double d_v, int n_ranks) : Ny(Ny), Nz(Nz), d_i(d_i), d_v(d_v), n_ranks(n_ranks) {


    // int size, rank;
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // MPI_Comm df_comm;

    n_cells = Ny * Nz;  
    timestep = 0;

    fake_grid();

    Iy = Vector(n_cells);
    Iz = Vector(n_cells);
    Nu_ys = vector<int>(n_cells);
    Nu_zs = vector<int>(n_cells);
    Nv_ys = vector<int>(n_cells);
    Nv_zs = vector<int>(n_cells);
    Nw_ys = vector<int>(n_cells);
    Nw_zs = vector<int>(n_cells);
    N_holder = vector<int>(6); 
    u_fluc = Vector(n_cells);
    v_fluc = Vector(n_cells);
    w_fluc = Vector(n_cells);
    u_fluc_old = Vector(n_cells);
    v_fluc_old = Vector(n_cells);
    w_fluc_old = Vector(n_cells);

    calculate_filter_properties();    
}

void DIGITAL_FILTER::generate_white_noise() {

    static pcg32 rng; 
    static normal_distribution<> dist(0.0, 1.0);

    auto fill_with_normals = [&](auto &vec) {
        for (auto &val : vec) {
            val = dist(rng);
        }
    };

    fill_with_normals(ru_ys);
    fill_with_normals(ru_zs);
    fill_with_normals(rv_ys);
    fill_with_normals(rv_zs);
    fill_with_normals(rw_ys);
    fill_with_normals(rw_zs);
}

/**
 *  This function uses prescribed integral length scales Iz and Iy to find the filter width
 *  for each cell as well as calculating the convolution coefficients and create the storage for 
 *  white noise vectors.
 */
void DIGITAL_FILTER::calculate_filter_properties() {

    /**
    *       Find filter half-width (N) for u' when filtering in the z-direction.
    */

    // Integral length scales
    Iz_out = 0.4 * d_i; 
    Iz_inn = 150 * d_v; 

    // Holders for maximum filter widths for creating ghost cells.
    Ny_max = 0; Nz_max = 0; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            int idx = j * Nz + k;

            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx] / d_i - 0.2) / 0.03));      // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[idx] / dz[idx]);                                                 // This finds the intermediate N value with double precision
            Nu_zs[idx] = 2 * static_cast<int>(n_int);                                                   // This casts the intermediate variable as an integer to be used in convolution coefficient calculation
            
            if (n_int > Nz_max) Nz_max =  2 * static_cast<int>(n_int);                                  // This line checks if the most recently calculated N is bigger than the last one to create boundaries for r_k
        }
    }

    N_holder[0] = Nz_max;       // For ease of storage, I pust Ny and Nz max for u, v, and w in a vector.
    ru_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);  // Create size of random data.

    // Calculate filter coefficients for u' in z-direction
    bu_z = Vector(2 * Nz_max + 1, 0.0); 

    double sum = 0.0;
    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max; 
        double val = exp(pi_c * abs(kk) / Nz_max);
        sum += val * val; 
    }
    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bu_z[i] = exp(pi_c * abs(kk) / Nz_max) / sum;
    }

    /**
     * Find N for u' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int idx = j * Nz + k;

            Iy[idx] = 0.67 * Iz[idx];
            double n_int = max(1.0, Iy[idx] / dy[idx]); 

            Nu_ys[idx] = static_cast<int>(n_int);
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
        double val = exp(pi_c *  abs(kk) / Ny_max); 
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bu_y[i] = exp(pi_c * abs(kk) / Ny_max) / sum;
    }

    

    /** /////////////////////////////////////////////////////////////////////////////////////////////
    *       Find N for v' filtering in the z-direction.
    */

    Iz_out = 0.3 * d_i;
    Iz_inn = 75 * d_v; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            int idx = j * Nz + k;

            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx] / d_i - 0.2) / 0.03));       // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[idx] / dz[idx]);                                                 // This finds the intermediate N value with double precision
            Nv_zs[idx] = static_cast<int>(n_int);                                                     // This casts the intermediate variable as an integer to be used in convolution coefficient calculation
            
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
        double val = exp(pi_c *  abs(kk) / Nz_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bv_z[i] = exp(pi_c * abs(kk) / Nz_max) / sum;
    }

    /**
     * Find N for v' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int idx = j * Nz + k;

            Iy[idx] = 0.67 * Iz[idx];

            double n_int = max(1.0, Iy[idx] / dy[idx]);  
            Nv_ys[idx] = static_cast<int>(n_int);

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
        double val = exp(pi_c *  abs(kk) / Ny_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bv_y[i] = exp(pi_c * abs(kk) / Ny_max) / sum;
    }

    /** /////////////////////////////////////////////////////////////////////////////////////////////////
     * Find N for w' filtering in the z-direction.
     */

    Iz_out = 0.4 * d_i;
    Iz_inn = 150 * d_v; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            int idx = j * Nz + k;

            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx]/ d_i - 0.2) / 0.03));       // This line computes the integral length scale dependant on y

            double n_int = max(1.0, Iz[idx] / dz[idx]);                                                 // This finds the intermediate N value with double precision
            Nw_zs[idx] = static_cast<int>(n_int);                                                     // This casts the intermediate variable as an integer to be used in convolution coefficient calculation

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
        double val = exp(pi_c *  abs(kk) / Nz_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Nz_max + 1; ++i) {
        int kk = i - Nz_max;
        bw_z[i] = exp(pi_c * abs(kk) / Nz_max) / sum;
    } 

    /**
     * Find N for w' filtering in the y-direction.
     */

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int idx = j * Nz + k;

            Iy[idx] = 0.67 * Iz[idx];

            double n_int = max(1.0, Iy[idx] / dy[idx]);  
            Nw_ys[idx] = static_cast<int>(n_int);

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
        double val = exp(pi_c *  abs(kk) / Ny_max);
        sum += val * val; 
    }

    sum = sqrt(sum); 

    for (int i = 0; i < 2 * Ny_max + 1; ++i) {
        int kk = i - Ny_max;
        bw_y[i] = exp(pi_c * abs(kk) / Ny_max) / sum;
    }

}
void DIGITAL_FILTER::filtering_sweeps() {

    int N_z = N_holder[0];
    int N_y = N_holder[1];    

    // Filter u' in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_y) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_z) + N_z + k); 
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nu_ys[idx] + 1; ++i) {

                int kk = i - Nu_ys[idx];
                sum += bu_y[N_y + kk] * ru_ys[r_idy + kk * Nz];
            }

            ru_zs[r_idz] = sum; 
        }
    }

    // Filter u' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_z) + N_z + k);
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nu_zs[idx] + 1; ++i) {
                int kk = i - Nu_zs[idx];
                sum += bu_z[N_z + kk] * ru_zs[r_idz + kk]; 
            }

            u_fluc[idx] = sum;
        }
    }

    N_z = N_holder[2];
    N_y = N_holder[3];

    // Filter v' in y-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_y) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_z) + N_z + k);
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nv_ys[idx] + 1; ++i) {
                int kk = i - Nv_ys[idx];
                sum += bv_y[N_y + kk] * rv_ys[r_idy + kk * Nz];
            }

            rv_zs[r_idz] = sum; 
        }
    }

    // Filter v' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_z) + N_z + k);
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nv_zs[idx] + 1; ++i) {
                int kk = i - Nv_zs[idx];
                sum += bv_z[N_z + kk] * rv_zs[r_idz + kk]; 
            }

            v_fluc[idx] = sum;
        }
    }

    N_z = N_holder[4];
    N_y = N_holder[5];

    // Filter w' in y-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idy = ((j + N_y) * Nz + k);
            int r_idz = (j * (Nz + 2 * N_z) + N_z + k); 
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nw_ys[idx] + 1; ++i) {
                int kk = i - Nw_ys[idx];
                sum += bw_y[N_y + kk] * rw_ys[r_idy + kk * Nz];
            }
     
            rw_zs[r_idz] = sum; 
        }
    }

    // Filter w' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            int r_idz = (j * (Nz + 2 * N_z) + N_z + k);
            int idx = j * Nz + k;

            double sum = 0.0;
            for (int i = 0; i < 2 * Nw_zs[idx] + 1; ++i) {

                int kk = i - Nw_zs[idx];
                sum += bw_z[N_z + kk] * rw_zs[r_idz + kk]; 

            }

            w_fluc[idx] = sum;
        }
    }

}

void DIGITAL_FILTER::correlate_fields() {
    
    // Correlate u'
    if (timestep > 0) {
        double Ix = 0.8 * d_i;

        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                int idx = j * Nz + k;
                u_fluc[idx] = u_fluc_old[idx] * 
            }
        }


        // Correlate v' 
        Ix = 0.3 * d_i;

        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                v_fluc
            }
        }

        // Correlate w'
        
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                v_fluc
            }
        }
    }

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            

        }
    }
}

void DIGITAL_FILTER::filter(double dt_input) {

    dt = dt_input;

    if (timestep == 0) {

        generate_white_noise();
        filtering_sweeps();
        correlate_fields(); 

    }

}


void DIGITAL_FILTER::test() {


    auto prog_start = NOW;

    // Time and run white noise generation
    auto start = NOW;
    generate_white_noise();
    auto end = NOW;
    auto elapsed = chrono::duration<double>(end - start); 
    cout << "White noise generation took: " << elapsed.count() << " seconds. " << endl;

    // Time and run filtering operations
    start = NOW;
    filtering_sweeps();
    end = NOW;
    elapsed = chrono::duration<double>(end - start);
    cout << "Filtering sweeps took: " << elapsed.count() << " seconds. " << endl;

    // Time and run Reynolds stress tensor operation
    start = NOW;
    correlate_fields();
    end = NOW;
    elapsed = chrono::duration<double>(end - start);
    cout << "Correlating fields took: " << elapsed.count() << " seconds. " << endl;

    auto prog_end = NOW;
    auto prog_elap = chrono::duration<double>(prog_end - prog_start);
    cout << "Total time: " << prog_elap.count() << " seconds." << endl;

    // Write fluctuations to a file.
    string filename = "velocity_fluc_contour_plot.dat";
    write_tecplot(filename);

    filename = "velocity_fluc_line_plot.dat";
    write_tecplot_line(filename);
}

void DIGITAL_FILTER::fake_grid() {
    double y_max = 0.01;
    double eta = 0.0;  
    double a = 3.0; 

    y = Vector((Ny + 1) * (Nz + 1), 0.0); 
    z = Vector((Ny + 1) * (Nz + 1), 0.0); 

    yc = Vector(Ny * Nz, 0.0);
    dy = Vector(Ny * Nz, 0.0); 
    dz = Vector(Ny * Nz, 0.0);


    for (int j = Ny; j >= 0; --j) {
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
void DIGITAL_FILTER::write_tecplot_line(const string &filename) {
    ofstream file(filename);

    int Z = static_cast<int>(Nz / 2);

    file << "VARIABLES = \"y\", \"u_fluc\", \"v_fluc\", \"w_fluc\"\n";
    file << "ZONE T=\"Line at z=" << Z << "\", I=" << Ny << ", F=POINT\n";
    file << "VARLOCATION=([2-4]=CELLCENTERED)\n";

    // Loop over y (rows) at fixed z = k_slice
    for (int j = 0; j < Ny; ++j) {
        int idx = j * Nz + Z;  // index in 1D row-major array
        file << y[j * (Nz + 1) + Z] << " "   // y-coordinate
             << u_fluc[idx] << " "
             << v_fluc[idx] << " "
             << w_fluc[idx] << endl;
    }

    file.close();
    cout << "Finished plotting line at z = " << Z << "." << endl;
}

void DIGITAL_FILTER::get_RST(const string filename, int N_values, double rho_e, double U_e,  double mu_e, Vector rhoy) {
    
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error opening input file\n";
        return;
    }

    string line;
    // skip first 142 lines (header)
    for (int i = 0; i < 142; ++i) {
        getline(fin, line);
    }


    Vector urms_tau(N_values), vrms_tau(N_values), wrms_tau(N_values);
    int count = 0;

    // read remaining lines (330 rows)
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }

        // columns 9, 10, 11 (1-based), so indices 8, 9, 10 (0-based)
        if (values.size() >= 11) {
            if (count < N_values) {
                urms_tau[count] = values[8];   // urms/utau
                vrms_tau[count] = values[9];   // vrms/utau
                wrms_tau[count] = values[10];  // wrms/utau
                count++;
            } else {
                break; // don't exceed vector size
            }
        }

        count++;        
    }

    double val, x_est, Re;
    Re = rho_e * U_e / mu_e;

    val = d_i / 0.16 * pow(Re, 1.0/7.0);
    x_est = pow(val, 7.0/6.0); 

    double Cf = 0.027 * pow(Re * x_est, -1.0/7.0);
    double tau_w = 0.5 * Cf * rho_e * U_e * U_e; 

    Vector u_rms(N_values), v_rms(N_values), w_rms(N_values); 

    for (int j = 0; j < N_values; ++j) {
        u_rms[j] = urms_tau[j] * sqrt(tau_w / rhoy[j]);
        v_rms[j] = vrms_tau[j] * sqrt(tau_w / rhoy[j]);
        w_rms[j] = wrms_tau[j] * sqrt(tau_w / rhoy[j]);
    }

    R11 = Vector(n_cells, 0.0);
    R21 = Vector(n_cells, 0.0);
    R22 = Vector(n_cells, 0.0);
    R33 = Vector(n_cells, 0.0);

    for (int j = 0; j < Ny; ++j) {

        

    }






}