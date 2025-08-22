 #include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(df_config config) : u(), v(), w() {

    d_i = config.d_i;
    rho_e = config.rho_e;
    U_e = config.U_e;
    mu_e = config.mu_e;
    grid_file = config.grid_file;
    vel_fluc_file = config.vel_fluc_file;
    vel_file_offset = config.vel_file_offset;
    vel_file_N_values = config.vel_file_N_values;

    rms_counter = 0;
   
    read_grid();        // Eventually this will be important. Currently it just makes up my own grid for testing.
    rho_y_test();       // This will be deleted after I can read in rhoy values and not make some up.  
    get_RST();          // Initializes the Reynold-stress tensor terms

    // Allocate data structures for the filter fields.
    allocate_data_structures(u);
    allocate_data_structures(v);
    allocate_data_structures(w);

    // Integral length scales and filter half-widths.
    u.Iz_out = 0.4 * d_i;
    u.Iz_inn = 150 * d_v;
    u.Lt = 0.8 * d_i / U_e;

    v.Iz_out = 0.3 * d_i;
    v.Iz_inn = 75 * d_v;
    v.Lt = 0.3 * d_i / U_e;

    w.Iz_out = 0.4 * d_i;
    w.Iz_inn = 150 * d_v;
    w.Lt = 0.3 * d_i / U_e;

    // Flow interpolation variables.
    My = Vector(Ny);
    Uy = Vector(Ny);
    Ty = Vector(Ny);
    rho_fluc = Vector(n_cells);
    T_fluc = Vector(n_cells);

    // Initialize coefficients and filter half-widths. 
    calculate_filter_properties(u); 
    calculate_filter_properties(v);
    calculate_filter_properties(w);

    // First timestep filtering.
    generate_white_noise();

    filtering_sweeps(u);
    filtering_sweeps(v);
    filtering_sweeps(w);

    // Only need to apply RST scaling for the first time step. No correlation is done.
    apply_RST_scaling();


    // get_rho_T_fluc(); 


}

void DIGITAL_FILTER::allocate_data_structures(FilterField& F) {
    F.N_ys = vector<int>(n_cells);
    F.N_zs = vector<int>(n_cells);
    F.fluc = Vector(n_cells);
    F.filt = Vector(n_cells);
    F.filt_old = Vector(n_cells);
    F.by_offsets = vector<int>(n_cells);
    F.bz_offsets = vector<int>(n_cells);
}

void DIGITAL_FILTER::allocate_rms_structures(FilterField F) {
    F.rms_added = Vector(n_cells, 0.0);
    F.rms = Vector(n_cells, 0.0);
}

/** 
 *  Somehow, this function reads the grid from a file and populates the y, yc, z, dy, dz vectors. 
 *  I do not know how to go about this part, but I do know that we need to get the vectors y, yc, z,
 *  dy, dz and intgeres Nz, Ny populated with the data from the file. 
*/
void DIGITAL_FILTER::read_grid() {

    Nz = 30;
    Ny = vel_file_N_values;
    n_cells = Nz * Ny;


    y = Vector((Ny + 1) * (Nz + 1), 0.0); 
    z = Vector((Ny + 1) * (Nz + 1), 0.0); 

    yc = Vector(n_cells, 0.0);
    dy = Vector(n_cells, 0.0); 
    dz = Vector(n_cells, 0.0);

    /**
     *  =============================: This part will be removed :=================================
     *  Is just a holder for y and z values since I cant read into grid yet. Needs to be replaced with real grid grabber. 
     */

    double y_max = 2 * d_i;
    double eta = 0.0;  
    double a = 2.0; 

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
            yc[j * Nz + k] = 0.25 * (y[j * (Nz + 1) + k] 
            + y[(j + 1) * (Nz + 1) + k] 
            + y[j * (Nz + 1) + k + 1] 
            + y[(j + 1) * (Nz + 1) + k + 1]);
        }
    }
}

/**
 *  This function is used to test the rho_y vector.
 */
void DIGITAL_FILTER::rho_y_test() {
    rhoy = Vector(Ny, 0.0);
    for (int j = 0; j < Ny; ++j) {
        rhoy[j] = rho_e * (0.7 * yc[j * Nz] / d_i + 0.6);
    }
}

/**
 *  This function generates white noise for the filtering. It uses the PCG random number generator 
 *  to generate normally distributed random numbers.
 */
void DIGITAL_FILTER::generate_white_noise() {

    static pcg32 rng{random_device{}()}; 
    static normal_distribution<> dist(0.0, 1.0);

    auto fill_with_normals = [&](auto &vec) {
        for (auto &val : vec) {
            val = dist(rng);
        }
    };

    fill_with_normals(u.r_ys);
    fill_with_normals(u.r_zs);
    fill_with_normals(v.r_ys);
    fill_with_normals(v.r_zs);
    fill_with_normals(w.r_ys);
    fill_with_normals(w.r_zs);
}

/**
 *  This function uses prescribed integral length scales Iz and Iy to find the filter width
 *  for each cell as well as calculating the convolution coefficients and create the storage for 
 *  white noise vectors. It is only called in the constructor and is not called again.
 */
void DIGITAL_FILTER::calculate_filter_properties(FilterField& F) {

    double n_int, sum, val, Iy; 
    int N, n_val, b_size;
    Vector Iz(n_cells);

    //===================================================================================================
    // Find filter half-width and convolution coefficients for velocity when filtering in the z-direction
    //===================================================================================================

    F.Nz_max = 0;
    F.Ny_max = 0;
    b_size = 0;

    for (int idx = 0; idx < n_cells; ++idx) {

        Iz[idx] = F.Iz_inn + (F.Iz_out - F.Iz_inn) * 0.5 * (1 + tanh((yc[idx] / d_i - 0.2) / 0.03));
        n_int = max(1.0, Iz[idx] / dz[idx]);
        n_val = 2 * static_cast<int>(n_int);
        F.N_zs[idx] = n_val;

        b_size += 2 * n_val + 1;
        F.bz_offsets[idx] = b_size - n_val;
        if (n_val > F.Nz_max) F.Nz_max = n_val;
    }

    // Allocate space for random data and filter coefficient arrays
    F.r_zs = Vector( (Nz + 2 * F.Nz_max) * Ny);
    F.bz = Vector(b_size);

    for (int idx = 0; idx < n_cells; ++idx) {
        sum = 0.0;
        N = F.N_zs[idx];

        // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
        for (int i = -N; i < N; ++i) {                     
            val = exp(pi_c * abs(i) / N);     
            sum += val *val;
        }
        sum = sqrt(sum);

        // Loop through and calculate final vector of filter coefficients.
        for (int i = -N; i < N; ++i) {                 
            F.bz[F.bz_offsets[idx] + i] = exp(pi_c * abs(i) / N) / sum; 
        }
    }            

    //===================================================================================================
    // Find filter half-width and convolution coefficients for velocity when filtering in the y-direction
    //===================================================================================================

    b_size =  0;

    for (int idx = 0; idx < n_cells; ++idx) {
        Iy = 0.67 * Iz[idx];
        n_int = max(1.0, Iy / dy[idx]);
        n_val = 2 * static_cast<int>(n_int);

        b_size += 2 * n_val + 1;
        F.by_offsets[idx] = b_size - n_val;
        F.N_ys[idx] = n_val;
        if (n_val > F.Ny_max) F.Ny_max = n_val;
    }

    F.r_ys = Vector(Nz * (2 * F.Ny_max + Ny)); // Allocate space for random data in y-sweep.
    F.by = Vector(b_size); // Allocate space for filter coefficients in y-sweep.
    // Allocate space for random data and filter coefficient arrays.
   
    for (int idx = 0; idx < n_cells; ++idx) {
        sum = 0.0;
        N = F.N_ys[idx];

        // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
        for (int i = -N; i < N; ++i) {                     
            val = exp(pi_c * abs(i) / N);     
            sum += val *val;
        }
        sum = sqrt(sum);

        // Loop through and calculate final vector of filter coefficients.
        for (int i = -N; i < N; ++i) {                 
            F.by[F.by_offsets[idx] + i] = exp(pi_c * abs(i) / N) / sum; 
        }
    }
}

/** 
 *  This functions filters the random data in y and z seeps using the filter coefficients and filter widths 
 *  previously calculated. The filtering is done in two sweeps, first in y and then in z direction.
 *  */

void DIGITAL_FILTER::filtering_sweeps(FilterField& F) {

    int r_idy, r_idz, idx;
    double sum;

    int Nz_pad = F.Nz_max;
    int Ny_pad = F.Ny_max;

    // Filter in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idy = ((j + Ny_pad) * Nz + k);                    // Index for r block in y-direction
            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);       // Index for r block in z-direction
            idx = j * Nz + k;                                   // Index for cell (i,j)

            sum = 0.0;
            for (int i = -F.N_ys[idx]; i < F.N_ys[idx]; ++i) {
                sum += F.by[F.by_offsets[idx] + i] * F.r_ys[r_idy + i * Nz];
            }

            F.r_zs[r_idz] = sum; 
        }
    }

    // Filter in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);
            idx = j * Nz + k;

            sum = 0.0;
            for (int i = -F.N_zs[idx]; i < F.N_zs[idx]; ++i) {
                sum += F.bz[F.bz_offsets[idx] + i] * F.r_zs[r_idz + i]; 
            }

            F.filt[idx] = sum;
        }
    }
}

/**
 *  This scales the fluctuations by the Reynolds stress tensor R_ij. It is only for the first time step 
 *  so it doesnt correlate the old fluctuations with the new ones. 
 */
void DIGITAL_FILTER::apply_RST_scaling() {
    
    double b;
    double max = 0.0;
    
    for (int j = 0; j < Ny; ++j) {

        if (R11[j] < 1e-10) {
            b = 0.0;
        }       
        else {
            b = R21[j] / sqrt(R11[j]); 
        } 

        for (int k = 0; k < Nz; ++k) {

            int idx = j * Nz + k;            

            u.fluc[idx] = sqrt(R11[j]) * u.filt[idx];
            v.fluc[idx] = b * u.filt[idx] + sqrt(R22[j] - b * b) * v.filt[idx];
            w.fluc[idx] = sqrt(R33[j]) * w.filt[idx];

        }
    }   

    set_old();
}

/**
 *  This scales the fluctuations by the Reynolds stress tensor R_ij. It then correclates the old fields 
 *  with the new ones by using an integral length/time scale.
 */
void DIGITAL_FILTER::correlate_fields(FilterField& F) {
    
    double b;
    double pi = 3.141592654;

    for (int idx = 0; idx < n_cells; ++idx) {
        F.filt[idx] = F.filt_old[idx] * exp(-pi * dt / (2.0 * F.Lt)) + F.filt[idx] * sqrt(1.0 - exp(-pi * dt / F.Lt));
    }

    set_old();
}

/**
 * This sets old fluctuations to new ones. 
 */
void DIGITAL_FILTER::set_old() {

    for (int idx = 0; idx < n_cells; ++idx) { 
        u.filt_old[idx] = u.filt[idx];
        v.filt_old[idx] = v.filt[idx];
        w.filt_old[idx] = w.filt[idx];        
    }

}

/** 
 *  This function uses the Strong Reynolds Analogy to find fluctuations for temperature and density assuming 
 *  pressure is constant in the boundary layer. It currently assumes constant gamma = 1.4
 */
void DIGITAL_FILTER::get_rho_T_fluc() {
    for (int j = 0; j < Ny; ++j) {

        double val = -(1.4 - 1.0) * My[j] * My[j] * Uy[j];

        for (int k = 0; k < Nz; ++k) {            

            double val = -(1.4 - 1.0) * My[j] * My[j] * Uy[j] * u.fluc[j * Nz + k];

            T_fluc[j * Nz + k] = val * u.fluc[j * Nz + k];
            rho_fluc[j * Nz + k] = -val * u.fluc[j * Nz + k] / Ty[j] * rhoy[j]; 

        }
    }
}

/**
 *  This runs all the filtering processes and updates the fluctuations. This is the only function that should 
 *  be called from the main program after initializing the DIGITAL_FILTER object. It generates white noise, filters 
 *  the velocity fluctuations in y and z sweeps, and correlates the fields. @param dt_input is the time step size 
 *  from the CFD simulation.
 */
void DIGITAL_FILTER::filter(double dt_input) {

    dt = dt_input;

    generate_white_noise();    
    filtering_sweeps(u);
    filtering_sweeps(v);
    filtering_sweeps(w);   
    correlate_fields(u);
    correlate_fields(v);
    correlate_fields(w); 
    apply_RST_scaling(); 

    display_data(u.fluc);
        
    string filename = "../files/cpp_vel_fluc.dat";
    write_tecplot(filename);
}

/**
 *  This function is used for testing purposes. It can be used to test the functionality of the DIGITAL_FILTER class.
 */
void DIGITAL_FILTER::test() {


    // generate_white_noise();
    // filtering_sweeps();
    // correlate_fields_ts1();

    // // Write fluctuations to a file.
    // string filename = "velocity_fluc_contour_plot.dat";
    // write_tecplot(filename);

    // filename = "velocity_fluc_line_plot.dat";
    // write_tecplot_line(filename);
}

/**
 *  This function finds the mean and variance of a vector and displays them. It was used to double check the RNG. 
 */
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

/**
 *  This function displays the data in the vector. It is used for debugging purposes.
 */
void DIGITAL_FILTER::display_data(Vector& v) {
    for (auto val : v) {
        std::cout << val << std::endl;
    }
}

/**
 *  This function reads the velocity fluctuation data from a file and calculates the Reynolds stress tensor R_ij. 
 *  It currently uses a number of turbulent boundary layer approximations to estimate the skin friction coefficient 
 *  and the friction velocity. The data is taken from a DNS paper and rescaled to the inflow conditions. It is only 
 *  called once in the constructor to initialize the R_ij tensor.
 */
void DIGITAL_FILTER::get_RST() {
    
    // Open the velocity fluctuation file.
    ifstream fin(vel_fluc_file);
    if (!fin) {
        cerr << "Error opening input file\n";
        return;
    }

    string line;
    // Skip the header lines.
    for (int i = 0; i < vel_file_offset; ++i) {
        getline(fin, line);
    }

    // Create data structures to input file data.
    Vector urms_us(vel_file_N_values), 
           vrms_us(vel_file_N_values), 
           wrms_us(vel_file_N_values), 
           uvrms_us(vel_file_N_values);
           

    int count = 0;

    // Read remaining lines and gather data.
    while (getline(fin, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }

        // columns 9, 10, 11 (1-based), so indices 8, 9, 10 (0-based)  
        if (count < vel_file_N_values) {
            urms_us[count] = values[8];   // urms_utau
            wrms_us[count] = values[9];   // vrms_utau
            vrms_us[count] = values[10];  // wrms_utau
            uvrms_us[count] = values[15]; // upwp_utausq
        }
        count++;       
    }

    
    double val, x_est, Re;      
    Re = rho_e * U_e / mu_e;                                // Re_x / x based on freestream conditions.

    val = d_i / 0.37 * pow(Re, 1.0/5.0);
    x_est = pow(val, 5.0/4.0);                              // Estimated x-location from inflow delta.
    cout << "x_est = " << x_est << endl;
    
    double Cf = 0.0576 / pow(Re * x_est, 1.0/5.0);          // Estimate skin friction coefficient using Prandtl's one-seventh-power law.
    double tau_w = 33.6;                                    // Skin friction
    double u_tau = sqrt(tau_w / 0.0264);
    
    cout << "C_f = " << Cf << endl;
    cout << "tau_w = " << tau_w << endl;
    cout << "u_tau = " << u_tau << endl;
    Vector u_rms(vel_file_N_values), v_rms(vel_file_N_values), w_rms(vel_file_N_values), uv_rms(vel_file_N_values); 

    // Scale u'_rms / u* to inflow u*
    for (int j = 0; j < vel_file_N_values; ++j) {
        u_rms[j] = urms_us[j] * u_tau;
        v_rms[j] = vrms_us[j] * u_tau;
        w_rms[j] = wrms_us[j] * u_tau;
        uv_rms[j] = uvrms_us[j] * u_tau * u_tau;
    }

    R11 = Vector(vel_file_N_values, 0.0);
    R21 = Vector(vel_file_N_values, 0.0);
    R22 = Vector(vel_file_N_values, 0.0);
    R33 = Vector(vel_file_N_values, 0.0);

    // Set Reynolds stress terms
    for (int j = 0; j < vel_file_N_values; ++j) {
        R11[j] = u_rms[j] * u_rms[j];
        R21[j] = uv_rms[j];
        R22[j] = v_rms[j] * v_rms[j];
        R33[j] = w_rms[j] * w_rms[j];
    }

    d_v = 0.0002 * d_i;     // This will need to be chaged.
}

/**
 *  This function writes the velocity fluctuations to a Tecplot file for visualization. It writes the z, y, u_fluc, 
 *  v_fluc, and w_fluc data. @param filename is the name of the file to write to.
 */
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
            file << u.fluc[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << v.fluc[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << w.fluc[idx] << endl;
        }
    }

    file.close();
    cout << "Finished plotting." << endl;
}

/**
 * This function writes the velocity fluctuations at a fixed z-slice to a Tecplot file for visualization.
 */
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
             << u.fluc[idx] << " "
             << v.fluc[idx] << " "
             << w.fluc[idx] << endl;
    }

    file.close();
    cout << "Finished plotting line at z = " << Z << "." << endl;
}

// Main function to calculate and plot the RMS of the velocity fluctuations.
void DIGITAL_FILTER::get_rms()  {
    
    allocate_rms_structures(u);
    allocate_rms_structures(v);
    allocate_rms_structures(w);

    dt = 1e-5;
    for (int i = 0; i < 500; ++i) {
        generate_white_noise();    
        filtering_sweeps(u);
        filtering_sweeps(v);
        filtering_sweeps(w);   
        correlate_fields(u); 
        correlate_fields(v); 
        correlate_fields(w); 
        apply_RST_scaling();
        rms_add();
    }

    plot_rms();
}

// This function adds the squares of the fluctuations to the RMS accumulators.
void DIGITAL_FILTER::rms_add()  {
    
    rms_counter++;

    for (int i = 0; i < n_cells; ++i) {
        u.rms_added[i] +=  u.fluc[i] * u.fluc[i];
        v.rms_added[i] +=  v.fluc[i] * v.fluc[i];
        w.rms_added[i] +=  w.fluc[i] * w.fluc[i];
    }    
}

// This function calculates the RMS of the velocity fluctuations and writes them to a file for visualization.
void DIGITAL_FILTER::plot_rms() {

    for (int i = 0; i < n_cells; ++i) {
        u.rms[i] = sqrt(u.rms_added[i] / rms_counter);
        v.rms[i] = sqrt(v.rms_added[i] / rms_counter);
        w.rms[i] = sqrt(w.rms_added[i] / rms_counter);
    }

    string filename = "../files/cpp_vel_fluc_rms.dat";

    ofstream file(filename);
    file << "VARIABLES = \"z\", \"y\", \"u'_rms\", \"v'_rms\", \"w'_rms\" \n";
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
            file << u.rms[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << v.rms[idx] << endl;
        }
    }


    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;  // row-major: y-major
            file << w.rms[idx] << endl;
        }
    }

    file.close();
    cout << "Finished plotting." << endl;
}

