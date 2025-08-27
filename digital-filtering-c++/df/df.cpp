 #include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(DFConfig config) : u(), v(), w() {

    // Set class members equal to config counterparts
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
    get_RST_in();       // Initializes the Reynold-stress tensor terms
    lerp_RST();         // Interpolates RST values to the CFD grid points.
    plot_RST_lerp();

    // Allocate data structures for the filter fields.
    for (FilterField* F : {&u, &v, &w}) {
        allocate_data_structures(*F);
    }

    // Integral length scales
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
    for (FilterField* F : {&u, &v, &w}) {
        calculate_filter_properties(*F);
    }    

    // First timestep filtering.

    generate_white_noise();
    for (FilterField* F : {&u, &v, &w}) {
        filtering_sweeps(*F);
    }

    apply_RST_scaling();


    // get_rho_T_fluc(); 
}


// ========: Filtering functions :========

void DIGITAL_FILTER::read_grid() {

    Nz = 400;
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

void DIGITAL_FILTER::allocate_data_structures(FilterField& F) {
    F.N_ys = vector<int>(n_cells);
    F.N_zs = vector<int>(n_cells);
    F.fluc = Vector(n_cells);
    F.filt = Vector(n_cells);
    F.filt_old = Vector(n_cells);
    F.by_offsets = vector<int>(n_cells);
    F.bz_offsets = vector<int>(n_cells);
}

void DIGITAL_FILTER::calculate_filter_properties(FilterField& F) {

    double n_int, sum, val, Iy; 
    int N, n_val, b_size, offset;
    Vector Iz(n_cells), temp;

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
        F.bz_offsets[idx] = b_size - n_val - 1;
        if (n_val > F.Nz_max) F.Nz_max = n_val;
    }

    // Allocate space for random data and filter coefficient arrays
    F.r_zs = Vector( (Nz + 2 * F.Nz_max) * Ny);
    F.bz = Vector(b_size);
    temp = Vector(F.Nz_max + 1);

    for (int idx = 0; idx < n_cells; ++idx) {
   
        N = F.N_zs[idx];
        offset = F.bz_offsets[idx];

        sum = 0.0;
        // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
        for (int i = 0; i <= N; ++i) {                     
            temp[i] = exp(pi_c * abs(i) / N);     
            sum += (i==0 ? 1.0 : 2.0) * temp[i] * temp[i];;
        }
        sum = sqrt(sum);

        // Loop through and calculate final vector of filter coefficients.
        for (int i = -N; i <= N; ++i) {                 
            F.bz[offset + i] = temp[abs(i)] / sum; 
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
        F.by_offsets[idx] = b_size - n_val - 1;
        F.N_ys[idx] = n_val;
        if (n_val > F.Ny_max) F.Ny_max = n_val;
    }

    F.r_ys = Vector(Nz * (2 * F.Ny_max + Ny)); // Allocate space for random data in y-sweep.
    F.by = Vector(b_size); // Allocate space for filter coefficients in y-sweep.
    temp = Vector(F.Ny_max + 1);

    for (int idx = 0; idx < n_cells; ++idx) {
        sum = 0.0;
        N = F.N_ys[idx];
        offset = F.by_offsets[idx];

        // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
        for (int i = 0; i <= N; ++i) {                     
            temp[i] = exp(pi_c * abs(i) / N);     
            sum += (i == 0? 1.0 : 2.0) * temp[i] * temp[i];
        }
        sum = sqrt(sum);

        // Loop through and calculate final vector of filter coefficients.
        for (int i = -N; i <= N; ++i) {                 
            F.by[offset + i] = temp[abs(i)] / sum; 
        }
    }
}

void DIGITAL_FILTER::get_RST_in() {
    
    y_d = Vector(vel_file_N_values);  // z locations from the fluctuation file
    R11_in = Vector(vel_file_N_values); 
    R21_in = Vector(vel_file_N_values);
    R22_in = Vector(vel_file_N_values);
    R33_in = Vector(vel_file_N_values);

    R11 = Vector(vel_file_N_values, 0.0);
    R21 = Vector(vel_file_N_values, 0.0);
    R22 = Vector(vel_file_N_values, 0.0);
    R33 = Vector(vel_file_N_values, 0.0);

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
            y_d[count] = values[1];      // z location
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

    double Cf = 0.0576 / pow(Re * x_est, 1.0/5.0);          // Estimate skin friction coefficient using Prandtl's one-seventh-power law.
    double tau_w = 33.6;                                    // Skin friction
    double u_tau = sqrt(tau_w / 0.0264);
    
    Vector u_rms(Ny), v_rms(Ny), w_rms(Ny), uv_rms(Ny); 

    // Scale u'_rms / u* to inflow u*
    for (int j = 0; j < Ny; ++j) {
        u_rms[j] = urms_us[j] * u_tau;
        v_rms[j] = vrms_us[j] * u_tau;
        w_rms[j] = wrms_us[j] * u_tau;
        uv_rms[j] = uvrms_us[j] * u_tau * u_tau;
    }

    // Set Reynolds stress terms
    for (int j = 0; j < Ny; ++j) {
        R11_in[j] = u_rms[j] * u_rms[j];
        R21_in[j] = uv_rms[j];
        R22_in[j] = v_rms[j] * v_rms[j];
        R33_in[j] = w_rms[j] * w_rms[j];
    }


    d_v = 0.0002 * d_i;     // This will need to be chaged.
}

void DIGITAL_FILTER::lerp_RST() {

    for (int j = 0; j < Ny; ++j) {
        // Pick the sample you want; here k=0 along z (cells layout Ny*Nz)
        const std::size_t idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(Nz);
        double yq = yc[idx] / d_i;

        // Below/above range: clamp (or set to zero if you prefer)
        if (yq <= y_d.front()) {
            R11[j] = R11_in.front();
            R21[j] = R21_in.front();
            R22[j] = R22_in.front();
            R33[j] = R33_in.front();
            continue;
        }
        if (yq >= y_d.back()) {
            // choose a policy: clamp or zero. Here we clamp:
            R11[j] = R11_in.back();
            R21[j] = R21_in.back();
            R22[j] = R22_in.back();
            R33[j] = R33_in.back();
            continue;
        }

            // Find first y_d[i1] >= yq
            auto it  = std::lower_bound(y_d.begin(), y_d.end(), yq);
            int i1   = static_cast<int>(it - y_d.begin());
            int i0   = i1 - 1;

            // Guard indices just in case
            if (i0 < 0)            { i0 = 0;   i1 = 1; }
            if (i1 >= Ny)           { i1 = Ny-1; i0 = Ny-2; }

        const double x0 = y_d[i0], x1 = y_d[i1];
        const double dx = x1 - x0;
        const double t  = (dx != 0.0) ? (yq - x0) / dx : 0.0;

        R11[j] = (1.0 - t)*R11_in[i0] + t*R11_in[i1];
        R21[j] = (1.0 - t)*R21_in[i0] + t*R21_in[i1];
        R22[j] = (1.0 - t)*R22_in[i0] + t*R22_in[i1];
        R33[j] = (1.0 - t)*R33_in[i0] + t*R33_in[i1];
    }
}

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

void DIGITAL_FILTER::filtering_sweeps(FilterField& F) {

    int r_idy, r_idz, idx, offset, N;
    double sum;

    int Nz_pad = F.Nz_max;
    int Ny_pad = F.Ny_max;

    // Filter in y-direction. 
    for (int j = 0; j < Ny; ++j) {

        r_idy = (j + Ny_pad) * Nz;
        r_idz = j * (Nz + 2 * Nz_pad) + Nz_pad;
        idx = j * Nz;   
        

        for (int k = 0; k < Nz; ++k) {

            offset = F.by_offsets[idx];
            N = F.N_ys[idx];

            sum = 0.0;
            for (int i = -N; i <= N; ++i) {
                sum += F.by[offset + i] * F.r_ys[r_idy + i * Nz];
            }

            F.r_zs[r_idz] = sum; 

            r_idy++;                // Index for r block in y-direction
            r_idz++;                // Index for r block in z-direction
            idx++;
        }
    }

    // Filter in z-direction.
    for (int j = 0; j < Ny; ++j) {

        r_idz = j * (Nz + 2 * Nz_pad) + Nz_pad;
        idx = j * Nz;

        for (int k = 0; k < Nz; ++k) {
            
            offset = F.bz_offsets[idx];
            N = F.N_zs[idx];

            sum = 0.0;
            for (int i = -N; i <= N; ++i) {
                sum += F.bz[offset + i] * F.r_zs[r_idz + i]; 
            }

            F.filt[idx] = sum;
            r_idz++;
            idx++;
        }
    }
}

void DIGITAL_FILTER::correlate_fields(FilterField& F) {
    
    double b;
    double pi = 3.141592654;
    double alpha = exp(-pi * dt / F.Lt); 

    for (int idx = 0; idx < n_cells; ++idx) {
        F.filt[idx] = F.filt_old[idx] * sqrt(alpha) + F.filt[idx] * sqrt(1.0 - alpha); 
    }
}

void DIGITAL_FILTER::apply_RST_scaling() {
    
    double b;
    
    for (int j = 0; j < Ny; ++j) {

        if (R11[j] < 1e-10) {
            b = 0.0;
        }       
        else {
            b = R21[j] / sqrt(R11[j]); 
        } 

        int idx = j * Nz;

        for (int k = 0; k < Nz; ++k) {         

            u.fluc[idx] = sqrt(R11[j]) * u.filt[idx];
            v.fluc[idx] = b * u.filt[idx] + sqrt(R22[j] - b * b) * v.filt[idx];
            w.fluc[idx] = sqrt(R33[j]) * w.filt[idx];

            for (FilterField* F : {&u, &v, &w}) {
                F->filt_old[idx] = F->filt[idx];
            }

            idx++;
        }
    }   
}

void DIGITAL_FILTER::filter(double dt_input) {

    dt = dt_input;
    auto start = NOW;
    generate_white_noise(); 

    for (FilterField* F : {&u, &v, &w}) {
        filtering_sweeps(*F);
        correlate_fields(*F);
    }

    apply_RST_scaling(); 
    auto end = NOW;
    auto elapsed = chrono::duration<double>(end - start);
    cout << "Filtering took " << elapsed.count() << " seconds." << endl;
    
    string filename = "../files/cpp_vel_fluc.dat";
    write_tecplot(filename);
}

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



// ========: Debugging functions : ========

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

void DIGITAL_FILTER::display_data(Vector& v) {
    for (auto val : v) {
        std::cout << val << std::endl;
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

void DIGITAL_FILTER::rho_y_test() {
    rhoy = Vector(Ny, 0.0);
    for (int j = 0; j < Ny; ++j) {
        rhoy[j] = rho_e * (0.7 * yc[j * Nz] / d_i + 0.6);
    }
}



// ========: RMS functions :=========

void DIGITAL_FILTER::allocate_rms_structures(FilterField& F) {
    F.rms_added = Vector(n_cells, 0.0);
    F.rms = Vector(n_cells, 0.0);
}

void DIGITAL_FILTER::rms_add()  {
    
    rms_counter++;

    for (int idx = 0; idx < n_cells; ++idx) {
        u.rms_added[idx] +=  u.fluc[idx] * u.fluc[idx];
        v.rms_added[idx] +=  v.fluc[idx] * v.fluc[idx];
        w.rms_added[idx] +=  w.fluc[idx] * w.fluc[idx];
    }    
}

void DIGITAL_FILTER::get_rms()  {
    
    allocate_rms_structures(u);
    allocate_rms_structures(v);
    allocate_rms_structures(w);

    dt = 1e-5;

    for (int i = 0; i < 20; ++i) {
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

void DIGITAL_FILTER::plot_rms() {

    for (int i = 0; i < n_cells; ++i) {
        u.rms[i] = sqrt(u.rms_added[i] / rms_counter);
        v.rms[i] = sqrt(v.rms_added[i] / rms_counter);
        w.rms[i] = sqrt(w.rms_added[i] / rms_counter);
    }

    string filename = "../files/cpp_vel_fluc_rms.csv";

    ofstream file(filename);
    // // file << "VARIABLES = \"z\", \"y\", \"u'_rms\", \"v'_rms\", \"w'_rms\" \n";
    // file << "ZONE T=\"Flow Field\", I=" << Nz + 1 << ", J=" << Ny + 1 << ", F=BLOCK\n";
    // file << "VARLOCATION=([3-5]=CELLCENTERED)\n";

    file << "z, y, u'_rms, v'_rms, w'_rms\n";

    // Loop over y (rows) and z (columns)
    for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
            int idx = j * (Nz) + k;  // row-major: y-major
            int iidx = j * (Nz + 1) + k; // for y and z which are one size larger
            file << z[iidx] << ", " << y[iidx] << ", " << u.rms[idx] << ", " << v.rms[idx] << ", " << w.rms[idx] << endl;
        }
    }


    // for (int j = 0; j < Ny + 1; ++j) {
    //     for (int k = 0; k < Nz + 1; ++k) {
    //         int idx = j * (Nz + 1) + k;  // row-major: y-major
    //         file << y[idx] << endl;
    //     }
    // }


    // for (int j = 0; j < Ny; ++j) {
    //     for (int k = 0; k < Nz; ++k) {
    //         int idx = j * Nz + k;  // row-major: y-major
    //         file << u.rms[idx] << endl;
    //     }
    // }


    // for (int j = 0; j < Ny; ++j) {
    //     for (int k = 0; k < Nz; ++k) {
    //         int idx = j * Nz + k;  // row-major: y-major
    //         file << v.rms[idx] << endl;
    //     }
    // }


    // for (int j = 0; j < Ny; ++j) {
    //     for (int k = 0; k < Nz; ++k) {
    //         int idx = j * Nz + k;  // row-major: y-major
    //         file << w.rms[idx] << endl;
    //     }
    // }

    file.close();
    cout << "Finished plotting to file: " << filename << endl;
}

void DIGITAL_FILTER::plot_RST_lerp() {

    string filename = "../files/RST.csv";

    ofstream file(filename);

    file << "y, y_d, R11_in, R21_in, R22_in, R33_in, R11, R21, R22, R33 \n";

    // Loop over y (rows) and z (columns)
    for (int j = 0; j < Ny; ++j) {
        file << yc[j * Nz] << ", " << y_d[j] * 0.0036 << ", " << R11_in[j] << ", " << R21_in[j] << ", " << R22_in[j]
                << ", " << R33_in[j] << ", " << R11[j] << ", " << R21[j] << ", " << R22[j] << ", " << R33[j]  << endl;
    }

    file.close();
    cout << "Finished plotting to file: " << filename << endl;
}



// ========: Plotting functions :=========

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
