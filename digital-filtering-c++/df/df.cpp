 #include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(DFConfig config) : u(), v(), w() {

    // Flow case specific!!
    d_i = 0.0013;
    rho_e = 0.044;
    U_e = 869.1;
    mu = 7.1212e-6;
    T_w = 97.5;
    gcon = 287.0;
    T_e = 55.2;
    P = rho_e * 287.0 * T_e;
    rho_w = 0.0249;
    line_file = "../line.dat";

    /**
     * Things that I need: 
     *  Access to line file (file name)
     *  Ny, Nz, y, dz, dy, yc
     */ 
    
    rms_counter = 0;
   
    read_grid();        // Eventually this will be important. Currently it just makes up my own grid for testing.
    get_RST_in();       // Initializes the Reynold-stress tensor terms and also fills in wall normal profiles for flow variables.

    // Allocate data structures for the filter fields.
    for (FilterField* F : {&u, &v, &w}) {
        allocate_data_structures(*F);
    }

    // Integral length scales (Currently an issue that needs to be resolved)
    u.Iz_out = 0.4 * d_i;
    u.Iz_inn = 150 * d_v;
    u.Lt = 0.8 * d_i / U_e;

    v.Iz_out = 0.3 * d_i;
    v.Iz_inn = 75 * d_v;
    v.Lt = 0.3 * d_i / U_e;

    w.Iz_out = 0.4 * d_i;
    w.Iz_inn = 150 * d_v;
    w.Lt = 0.3 * d_i / U_e;

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


    // // get_rho_T_fluc(); 
}


// ========: Filtering functions :========

void DIGITAL_FILTER::read_grid() {

    Nz = 400;
    Ny = 560; 
    n_cells = Nz * Ny;

    y = Vector((Ny + 1) * (Nz + 1)); 
    z = Vector((Ny + 1) * (Nz + 1)); 

    yc = Vector(n_cells);
    yc_d = Vector(n_cells);
    dy = Vector(n_cells); 
    dz = Vector(n_cells);
    ydline = Vector(Ny);
    yline = Vector(Ny);

    /**
     *  =============================: This part will be removed :=================================
     *  Is just a holder for y and z values since I cant read into grid yet. Needs to be replaced with real grid grabber. 
     */

    double y_max = 3 * d_i;
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
            int idx = j * Nz + k;
            dy[idx] = y[(j + 1) * (Nz + 1) + k] - y[j * (Nz + 1) + k];
            dz[idx] = 0.000133;        // This is just a holder for dz
            yc[idx] = 0.25 * (y[j * (Nz + 1) + k] 
            + y[(j + 1) * (Nz + 1) + k] 
            + y[j * (Nz + 1) + k + 1] 
            + y[(j + 1) * (Nz + 1) + k + 1]);
            yc_d[idx] = yc[idx] / d_i;
        }
        ydline[j] = yc_d[j * Nz];
        yline[j] = yc[j * Nz];
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
    

    // Open RST file. 
    ifstream fin("../files/RST.dat"); 
    if (!fin) {
        cerr << "Error opening input file\n";
        return;
    }
    string line;

    // Skip header line
    getline(fin, line); 
   
    // Find number of values to read in
    getline(fin, line);
    size_t pos = line.find("i=");
    if (pos != std::string::npos) {
        std::istringstream iss(line.substr(pos + 2));
        iss >> N_in;
    }

    // Vector allocations
    y_in = Vector(N_in);
    yin_d = Vector(N_in);
    R11_in = Vector(N_in); 
    R21_in = Vector(N_in);
    R22_in = Vector(N_in);
    R33_in = Vector(N_in);

    Vector urms_us(N_in), 
           vrms_us(N_in), 
           wrms_us(N_in), 
           uvrms_us(N_in);           

    int count = 0;

    // Read lines and put values into data structures
    while (getline(fin, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }

        // columns 9, 10, 11 (1-based), so indices 8, 9, 10 (0-based)  
        if (count < N_in) {
            y_in[count] = values[0];
            yin_d[count] = values[1];      // z location
            urms_us[count] = values[2];   // urms_utau
            vrms_us[count] = values[3];   // vrms_utau
            wrms_us[count] = values[4];  // wrms_utau
            uvrms_us[count] = values[5]; // upwp_utausq
        }
        count++;       
    }

    // This snippet resizes Ny to make sure that it stays within y/delta from RST file. What this 
    // means is that we only add perturbations within ~ 2.5 * d_i.
    int new_Ny = 0;
    int j = 0;
    while (yc_d[j * Nz] <= yin_d[N_in - 1] && j < Ny) {
        new_Ny++;
        j++;
    }
    Ny = new_Ny;


    // Allocating RST vectors
    R11 = Vector(Ny);
    R21 = Vector(Ny);
    R22 = Vector(Ny);
    R33 = Vector(Ny);

    // Resize geometry vectors
    n_cells = Ny * Nz;
    y.resize((Ny + 1) * (Nz + 1));
    yc.resize(n_cells);
    yc_d.resize(n_cells);
    dy.resize(n_cells);
    ydline.resize(Ny);
    yline.resize(Ny);

    // Read in line.dat
    read_line_file();

    Vector u_rms(N_in), v_rms(N_in), w_rms(N_in), uv_rms(N_in); 

    // Scale u'_rms / u* to inflow u*
    for (int j = 0; j < N_in; ++j) {
        R11_in[j] = urms_us[j] * urms_us[j] * u_tau * u_tau;
        R22_in[j] = vrms_us[j] * vrms_us[j] * u_tau * u_tau;
        R33_in[j] = wrms_us[j] * wrms_us[j] * u_tau * u_tau;
        R21_in[j] = uvrms_us[j] * u_tau * u_tau;
    }

    // Interpolate values to our grid values
    R11 = linear_interpolate(yin_d, R11_in, ydline);
    R22 = linear_interpolate(yin_d, R22_in, ydline);
    R21 = linear_interpolate(yin_d, R21_in, ydline);
    R33 = linear_interpolate(yin_d, R33_in, ydline);

    // d_v = mu / (rhos[0] * u_tau); 
    d_v = d_i / 4500; 
    cout << "rho_w = " <<  rhos[0] << endl;
    cout << "d_v = " << d_v << endl;
    cout << "d_i / d_v = " << d_i/d_v << endl;
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
    get_rho_T_fluc();
    auto end = NOW;
    auto elapsed = chrono::duration<double>(end - start);
    cout << "Filtering took " << elapsed.count() << " seconds." << endl;
    
    string filename = "../files/cpp_vel_fluc.csv";
    write_csv(filename);
}

void DIGITAL_FILTER::get_rho_T_fluc() { 

    for (int j = 0; j < Ny; ++j) {

        double temp1 = - 0.5 * (1.4 - 1) * Ms[j] * Ms[j] / Us[j];

        for (int k = 0; k < Nz; ++k) {            

            double temp2 = temp1 * u.fluc[j * Nz + k];

            T_fluc[j * Nz + k] = temp2 * Ts[j];
            rho_fluc[j * Nz + k] = -temp2 * rhos[j]; 

        }
    }
}

void DIGITAL_FILTER::read_line_file() {
    

    // Open RST file. 
    ifstream fin(line_file); 
    if (!fin) {
        cerr << "Error opening input file\n";
        return;
    }
    string line;

    getline(fin, line); // Skip header line
    int N_line;
    // Find number of values
    getline(fin, line);
    size_t pos = line.find("i=");
    if (pos != string::npos) {
        istringstream iss(line.substr(pos + 2));
        iss >> N_line;
    }

    Vector u_file(N_line), p_file(N_line), rho_file(N_line), T_file(N_line), y_file(N_line);
    Us = Vector(Ny);
    Ts = Vector(Ny);
    Ps = Vector(Ny);
    rhos = Vector(Ny);    
    Ms = Vector(Ny); 

    int count = 0;

    // Read remaining lines and gather data.
    while (getline(fin, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        Vector values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }

        // columns 9, 10, 11 (1-based), so indices 8, 9, 10 (0-based)  
        if (count < N_line) {
            y_file[count] = values[1];
            rho_file[count] = values[4];   // urms_utau
            u_file[count] = values[5];   // vrms_utau
            T_file[count] = values[8];  // wrms_utau
            p_file[count] = values[9];  // wrms_utau
        }
        count++;       
    }

    Us = linear_interpolate(y_file, u_file, yline);
    Ps = linear_interpolate(y_file, p_file, yline);
    Ts = linear_interpolate(y_file, T_file, yline);
    rhos = linear_interpolate(y_file, rho_file, yline);
    for (int j = 0; j < Ny; ++j) {
        Ms[j] = Us[j] / sqrt(1.4 * gcon * Ts[j]);
    }

    double dy = y_file[1] - y_file[0];
    double du = Us[1] - Us[0];
    tau_w = mu * du/dy;
    u_tau = sqrt(tau_w / rhos[0]);
    cout << "tau_w = " << tau_w << endl;
    cout << "u_tau = " << u_tau << endl;
}

// ========: Debugging functions : ========

void DIGITAL_FILTER::display_data(Vector& v) {
    for (auto val : v) {
        std::cout << val << std::endl;
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
        T_rms_added[idx] += T_fluc[idx] * T_fluc[idx];
        rho_rms_added[idx] += rho_fluc[idx] * rho_fluc[idx];
    }    
}

void DIGITAL_FILTER::get_rms()  {
    
    allocate_rms_structures(u);
    allocate_rms_structures(v);
    allocate_rms_structures(w);
    T_rms = Vector(n_cells, 0.0);
    rho_rms = Vector(n_cells, 0.0); 
    T_rms_added = Vector(n_cells, 0.0);
    rho_rms_added = Vector(n_cells, 0.0); 

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
        get_rho_T_fluc();
        rms_add();
    }


    plot_rms();
}

void DIGITAL_FILTER::plot_rms() {

    for (int i = 0; i < n_cells; ++i) {
        u.rms[i] = sqrt(u.rms_added[i] / rms_counter);
        v.rms[i] = sqrt(v.rms_added[i] / rms_counter);
        w.rms[i] = sqrt(w.rms_added[i] / rms_counter);
        T_rms[i] = sqrt(T_rms_added[i] / rms_counter);
        rho_rms[i] = sqrt(rho_rms_added[i] / rms_counter);
    }

    string filename = "../files/cpp_vel_fluc_rms.csv";

    ofstream file(filename);
    // // file << "VARIABLES = \"z\", \"y\", \"u'_rms\", \"v'_rms\", \"w'_rms\" \n";
    // file << "ZONE T=\"Flow Field\", I=" << Nz + 1 << ", J=" << Ny + 1 << ", F=BLOCK\n";
    // file << "VARLOCATION=([3-5]=CELLCENTERED)\n";

    file << "z, y, u'_rms, v'_rms, w'_rms, T'_rms, rho'_rms \n";

    // Loop over y (rows) and z (columns)
    for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
            int idx = j * (Nz) + k;  // row-major: y-major
            int iidx = j * (Nz + 1) + k; // for y and z which are one size larger
            file << z[iidx] << ", " << y[iidx] << ", " << u.rms[idx] << ", " << v.rms[idx] << ", " << w.rms[idx] << ", " << T_rms[idx] << ", " << rho_rms[idx] << endl;
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

    ofstream file("../files/myRST.csv");

    file << "y, R11, R21, R22, R33 \n";

    // Loop over y (rows) and z (columns)
    for (int j = 0; j < Ny; ++j) {
        file << yc[j * Nz] << ", " << ", " << R11[j] << ", " << R21[j] << ", " << R22[j] << ", " << R33[j]  << endl;
    }

    file.close();


    ofstream file1("../files/duanRST.csv");
        

    file1 << "y_d, R11_in, R21_in, R22_in, R33_in \n";

    // Loop over y (rows) and z (columns)
    for (int j = 0; j < N_in; ++j) {
        file1 << y_in[j] << ", " << R11_in[j] << ", " << R21_in[j] << ", " << R22_in[j]
                << ", " << R33_in[j] << endl;
    }

    file1.close();


    cout << "Finished plotting RST to file. " << endl; 
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

void DIGITAL_FILTER::write_csv(const std::string &filename) {
    std::ofstream file(filename, std::ios::trunc);
    if (!file) {
        std::cerr << "Error: cannot open " << filename << " for writing.\n";
        return;
    }

    // Header
    file << "z,y,u_fluc,v_fluc,w_fluc,T_fluc,rho_fluc\n";
    file << std::setprecision(15) << std::fixed;

    // Loop over cells (Ny x Nz)
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            // Node indices for the four corners of cell (j,k)
            int n00 =  j      * (Nz + 1) +  k;
            int n01 =  j      * (Nz + 1) + (k + 1);
            int n10 = (j + 1) * (Nz + 1) +  k;
            int n11 = (j + 1) * (Nz + 1) + (k + 1);

            // Cell-center coordinates (average of 4 nodes)
            double yc = 0.25 * (y[n00] + y[n01] + y[n10] + y[n11]);
            double zc = 0.25 * (z[n00] + z[n01] + z[n10] + z[n11]);

            // Cell-centered data index
            int cidx = j * Nz + k;

            // CSV row
            file << zc << "," << yc << ","
                 << u.fluc[cidx] << ","
                 << v.fluc[cidx] << ","
                 << w.fluc[cidx] << "," 
                 << T_fluc[cidx] << ","
                 << rho_fluc[cidx] << "\n";
        }
    }

    file.close();
    std::cout << "CSV written to " << filename << "\n";
}

Vector DIGITAL_FILTER::linear_interpolate(
    const vector<double>& y_data,
    const vector<double>& f_data,
    const vector<double>& y_new)
{
    if (y_data.size() != f_data.size()) {
        throw invalid_argument("y_data and f_data must be the same size.");
    }
    if (y_data.size() < 2) {
        throw invalid_argument("Need at least two data points to interpolate.");
    }

    vector<double> f_new(y_new.size());

    for (size_t j = 0; j < y_new.size(); ++j) {
        double y = y_new[j];

        // Handle out-of-bounds with clamping (or could extrapolate if desired)
        if (y <= y_data.front()) {
            f_new[j] = f_data.front();
            continue;
        }
        if (y >= y_data.back()) {
            f_new[j] = f_data.back();
            continue;
        }

        // Find interval [y_data[i], y_data[i+1]]
        size_t i = 0;
        while (i + 1 < y_data.size() && y > y_data[i+1]) {
            ++i;
        }

        double x0 = y_data[i];
        double x1 = y_data[i+1];
        double f0 = f_data[i];
        double f1 = f_data[i+1];

        // Linear interpolation formula
        f_new[j] = f0 + (f1 - f0) * ( (y - x0) / (x1 - x0) );
    }

    return f_new;
}
