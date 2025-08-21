#include "df.hpp"

// Constructor
DIGITAL_FILTER::DIGITAL_FILTER(df_config config) {

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


    // Data structures.
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
    u_filt = Vector(n_cells);
    v_filt = Vector(n_cells);
    w_filt = Vector(n_cells);
    u_filt_old = Vector(n_cells);
    v_filt_old = Vector(n_cells);
    w_filt_old = Vector(n_cells);
    buz_offsets = Vector(n_cells); 
    buy_offsets = Vector(n_cells); 
    bvz_offsets = Vector(n_cells); 
    bvy_offsets = Vector(n_cells); 
    bwz_offsets = Vector(n_cells); 
    bwy_offsets = Vector(n_cells); 

    My = Vector(Ny);
    Uy = Vector(Ny);
    Ty = Vector(Ny);
    rho_fluc = Vector(n_cells);
    T_fluc = Vector(n_cells);
    rms_added = Vector(n_cells, 0.0);
    rms = Vector(n_cells, 0.0);

    calculate_filter_properties(); // Initialize coefficients and filter half-widths. 

    cout << "Calculated filter properties" << endl;
    // First timestep filtering.
    generate_white_noise();
    filtering_sweeps();
    correlate_fields_ts1();
    // string filename = "../files/cpp_vel_fluc_init.dat";
    // write_tecplot(filename);

    // get_rho_T_fluc(); 


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
 *  white noise vectors. It is only called in the constructor and is not called again.
 */
void DIGITAL_FILTER::calculate_filter_properties() {


    // Holders for maximum filter widths for creating ghost cells.
    double Ny_max = 0; Nz_max = 0;
    double n_int, sum, val; 
    int idx, kk, n_val, b_size = 0;

    //=============================================================================================
    // Find filter half-width and convolution coefficients for u' when filtering in the z-direction
    //=============================================================================================


    // Integral length scales
    Iz_out = 0.4 * d_i; 
    Iz_inn = 150 * d_v; 

    // Loop to find integral length scales per cell and get half-widths. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            idx = j * Nz + k;
            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx] / d_i - 0.2) / 0.03));     

            n_int = max(1.0, Iz[idx] / dz[idx]);  
            n_val = 2 * static_cast<int>(n_int);
            Nu_zs[idx] = n_val;      

            b_size += 2 * n_val + 1;
            buz_offsets[idx] = b_size - n_val;

            if (n_val > Nz_max) Nz_max =  n_val;
        }
    }

    // Allocate space for random data and filter coefficient arrays
    N_holder[0] = Nz_max; 
    ru_zs = Vector( (Nz + 2 * Nz_max) * Ny);
    bu_z = Vector(b_size); 

    //This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;       // Cell index.
            
            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nu_zs[idx] + 1; ++i) {
                kk = i - Nu_zs[idx];                        // Offset i-value since sum is from -N to N.
                val = exp(pi_c * abs(kk) / Nu_zs[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nu_zs[idx] + 1; ++i) {
                /**
                 *  Since filter width varies per cell, this index brings us to the actual start of the b vector 
                 *  for this cell. We go up to the start of the vector N_idx, but since each cell gets 2 * Nz_max + 1,
                 *  and most cells dont need this large of a vector, it goes to the center of the b-vector ( + Nz_max).
                 */                 
                kk = i - Nu_zs[idx];    // Offset i-value since sum is from -N to N.
                bu_z[buz_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nu_zs[idx]) / sum; // Final filter coefficient for (i,j).
            }
        }
    }



    //=============================================================================================
    // Find filter half-width and convolution coefficients for u' when filtering in the y-direction
    //=============================================================================================

    b_size =  0;
    // Loop to find integral length scales per cell and get half-widths. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;

            Iy[idx] = 0.67 * Iz[idx];
            n_int = max(1.0, Iy[idx] / dy[idx]); 
            n_val = 2 * static_cast<int>(n_int);

            b_size += 2 * n_val + 1;
            buy_offsets[idx] = b_size - n_val;

            Nu_ys[idx] = n_val;
            if (n_val > Ny_max) Ny_max = n_val;
        }
    }

    // Allocate space for random data and filter coefficient arrays.
    N_holder[1] = Ny_max; 
    ru_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);
    bu_y = Vector(b_size, 0.0);

    //This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;       

            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nu_ys[idx] + 1; ++i) {
                kk = i - Nu_ys[idx];                        
                val = exp(pi_c * abs(kk) / Nu_ys[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nu_ys[idx] + 1; ++i) {
                                
                kk = i - Nu_ys[idx];    
                bu_y[buy_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nu_ys[idx]) / sum; 
            }
        }
    }

    //=============================================================================================
    // Find filter half-width and convolution coefficients for v' when filtering in the z-direction
    //=============================================================================================

    // Reset max values for filter widths.
    Nz_max = 0, Ny_max = 0;
    b_size = 0;

    // Define integral length scales
    Iz_out = 0.3 * d_i;
    Iz_inn = 75 * d_v; 

    // Loop that calculates integral length scales per cell and filter half-widths
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            idx = j * Nz + k;
            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx] / d_i - 0.2) / 0.03));       

            n_int = max(1.0, Iz[idx] / dz[idx]);  
            n_val = 2 * static_cast<int>(n_int);                                               
            Nv_zs[idx] = n_val;
            
            b_size += 2 * n_val + 1;
            bvz_offsets[idx] = b_size - n_val;
            
            if (n_val > Nz_max) Nz_max = n_val;                                           
        }
    }

    // Allocate space for random data and filter coefficient arrays.
    N_holder[2] = Nz_max; 
    rv_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0);  
    bv_z = Vector(b_size, 0.0);

    //This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;       
            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nv_zs[idx] + 1; ++i) {
                kk = i - Nv_zs[idx];                        
                val = exp(pi_c * abs(kk) / Nv_zs[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nv_zs[idx] + 1; ++i) {
                              
                kk = i - Nv_zs[idx];    
                bv_z[bvz_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nv_zs[idx]) / sum; 
            }
        }
    }

    //=============================================================================================
    // Find filter half-width and convolution coefficients for v' when filtering in the y-direction
    //=============================================================================================
    b_size = 0;
    // Loop to find integral length scales per cell and get half-widths.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;
            Iy[idx] = 0.67 * Iz[idx];

            n_int = max(1.0, Iy[idx] / dy[idx]);  
            n_val = 2 * static_cast<int>(n_int);
            Nv_ys[idx] = n_val;

            b_size += 2 * n_val + 1;
            bvy_offsets[idx] = b_size - n_val; 

            if (n_val > Ny_max) Ny_max = n_val;
        }
    }

    // Allocate space for random data and filter coefficient arrays.
    N_holder[3] = Ny_max; 
    rv_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);
    bv_y = Vector(b_size, 0.0);

    // This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;       
            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nv_ys[idx] + 1; ++i) {
                kk = i - Nv_ys[idx];                        
                val = exp(pi_c * abs(kk) / Nv_ys[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nv_ys[idx] + 1; ++i) {
                            
                kk = i - Nv_ys[idx];    
                bv_y[bvy_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nv_ys[idx]) / sum; 
            }
        }
    }

    //=============================================================================================
    // Find filter half-width and convolution coefficients for w' when filtering in the z-direction
    //=============================================================================================

    // Reset max filter widths.
    Nz_max = 0, Ny_max = 0;
    b_size = 0;

    // Define integral length scales.
    Iz_out = 0.4 * d_i;
    Iz_inn = 150 * d_v; 

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) { 

            idx = j * Nz + k;
            Iz[idx] = Iz_inn + (Iz_out - Iz_inn) * 0.5 * (1 + tanh((yc[idx]/ d_i - 0.2) / 0.03));       

            n_int = max(1.0, Iz[idx] / dz[idx]); 
            n_val = 2 * static_cast<int>(n_int);                                             
            Nw_zs[idx] = n_val;             
            
            b_size += 2 * n_val + 1;
            bwz_offsets[idx] = b_size - n_val;

            if (n_val > Nz_max) Nz_max = n_val;                                        
        }
    }
    
    // Allocate space for random data and filter coefficient arrays.
    N_holder[4] = Nz_max;
    rw_zs = Vector( (Nz + 2 * Nz_max) * Ny, 0.0); 
    bw_z = Vector(b_size, 0.0);

    // This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k; 
            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nw_zs[idx] + 1; ++i) {
                kk = i - Nw_zs[idx];                        
                val = exp(pi_c * abs(kk) / Nw_zs[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nw_zs[idx] + 1; ++i) {
                
                kk = i - Nw_zs[idx];    
                bw_z[bwz_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nw_zs[idx]) / sum; 
            }
        }
    }

    //=============================================================================================
    // Find filter half-width and convolution coefficients for w' when filtering in the y-direction
    //=============================================================================================
    b_size = 0;
    // Loop that calculates integral length scales per cell and filter half-widths
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;
            Iy[idx] = 0.67 * Iz[idx];

            n_int = max(1.0, Iy[idx] / dy[idx]);  
            n_val = 2 * static_cast<int>(n_int);
            Nw_ys[idx] = n_val;

            b_size += 2 * n_val + 1;
            bwy_offsets[idx] = b_size - n_val;

            if (n_val > Ny_max) Ny_max = n_val;
        }
    }

    // Allocate space for random data and filter coefficient arrays.
    N_holder[5] = Ny_max;
    rw_ys = Vector(Nz * (2 * Ny_max + Ny), 0.0);
    bw_y = Vector(b_size, 0.0);


    // This loop calculates the vector of filter coefficients for each cell (i,j).
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            idx = j * Nz + k;  
            sum = 0.0;

            // Loop through filter width of cell to compute root sum of intermediate filter coefficients.
            for (int i = 0; i < 2 * Nw_ys[idx] + 1; ++i) {
                kk = i - Nw_ys[idx];                        
                val = exp(pi_c * abs(kk) / Nw_ys[idx]);     
                sum += val *val;
            }
            sum = sqrt(sum);

            // Loop through and calculate final vector of filter coefficients.
            for (int i = 0; i < 2 * Nw_ys[idx] + 1; ++i) {                            
                kk = i - Nw_ys[idx];  
                bw_y[bwy_offsets[idx] + kk] = exp(pi_c * abs(kk) / Nw_ys[idx]) / sum;
            }
        }
    }
}

/** 
 *  This functions filters the random data in y and z seeps using the filter coefficients and filter widths 
 *  previously calculated. The filtering is done in two sweeps, first in y and then in z direction.
 *  */

void DIGITAL_FILTER::filtering_sweeps() {

    int r_idy, r_idz, idx, kk;
    double sum;

    int Nz_pad = N_holder[0];
    int Ny_pad = N_holder[1]; 

    // Filter u' in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idy = ((j + Ny_pad) * Nz + k);                    // Index for r block in y-direction
            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);       // Index for r block in z-direction
            idx = j * Nz + k;                                   // Index for cell (i,j)

            sum = 0.0;
            for (int i = 0; i < 2 * Nu_ys[idx] + 1; ++i) {
                kk = i - Nu_ys[idx];
                sum += bu_y[buy_offsets[idx] + kk] * ru_ys[r_idy + kk * Nz];
            }

            ru_zs[r_idz] = sum; 
        }
    }

    // Filter u' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);
            idx = j * Nz + k;

            sum = 0.0;
            for (int i = 0; i < 2 * Nu_zs[idx] + 1; ++i) {
                kk = i - Nu_zs[idx];
                sum += bu_z[buz_offsets[idx] + kk] * ru_zs[r_idz + kk]; 
            }

            u_filt[idx] = sum;
        }
    }


    Nz_pad = N_holder[2];
    Ny_pad = N_holder[3];

    // Filter v' in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idy = ((j + Ny_pad) * Nz + k);                    // Index for r block in y-direction
            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);       // Index for r block in z-direction
            idx = j * Nz + k;                                   // Index for cell (i,j)

            sum = 0.0;
            for (int i = 0; i < 2 * Nv_ys[idx] + 1; ++i) {
                kk = i - Nv_ys[idx];
                sum += bv_y[bvy_offsets[idx] + kk] * rv_ys[r_idy + kk * Nz];
            }

            rv_zs[r_idz] = sum; 
        }
    }

    // Filter v' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);
            idx = j * Nz + k;

            sum = 0.0;
            for (int i = 0; i < 2 * Nv_zs[idx] + 1; ++i) {
                kk = i - Nv_zs[idx];
                sum += bv_z[bvz_offsets[idx] + kk] * rv_zs[r_idz + kk]; 
            }

            v_filt[idx] = sum;
        }
    }


    Nz_pad = N_holder[4];
    Ny_pad = N_holder[5]; 

    // Filter w' in y-direction. 
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idy = ((j + Ny_pad) * Nz + k);                    // Index for r block in y-direction
            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);       // Index for r block in z-direction
            idx = j * Nz + k;                                   // Index for cell (i,j)

            sum = 0.0;
            for (int i = 0; i < 2 * Nw_ys[idx] + 1; ++i) {
                kk = i - Nw_ys[idx];
                sum += bw_y[bwy_offsets[idx] + kk] * rw_ys[r_idy + kk * Nz];
            }

            rw_zs[r_idz] = sum; 
        }
    }

    // Filter w' in z-direction.
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {

            r_idz = (j * (Nz + 2 * Nz_pad) + Nz_pad + k);
            idx = j * Nz + k;

            sum = 0.0;
            for (int i = 0; i < 2 * Nw_zs[idx] + 1; ++i) {
                kk = i - Nw_zs[idx];                
                sum += bw_z[bwz_offsets[idx] + kk] * rw_zs[r_idz + kk]; 
            }

            w_filt[idx] = sum;
        }
    }
}

/**
 *  This scales the fluctuations by the Reynolds stress tensor R_ij. It is only for the first time step 
 *  so it doesnt correlate the old fluctuations with the new ones. 
 */
void DIGITAL_FILTER::correlate_fields_ts1() {
    
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

            u_fluc[idx] = sqrt(R11[j]) * u_filt[idx];
            v_fluc[idx] = b * u_filt[idx] + sqrt(R22[j] - b * b) * v_filt[idx];
            w_fluc[idx] = sqrt(R33[j]) * w_filt[idx];

        }
    }   

    set_old();
}

/**
 *  This scales the fluctuations by the Reynolds stress tensor R_ij. It then correclates the old fields 
 *  with the new ones by using an integral length/time scale.
 */
void DIGITAL_FILTER::correlate_fields() {
    
    double b;
    double pi = 3.141592654;

    // Timestep correlation for u'
    double Ix = 0.8 * d_i;
    double lt = Ix / U_e;

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            u_filt[idx] = u_filt_old[idx] * exp(-pi * dt / (2.0 * lt)) + u_filt[idx] * sqrt(1.0 - exp(-pi * dt / lt));
        }
    }

    // Timestep correlation for v'
    Ix = 0.3 * d_i;
    lt = Ix / U_e;

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            v_filt[idx] = v_filt_old[idx] * exp(-pi * dt / (2.0 * lt)) + v_filt[idx] * sqrt(1.0 - exp(-pi * dt / lt));
        }
    }

    // Timestep correlation fir w'
    Ix = 0.3 * d_i;
    lt = Ix / U_e;

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            w_filt[idx] = w_filt_old[idx] * exp(-pi * dt / (2.0 * lt)) + w_filt[idx] * sqrt(1.0 - exp(-pi * dt / lt));
        }
    }


    // Scale by Reynolds-stress tensor
    
    for (int j = 0; j < Ny; ++j) {

        if (R11[j] < 1e-10) {
            b = 0.0;
        }       
        else {
            b = R21[j] / sqrt(R11[j]); 
        } 

        for (int k = 0; k < Nz; ++k) {

            int idx = j * Nz + k;            

            u_fluc[idx] = sqrt(R11[j]) * u_filt[idx];
            v_fluc[idx] = b * u_filt[idx] + sqrt(R22[j] - b * b) * v_filt[idx];
            w_fluc[idx] = sqrt(R33[j]) * w_filt[idx];
        }
    }    

    set_old();
}

/**
 * This sets old fluctuations to new ones. 
 */
void DIGITAL_FILTER::set_old() {

    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int idx = j * Nz + k;
            u_filt_old[idx] = u_filt[idx];
            v_filt_old[idx] = v_filt[idx];
            w_filt_old[idx] = w_filt[idx];
        }
    }

}

/** 
 *  This function uses the Strong Reynolds Analogy to find fluctuations for temperature and density assuming 
 *  pressure is constant in the boundary layer. It currently assumes constant gamma = 1.4
 */
void DIGITAL_FILTER::get_rho_T_fluc() {
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {            

            double val = -(1.4 - 1.0) * My[j] * My[j] * Uy[j] * u_fluc[j * Nz + k];
            T_fluc[j * Nz + k] = val; 

            rho_fluc[j * Nz + k] = -val / Ty[j] * rhoy[j]; 

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
    filtering_sweeps();   
    correlate_fields(); 

        
    string filename = "../files/cpp_vel_fluc.dat";
    write_tecplot(filename);
  
    // get_rho_T_fluc(); 
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
             << u_fluc[idx] << " "
             << v_fluc[idx] << " "
             << w_fluc[idx] << endl;
    }

    file.close();
    cout << "Finished plotting line at z = " << Z << "." << endl;
}


void DIGITAL_FILTER::rms_add(Vector& V)  {
    
    rms_counter++;

    for (int i = 0; i < n_cells; ++i) {
        rms_added[i] +=  V[i] * V[i];
    }    
}


void DIGITAL_FILTER::get_rms()  {
    
    for (int i = 0; i < n_cells; ++i) {
        rms[i] = sqrt(rms_added[i] / rms_counter);
    }

}
