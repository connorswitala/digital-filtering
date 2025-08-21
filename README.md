# Information:

This library contains the class "DIGITAL_FILTER" which was created to use digital filtering to seed inflow with velocity, temperature, and density fluctuations. The code is currently "complete" in the sense that everything is set up to create velocity fluctuations. There are currently data structures that are used that have no values attached to them (or test values added in). Once these are taken care of it should work fine. 

## Directories 
1) "digital-filtering-c++" is the c++ directory for the DIGITAL_FILTERING class that has everything necessary to run digital filtering. Within it are three subdirectories and a makefile. 'df' contains the .cpp and .hpp file. 'pcg-cpp' contains the library that is used for random number generation. 'test' contains a test file for debugging called 'cpp-main.cpp' that should also show how the library is used. There is also a Makefile to build the program in c++
2) "digital-filtering-fortran" is the fortran directory and contains the the DIGITAL_FILTERING module and types. Like above, containts a 'df' and 'test' directory along with a Makefile to build the Fortran files. 
3) "files" is where data files are stored for plots or files to be read in.
4) "programs" is where the built executable files go for both Fortran and C++.
5) "build" is where the .o and .mod files go. There are separate directories for both C++ and Fortran build files.
6) There is a master Makefile that is used in the master directory to build both codes.
7) us3d-user.f90 is an example Fortan code given to me by Graham B.


## Necessary inputs
In its current state, it only takes structured grids with perfect quadrilaterals. When initializing the constructor you will need a number of things. These are taken in a a 'dc_config' type.

First, create the configuration type in Fortran or C++. Here, it will be named 'config'. Then fill in:
1) The inflow boundary layer height: config.d_i 
2) The name of the grid file you want to read it (might end up being deprecated): config.grid_file
3) The name of the file you use to grab velocity fluctuations in order to create the Reynolds stress tensor (RST) to scale the filtered fluctuations: config.vel_fluc_file
4) An integer offset to skip 'n' number of header lines in said fluctuation file: config.vel_file_offset
5) An integer to read in 'n' lines from said fluctuation file: config.vel_file_N_values
6) Freestream values for viscosity, density, and velocity: config.rho_e, .U_e, .mu_e
7) Mean profiles for density, velocity, mach number, and temperature (not finished). 

Once the constructor is set up, a single function needs to be used, "filter(dt_input)". This function needs the timestep for the current CFD iteration to update the fluctuations appropriately. This function does everything necessary to create all fluctuations. The constructor runs a version of this to create fluctuations at the first time step, so it does not need to be called for the intial timestep.


## Current code issues:

Currently, a number of things need to be fixed or double checked. 

1) Mean profiles need to be inserted for boundary layer flow variables (see [7] above). Reading in a .dat file with the mean profiles is probably the best method. 

2) Somehow, this code needs to read in a grid. At the moment, I do not intend to fix this problem as the xCFD code will be run after the US3D test. This will be fixed in US3D first and then we will figure out how to read grid files for the xCFD code. If there are any ideas for how to ready the grid data in Fortran (read in Fluent .cas file, take from US3D, something else), that would be helpful.

3) In the function "get_RST()" I use a number of estimates for turbulent boundary layers. This has an effect on the RST scaling of the fluctuations. I first form a Reynolds number based on free stream without the length, Re = rho_e * U_e / mu_e, where subscipt 'e' refers to freestream edge values. Using this reynolds number I use the formula: 
delta(x) = 0.37 * x / (Re ^ 1/5). Using my input inflow boundary layer height, d_i, I find an estimated x-location that corresponds to this d_i. I then use this x value to find an estimate for the skin friction coefficient, Cf, using Prandtl's one-seventh-power law: Cf = 0.0576 / ( [Re_x]^{1/5}). Cf is then used to find the skin friction 
tau_w = 0.5 * Cf * rho_e * U_e^2. This coefficient is used to find the Morkovin velocity scale of our inflow in order to re-normalize Duan's fluctuations to our particular inflow delta, since they normalize their u'_rms values by their Morkovin velocity scale. Double checking that these methods are valid is important.

4) Once the fluctuations are calculated, I don't know if we want to use this plugin to set inflow values, or just to add the fluctuations to the inflow. For example, is this function called to reset the entire inflow (which would mean that these fluctuations need to be added to the mean values for the inflow and then attach into US3D / xCFD to reset the entire inflow)? Or is this function called and then US3D/xCFD takes the fluctuations and adds to its own mean velocity inflow?

5) Current temperature and density fluctuations invoke the Strong Reynolds Analogy: T' / T = -(gamma - 1) * Ma^2 * u' / U, where Ma^2 = M^2 * U^2 / T. This is taken from page 105 of this document (along with many other of my methods) : https://link.springer.com/article/10.1007/s00162-009-0103-z. I believe that this is all the same as written in Touber's thesis.  I cannot tell if the U and M in this equation is the static inflow velocity U_infty, or mean profile dependant on y in the boundary layer. Using algebra, I believe T' = -(gamma - 1) * M^2 * U * u'. This U shows up in the Lagrangian time scale, so it is important to know if U is a constant or wall-normal dependant. 

6) This code only creates fluctuations up the point that Duan stopped collecting data. That means that fluctuations cannot be created past y/delta = 1.5. I don't know enough about turbulence to state if this is physical, but if fluctuations are needed past this point, we will need to figure out a way to add more fluctuations past this y-point. 

7) I have just updated the way the filter works. I needed to create individual vectors for the filter coeffients **for every cell (i,j)**. I currently have it implemented so that I find the maximum filter width and then the size of the vector for coefficients is (2 * N_max + 1) * (Ny * Nz). However, in some cases, the filter width is quite large. This means that for the cells with large filter widths, it fills in the entire B vector for the cell (i,j) ( cell (i,j) holds 2 * N_max + 1 values in vector B as B is a 1D flattened vector). Unfortunately, this means that whereever a cell has a small filter width and therefore only need ~ 5 filter coefficients calculated, all memory spots in the vector B around these 5 coefficients are set to 0. This is a waste of good memory. It was relatively easy (not really, but **easier**) to do it this way then to dynamically create a B vector of size { (Ny * Nz) * (cell[i,j] filter width) + (Ny * Nz) * (cell[i + 1,j] filter width) + ...}. This would suck to do, but I have ideas on how to fix it if necessary.

8) The viscous length scale d_v needs to be calculated. I have been using an estimation for now, but I should be able to get this value easily once I have all the mean wall-normal flow variables read in.

9) MPI support to distribute result to different cores.

## How the code works:

As previously stated, this library uses digital filtering to simulate inflow turbulence. By inputting everything necessary, it creates turbulence mostly within the specified inflow boundary layer height, d_i. Calling the constructor initializes a number of variables that are used everytime the filtering process is called. It first reads in the velocity fluctuation rms data (currently using Zhang, Duan, and Choudhari's .dat files of which our first case if the M6Tw025 run) in order to form the Reynolds-stress tensor that is used for scaling the filtered data. The constructor also calculates the filter half-width and convolution coefficients. These values stay the same forever, so they do not need to be recalculated after the initialization. 

After reading this and setting the number of points in the y-direction for filtering equal to the number of velocity flucation data points read in from the file, the constructor creates a number of data structures of size Ny * Nz (all cells in spanwise direction are filtered). Somehow, the grid is then read and vertices / cell width and heights are stored in vectors. Note, all vectors are 1D flattened array of size Ny * Nz. 

Six different vectors of convolution coefficients 'b' are created. 2 for each velocity. Each velocity gets one filter width and vector of filter coefficients for each filter direction, y and z. The white noise generation fills in random values with 0 mean and a variance of 1. It then filters these for u', v', w' in the y-direction first with 
r*(k) = sum_(j = -N) ^ (N) b_j *  r_(k + j) for all of the cells, where 'k' denotes cell number. It then sweeps in the z-direction using v(k) = sum_(j = -N) ^ (N) b_j *  
r*_(k + j). The random data is now filtered.


Once the data is filtered, the data needs to be correlated and scaled. It first correlates the fluctuations with the fluctuations from the last time step. Then it scales these fluctuations using values from the RST. Finally, the temperature and density fluctuations are found using these velocity fluctuations using the Strong Reynolds Analogy. 


