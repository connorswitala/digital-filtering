! digital filter_module

module DIGITAL_FILTERING
    implicit none
    private
    public :: digital_filter_type, create_digital_filter, test, filter, df_config

    INTEGER,PARAMETER :: dp = selected_real_kind(15)
    real(kind=dp), PARAMETER :: pi = acos(-1.0_dp)

    ! Define a type (like a class)
    type :: digital_filter_type
        ! Private variables
        character(len=:), allocatable :: grid_file
        character(len=:), allocatable :: vel_fluc_file
        integer :: vel_file_offset
        integer :: vel_file_N_values

        integer :: Ny, Nz, n_cells
        integer :: Ny_max, Nz_max
        real(kind=dp) :: rand1, rand2
        real(kind=dp) :: d_i, d_v
        real(kind=dp) :: rho_e, U_e, mu_e

        real(kind=dp), allocatable :: rhoy(:), Uy(:), My(:), Ty(:)
        real(kind=dp), allocatable :: ru_ys(:), rv_ys(:), rw_ys(:)
        real(kind=dp), allocatable :: ru_zs(:), rv_zs(:), rw_zs(:)
        integer, allocatable :: Nu_ys(:), Nv_ys(:), Nw_ys(:)
        integer, allocatable :: Nu_zs(:), Nv_zs(:), Nw_zs(:)
        real(kind=dp), allocatable :: bu_y(:), bv_y(:), bw_y(:)
        real(kind=dp), allocatable :: bu_z(:), bv_z(:), bw_z(:)

        real(kind=dp), allocatable :: rho_fluc(:), T_fluc(:)
        real(kind=dp), allocatable :: u_fluc(:), v_fluc(:), w_fluc(:)
        real(kind=dp), allocatable :: u_filt(:), v_filt(:), w_filt(:)
        real(kind=dp), allocatable :: u_filt_old(:), v_filt_old(:), w_filt_old(:)
        real(kind=dp), allocatable :: R11(:), R21(:), R22(:), R33(:)

        real(kind=dp) :: dt

        integer, allocatable :: N_holder(:)
        real(kind=dp), allocatable :: y(:), yc(:), z(:), dy(:), dz(:)
        real(kind=dp) :: Iz_inn, Iz_out
        real(kind=dp), allocatable :: Iy(:), Iz(:)
    
    end type digital_filter_type



    type :: df_config
        real(kind=dp) :: d_i, rho_e, U_e, mu_e
        integer :: vel_file_offset, vel_file_N_values
        character(len=256) :: grid_file, vel_fluc_file
    end type df_config

contains

    !   Constructor subroutine

    function create_digital_filter(config) result(obj)
        implicit none
        type(df_config), intent(in) :: config
        character(len=50) :: filename
        type(digital_filter_type) :: obj
        type(digital_filter_type) :: df

        df%d_i = config%d_i
        df%rho_e = config%rho_e
        df%U_e = config%U_e
        df%mu_e = config%mu_e
        df%grid_file = config%grid_file
        df%vel_fluc_file = config%vel_fluc_file
        df%vel_file_offset = config%vel_file_offset
        df%vel_file_N_values = config%vel_file_N_values

        call read_grid(df)
        call rhoy_test(df)
        call get_RST(df)

        allocate(df%Iy(df%n_cells))
        allocate(df%Iz(df%n_cells))
        allocate(df%Nu_ys(df%n_cells))
        allocate(df%Nu_zs(df%n_cells))
        allocate(df%Nv_ys(df%n_cells))
        allocate(df%Nv_zs(df%n_cells))
        allocate(df%Nw_ys(df%n_cells))
        allocate(df%Nw_zs(df%n_cells))
        allocate(df%N_holder(6))
        allocate(df%u_fluc(df%n_cells))
        allocate(df%v_fluc(df%n_cells))
        allocate(df%w_fluc(df%n_cells))
        allocate(df%u_filt(df%n_cells))
        allocate(df%v_filt(df%n_cells))
        allocate(df%w_filt(df%n_cells))
        allocate(df%u_filt_old(df%n_cells))
        allocate(df%v_filt_old(df%n_cells))
        allocate(df%w_filt_old(df%n_cells))
        allocate(df%My(df%Ny))
        allocate(df%Uy(df%Ny))
        allocate(df%Ty(df%Ny))
        allocate(df%rho_fluc(df%n_cells))
        allocate(df%T_fluc(df%n_cells))

     
        call calculate_filter_properties(df)

        ! First timestep
        call generate_white_noise(df)
        call filtering_sweeps(df)
        call correlate_fields_ts1(df)
        filename = "../files/initial-velocity-fluctuations.csv"
        call write_csv(df, filename)


        obj = df
    end function create_digital_filter


    !   This subroutine displays the array that is input. Used for bug checking.

    subroutine display_data(array)
            implicit none
            real(kind=dp), intent(in) :: array(:)       ! the array to display
            integer :: i

            do i = 1, size(array)
                print *, array(i)
            end do
    end subroutine display_data


    !   Somehow, this subroutine reads the grid from a file and populates the y, yc, z, dy, dz vectors. 
    !   I do not know how to go about this part, but I do know that we need to get the vectors y, yc, z,
    !   dy, dz and intgeres Nz, Ny populated with the data from the file. 
 
    subroutine read_grid(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp) :: y_max, eta, a
        integer :: j, k, idx

        ! set sizes
        df%Nz = 80
        df%Ny = df%vel_file_N_values
        df%n_cells = df%Nz * df%Ny

        !allocate arrays
        allocate(df%y( (df%Ny + 1) * (df%Nz + 1)))
        allocate(df%z( (df%Ny + 1) * (df%Nz + 1)))
        allocate(df%yc( (df%n_cells)))
        allocate(df%dy( (df%n_cells)))
        allocate(df%dz( (df%n_cells)))

        y_max = 0.01
        eta = 0.0
        a = 3.0
        

        !--- First loop: y and z
        do j = df%Ny, 1, -1
            do k = 1, df%Nz + 1
                eta = real(j, dp) / real(df%Ny + 1, dp)
                df%y((df%Ny - j) * (df%Nz + 1) + k) = &
                    y_max * (1.0_dp - tanh(a * eta) / tanh(a))
                df%z((j - 1) * (df%Nz + 1) + k) = k * 0.00133_dp
            end do
        end do

        !--- Second loop: dy, dz, yc
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                df%dy(idx) = df%y((j + 0) * (df%Nz + 1) + k) - df%y((j - 1) * (df%Nz + 1) + k)
                df%dz(idx) = 0.000133_dp
                df%yc(idx) = 0.25_dp * ( df%y((j - 1)*(df%Nz + 1) + k) + df%y(j*(df%Nz + 1) + k) &
                                    + df%y((j - 1)*(df%Nz + 1) + k + 1) + df%y(j*(df%Nz + 1) + k + 1) )
            end do
        end do
    end subroutine read_grid


    !   This subroutine is used to test the rho_y vector.

    subroutine rhoy_test(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        integer j

        allocate(df%rhoy(df%Ny))       

        do j = 1, df%Ny
            df%rhoy(j) = df%rho_e * (0.7_dp * df%yc(j * df%Nz) / df%d_i + 0.5_dp)
        end do
    end subroutine rhoy_test


    !   This subroutine generates white noise for the filtering. It uses the PCG random number generator 
    !   to generate normally distributed random numbers.

    subroutine generate_white_noise(df)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        type(digital_filter_type), intent(inout) :: df

        call fill_with_noise(df%ru_ys)
        call fill_with_noise(df%ru_zs)
        call fill_with_noise(df%rv_ys)
        call fill_with_noise(df%rv_zs)
        call fill_with_noise(df%rw_ys)
        call fill_with_noise(df%rw_zs)
 
        contains
            subroutine fill_with_noise(vec)
                real(kind=dp), intent(inout) :: vec(:)
                integer :: j
                real(kind=dp) :: u1, u2, r, theta, z1, z2

                j = 1
                do while (j <= size(vec))
                    call random_number(u1)
                    call random_number(u2)
                    if (u1 == 0.0_dp) u1 = 1.0e-12_dp 

                    r = sqrt(-2.0_dp * log(u1))
                    theta = 2.0_dp * acos(-1.0_dp) * u2
                    z1 = r * cos(theta)
                    z2 = r * sin(theta)

                    vec(j) = z1
                    if (j + 1 <= size(vec)) vec(j + 1) = z2
                    j = j + 2
                end do
            end subroutine fill_with_noise

    end subroutine generate_white_noise
    

    !   This subroutine uses prescribed integral length scales Iz and Iy to find the filter width
    !   for each cell as well as calculating the convolution coefficients and create the storage for 
    !   white noise vectors. It is only called in the constructor and is not called again.

    subroutine calculate_filter_properties(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
    

        integer Ny_max, Nz_max, i, j, k, idx, N_idx, kk, offset, depth, n_val
        real(kind=dp) n_int, sum, val, Iz_out, Iz_inn


        !=============================================================================================
        ! Find filter half-width and convolution coefficients for u' when filtering in the z-direction
        !=============================================================================================

        ! Define integral length scales
        Iz_out = 0.4_dp * df%d_i
        Iz_inn = 150.0_dp * df%d_v

        ! Holdering for maximum filter widths. Used for creating ghost cells.
        Ny_max = 0
        Nz_max = 0 

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%Iz(idx) = Iz_inn + (Iz_out - Iz_inn) * 0.5_dp * (1.0_dp + tanh((df%yc(idx) / df%d_i - 0.2_dp) / 0.03_dp))

                n_int = max(1.0_dp, df%Iz(idx) / df%dz(idx))
                n_val = 2 * int(n_int);
                df%Nu_zs(idx) = n_val

                if (n_val > Nz_max) Nz_max = n_val

            end do  
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(1) = Nz_max
        depth = 2 * Nz_max + 1
        allocate(df%ru_zs((df%Nz + 2 * Nz_max) * df%Ny))
        allocate(df%bu_z(depth * df%n_cells))              

        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp
                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nu_zs(idx) + 1
                    kk = (i - 1) - df%Nu_zs(idx)          ! Offset i-value since sum is from -N to N.
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nu_zs(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nu_zs(idx) + 1

                 ! Since filter width varies per cell, this index brings us to the actual start of the b vector 
                 ! for this cell. We go up to the start of the vector N_idx, but since each cell gets 2 * Nz_max + 1,
                 ! and most cells dont need this large of a vector, it goes to the center of the b-vector ( + Nz_max).
                 offset = N_idx + Nz_max
                 kk = (i - 1) - df%Nu_zs(idx)
                 df%bu_z(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nu_zs(idx), dp)) / sum
                end do
            end do
        end do


        !=============================================================================================
        ! Find filter half-width and convolution coefficients for u' when filtering in the y-direction
        !=============================================================================================

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k

                df%Iy(idx) = 0.67_dp * df%Iz(idx)
                n_int = max(1.0_dp, df%Iy(idx) / df%dy(idx))
                n_val = 2 * int(n_int)

                df%Nu_ys(idx) = n_val
                if (n_val > Ny_max) Ny_max = n_val
            end do
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(2) = Ny_max
        depth = 2 * Ny_max + 1
        allocate(df%ru_ys(df%Nz * (2 * Ny_max + df%Ny)))
        allocate(df%bu_y(depth * df%n_cells))


        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df% Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp

                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nu_ys(idx) + 1
                    kk = (i - 1) - df%Nu_ys(idx)
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nu_ys(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nu_ys(idx) + 1
                    offset = N_idx + Ny_max
                    kk = (i - 1) - df%Nu_ys(idx)
                    df%bu_y(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nu_ys(idx), dp)) / sum
                end do
            end do
        end do

        !=============================================================================================
        ! Find filter half-width and convolution coefficients for v' when filtering in the z-direction
        !=============================================================================================

        ! Define integral length scales
        Iz_out = 0.3_dp * df%d_i
        Iz_inn = 75.0_dp * df%d_v

        ! Holdering for maximum filter widths. Used for creating ghost cells.
        Ny_max = 0
        Nz_max = 0 

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%Iz(idx) = Iz_inn + (Iz_out - Iz_inn) * 0.5_dp * (1.0_dp + tanh((df%yc(idx) / df%d_i - 0.2_dp) / 0.03_dp))

                n_int = max(1.0_dp, df%Iz(idx) / df%dz(idx))
                n_val = 2 * int(n_int);
                df%Nv_zs(idx) = n_val

                if (n_val > Nz_max) Nz_max = n_val

            end do  
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(3) = Nz_max
        depth = 2 * Nz_max + 1
        allocate(df%rv_zs((df%Nz + 2 * Nz_max) * df%Ny))
        allocate(df%bv_z(depth * df%n_cells))              

        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp
                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nv_zs(idx) + 1
                    kk = (i - 1) - df%Nv_zs(idx)          ! Offset i-value since sum is from -N to N.
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nv_zs(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nv_zs(idx) + 1

                 ! Since filter width varies per cell, this index brings us to the actual start of the b vector 
                 ! for this cell. We go up to the start of the vector N_idx, but since each cell gets 2 * Nz_max + 1,
                 ! and most cells dont need this large of a vector, it goes to the center of the b-vector ( + Nz_max).
                 offset = N_idx + Nz_max
                 kk = (i - 1) - df%Nv_zs(idx)
                 df%bv_z(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nv_zs(idx), dp)) / sum
                end do
            end do
        end do


        !=============================================================================================
        ! Find filter half-width and convolution coefficients for v' when filtering in the y-direction
        !=============================================================================================

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k

                df%Iy(idx) = 0.67_dp * df%Iz(idx)
                n_int = max(1.0_dp, df%Iy(idx) / df%dy(idx))
                n_val = 2 * int(n_int)

                df%Nv_ys(idx) = n_val
                if (n_val > Ny_max) Ny_max = n_val
            end do
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(4) = Ny_max
        depth = 2 * Ny_max + 1
        allocate(df%rv_ys(df%Nz * (2 * Ny_max + df%Ny)))
        allocate(df%bv_y(depth * df%n_cells))


        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df% Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp

                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nv_ys(idx) + 1
                    kk = (i - 1) - df%Nv_ys(idx)
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nv_ys(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nv_ys(idx) + 1
                    offset = N_idx + Ny_max
                    kk = (i - 1) - df%Nv_ys(idx)
                    df%bv_y(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nv_ys(idx), dp)) / sum
                end do
            end do
        end do


        !=============================================================================================
        ! Find filter half-width and convolution coefficients for w' when filtering in the z-direction
        !=============================================================================================

        ! Define integral length scales
        Iz_out = 0.4_dp * df%d_i
        Iz_inn = 150.0_dp * df%d_v

        ! Holdering for maximum filter widths. Used for creating ghost cells.
        Ny_max = 0
        Nz_max = 0 

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%Iz(idx) = Iz_inn + (Iz_out - Iz_inn) * 0.5_dp * (1.0_dp + tanh((df%yc(idx) / df%d_i - 0.2_dp) / 0.03_dp))

                n_int = max(1.0_dp, df%Iz(idx) / df%dz(idx))
                n_val = 2 * int(n_int);
                df%Nw_zs(idx) = n_val

                if (n_val > Nz_max) Nz_max = n_val

            end do  
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(5) = Nz_max
        depth = 2 * Nz_max + 1
        allocate(df%rw_zs((df%Nz + 2 * Nz_max) * df%Ny))
        allocate(df%bw_z(depth * df%n_cells))              

        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp
                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nw_zs(idx) + 1
                    kk = (i - 1) - df%Nw_zs(idx)          ! Offset i-value since sum is from -N to N.
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nw_zs(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nw_zs(idx) + 1

                 ! Since filter width varies per cell, this index brings us to the actual start of the b vector 
                 ! for this cell. We go up to the start of the vector N_idx, but since each cell gets 2 * Nz_max + 1,
                 ! and most cells dont need this large of a vector, it goes to the center of the b-vector ( + Nz_max).
                 offset = N_idx + Nz_max
                 kk = (i - 1) - df%Nw_zs(idx)
                 df%bw_z(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nw_zs(idx), dp)) / sum
                end do
            end do
        end do



        !=============================================================================================
        ! Find filter half-width and convolution coefficients for w' when filtering in the y-direction
        !=============================================================================================

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k

                df%Iy(idx) = 0.67_dp * df%Iz(idx)
                n_int = max(1.0_dp, df%Iy(idx) / df%dy(idx))
                n_val = 2 * int(n_int)

                df%Nw_ys(idx) = n_val
                if (n_val > Ny_max) Ny_max = n_val
            end do
        end do

        ! Allocate space for random data and filter coefficient arrays.
        df%N_holder(6) = Ny_max
        depth = 2 * Ny_max + 1
        allocate(df%rw_ys(df%Nz * (2 * Ny_max + df%Ny)))
        allocate(df%bw_y(depth * df%n_cells))


        ! This loop calculates the vector of filter coefficients for each cell (i,j).
        do j = 1, df%Ny
            do k = 1, df% Nz

                idx = (j - 1) * df%Nz + k       ! Cell index
                N_idx = (idx - 1) * depth + 1   ! Starting index for b_vector of cell (i,j)

                sum = 0.0_dp

                ! Loop through filter width of cell to compute root sum of intermediate filter coefficients.
                do i = 1, 2 * df%Nw_ys(idx) + 1
                    kk = (i - 1) - df%Nw_ys(idx)
                    val = exp(-2.0_dp * pi * abs(kk) / real(df%Nw_ys(idx), dp))
                    sum = sum + val * val
                end do
                sum = sqrt(sum)

                ! Loop through and calculate final vector of filter coefficients.
                do i = 1, 2 * df%Nw_ys(idx) + 1
                    offset = N_idx + Ny_max
                    kk = (i - 1) - df%Nw_ys(idx)
                    df%bw_y(offset + kk) = exp(-2.0_dp * pi * abs(kk) / real(df%Nw_ys(idx), dp)) / sum
                end do
            end do
        end do
       
    end subroutine calculate_filter_properties


    !  This scales the fluctuations by the Reynolds stress tensor R_ij. It is only for the first time step 
    !  so it doesnt correlate the old fluctuations with the new ones. 

    subroutine correlate_fields_ts1(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp) b
        integer :: k, j, idx

        do j = 1, df%Ny

            if (df%R11(j) < 1e-10) then
                b = 0.0
            else
                b = df%R21(j) / sqrt(df%R11(j))
            end if

            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k

                df%u_fluc(idx) = sqrt(df%R11(j)) * df%u_filt(idx)
                df%v_fluc(idx) = b * df%u_filt(idx) + sqrt(df%R22(j) - b**2) * df%v_filt(idx)
                df%w_fluc(idx) = sqrt(df%R33(j)) * df%w_filt(idx)
            end do
        end do
        
        call set_old(df)
    end subroutine correlate_fields_ts1


    !  This scales the fluctuations by the Reynolds stress tensor R_ij. It then correclates the old fields 
    !  with the new ones by using an integral length/time scale.

    subroutine correlate_fields(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp) b, Ix, lt
        integer :: j, k, idx

        ! Timestep correlation for u'
        Ix = 0.8_dp * df%d_i
        lt = Ix / df%U_e

        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%u_filt(idx) = df%u_filt_old(idx) * exp(-pi * df%dt / (2.0_dp * lt)) + df%u_filt(idx) * sqrt(1.0_dp - exp(-pi * df%dt / lt))
            end do
        end do


        ! Timestep correlation for v'
        Ix = 0.3_dp * df%d_i
        lt = Ix / df%U_e

        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%v_filt(idx) = df%v_filt_old(idx) * exp(-pi * df%dt / (2.0_dp * lt)) + df%v_filt(idx) * sqrt(1.0_dp - exp(-pi * df%dt / lt))
            end do
        end do


        ! Timestep correlation for w'
        Ix = 0.8_dp * df%d_i
        lt = Ix / df%U_e

        do j = 1, df%Ny
            do k = 1, df%Nz

                idx = (j - 1) * df%Nz + k
                df%w_filt(idx) = df%w_filt_old(idx) * exp(-pi * df%dt / (2.0_dp * lt)) + df%w_filt(idx) * sqrt(1.0_dp - exp(-pi * df%dt / lt))
            end do
        end do


        ! Scale by Reynolds stress tensor
        do j = 1, df%Ny

            if (df%R11(j) < 1e-10) then
                b = 0.0
            else
                b = df%R21(j) / sqrt(df%R11(j))
            end if

            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k

                df%u_fluc(idx) = sqrt(df%R11(j)) * df%u_filt(idx)
                df%v_fluc(idx) = b * df%u_filt(idx) + sqrt(df%R22(j) - b**2) * df%v_filt(idx)
                df%w_fluc(idx) = sqrt(df%R33(j)) * df%w_filt(idx)
            end do
        end do
        
        call set_old(df)

    end subroutine correlate_fields

    !  This functions filters the random data in y and z seeps using the filter coefficients and filter widths 
    !  previously calculated. The filtering is done in two sweeps, first in y and then in z direction.

    subroutine filtering_sweeps(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        integer :: r_idy, r_idz, idx, N_idx, offset, kk, Ny_pad, Nz_pad, i, j, k
        real(kind=dp) :: sum

        ! ========================================
        ! ------: Filter u' in y-direction :------
        ! ========================================
        
        Nz_pad = df%N_holder(1)
        Ny_pad = df%N_holder(2)

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idy = ((j - 1) + Ny_pad) * df%Nz + k
                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Ny_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nu_ys(idx) + 1
                    kk = (i - 1) - df%Nu_ys(idx)
                    offset = N_idx + Ny_pad
                    sum = sum + df%bu_y(offset + kk) * df%ru_ys(r_idy + kk * df%Nz)
                end do

                df%ru_zs(r_idz) = sum

            end do
        end do


        ! ========================================
        ! ------: Filter u' in z-direction :------
        ! ========================================

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Nz_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nu_zs(idx) + 1
                    kk = (i - 1) - df%Nu_zs(idx)
                    offset = N_idx + Nz_pad
                    sum = sum + df%bu_z(offset + kk) * df%ru_zs(r_idz + kk)
                end do

                df%u_filt(idx) = sum

            end do
        end do


        ! ========================================
        ! ------: Filter v' in y-direction :------
        ! ========================================
        
        Nz_pad = df%N_holder(3)
        Ny_pad = df%N_holder(4)

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idy = ((j - 1) + Ny_pad) * df%Nz + k
                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Ny_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nv_ys(idx) + 1
                    kk = (i - 1) - df%Nv_ys(idx)
                    offset = N_idx + Ny_pad
                    sum = sum + df%bv_y(offset + kk) * df%rv_ys(r_idy + kk * df%Nz)
                end do

                df%rv_zs(r_idz) = sum

            end do
        end do
        

        ! ========================================
        ! ------: Filter v' in z-direction :------
        ! ========================================

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Nz_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nv_zs(idx) + 1
                    kk = (i - 1) - df%Nv_zs(idx)
                    offset = N_idx + Nz_pad
                    sum = sum + df%bv_z(offset + kk) * df%rv_zs(r_idz + kk)
                end do

                df%v_filt(idx) = sum

            end do
        end do


        ! ========================================
        ! ------: Filter w' in y-direction :------
        ! ========================================
        
        Nz_pad = df%N_holder(5)
        Ny_pad = df%N_holder(6)

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idy = ((j - 1) + Ny_pad) * df%Nz + k
                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Ny_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nw_ys(idx) + 1
                    kk = (i - 1) - df%Nw_ys(idx)
                    offset = N_idx + Ny_pad
                    sum = sum + df%bw_y(offset + kk) * df%rw_ys(r_idy + kk * df%Nz)
                end do

                df%rw_zs(r_idz) = sum

            end do
        end do

        ! ========================================
        ! ------: Filter w' in z-direction :------
        ! ========================================

        do j = 1, df%Ny
            do k = 1, df%Nz

                r_idz = (j - 1) * (df%Nz + 2 * Nz_pad) + Nz_pad + k
                idx = (j - 1) * df%Nz + k
                N_idx = (idx - 1) * (2 * Nz_pad + 1) + 1

                sum = 0.0
                do i = 1, 2 * df%Nw_zs(idx) + 1
                    kk = (i - 1) - df%Nw_zs(idx)
                    offset = N_idx + Nz_pad
                    sum = sum + df%bw_z(offset + kk) * df%rw_zs(r_idz + kk)
                end do

                df%w_filt(idx) = sum

            end do
        end do

    end subroutine filtering_sweeps


    !  This runs all the filtering processes and updates the fluctuations. This is the only function that should 
    !  be called from the main program after initializing the DIGITAL_FILTER object. It generates white noise, filters 
    !  the velocity fluctuations in y and z sweeps, and correlates the fields. @param dt_input is the time step size 
    !  from the CFD simulation.

    subroutine filter(df, dt_input)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp), intent(in) :: dt_input
        real(kind=dp) :: start, end, elapsed
        character(len=50) :: filename
        integer :: i

        df%dt = dt_input

        call generate_white_noise(df)   
        call filtering_sweeps(df)
        call correlate_fields(df)     
        filename = "../files/velocity-fluctuations_+dt.csv"
        call write_csv(df, filename)


    end subroutine filter


    !  This function uses the Strong Reynolds Analogy to find fluctuations for temperature and density assuming 
    !  pressure is constant in the boundary layer. It currently assumes constant gamma = 1.4

    subroutine get_rho_T_fluc(df)
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp) :: val
        integer :: j, k, idx

        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                val = (-1.4_dp - 1.0_dp) * df%My(j)**2 * df%Uy(j) * df%u_fluc(idx)
                df%T_fluc(idx) = val
                df%rho_fluc(idx) = -val / df%Ty(j) * df%rhoy(j)
            end do
        end do
        
    end subroutine get_rho_T_fluc

    !   This sets old fluctuations to new ones. 
    subroutine set_old(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        integer :: j, k, idx

        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k

                df%u_filt_old(idx) = df%u_filt(idx)
                df%v_filt_old(idx) = df%v_filt(idx)
                df%w_filt_old(idx) = df%w_filt(idx)
            end do  
        end do

    end subroutine set_old

    subroutine test(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        character(len=50) :: filename
        filename = "../../files/fluc.csv"
        call generate_white_noise(df)
        call filtering_sweeps(df)
        call correlate_fields_ts1(df)
        call write_csv(df, filename)    
    end subroutine test

    subroutine find_mean_variance(df, v)
        implicit none
        type(digital_filter_type), intent(inout) :: df

        real(kind=dp), intent(in) :: v(df%n_cells)
        integer i  
        real(kind=dp) :: mean_sum, var_sum, mean, variance
    
        mean_sum = 0.0_dp
        var_sum = 0.0_dp

        do i = 1, df%n_cells
            mean_sum = mean_sum + v(i)
        end do

        mean = mean_sum / real(df%n_cells, dp)

        do i = 1, df%n_cells
            var_sum = var_sum + (v(i) - mean)**2
        end do
        
        variance = var_sum / real(df%n_cells, dp)

        print *, "Mean:", mean, " Variance:", variance

    end subroutine find_mean_variance

    subroutine get_RST(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        integer :: i, j, count, ios
        character(len=500) :: line
        real(kind=dp), allocatable :: urms_us(:), vrms_us(:), wrms_us(:), uvrms_us(:), u_rms(:), v_rms(:), w_rms(:), uv_rms(:)
        real(kind=dp) :: values(28)
        real(kind=dp) :: val, x_est, Re, Cf, tau_w, u_tau, u_Mor

        ! Allocate vectors for arrays
        allocate(urms_us(df%vel_file_N_values))
        allocate(vrms_us(df%vel_file_N_values))
        allocate(wrms_us(df%vel_file_N_values))
        allocate(uvrms_us(df%vel_file_N_values))
        allocate(u_rms(df%vel_file_N_values))
        allocate(v_rms(df%vel_file_N_values))
        allocate(w_rms(df%vel_file_N_values))
        allocate(uv_rms(df%vel_file_N_values))
        allocate(df%R11(df%vel_file_N_values))
        allocate(df%R21(df%vel_file_N_values))
        allocate(df%R22(df%vel_file_N_values))
        allocate(df%R33(df%vel_file_N_values))

        open(unit = 10, file=df%vel_fluc_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening input file'
            stop 
        end if

        ! Skip header lines
        do i = 1, df%vel_file_offset
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
        end do

        count = 1
        
        ! Read lines to gather data
        do  
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit

            if (len_trim(line) == 0) cycle

            read(line, *, iostat=ios) (values(j), j = 1,28)
            if (ios /= 0) cycle

            if (count <= df%vel_file_N_values) then
                urms_us(count) = values(9)
                vrms_us(count) = values(11)
                wrms_us(count) = values(10)
                uvrms_us(count) = values(16)
            end if
            count = count + 1
        end do

        close(10)

        ! Turbulent boundary layer estimates
        Re = df%rho_e * df%U_e / df%mu_e
        val = df%d_i / 0.37_dp * Re**(1.0_dp/5.0_dp)
        x_est = val**(5.0_dp / 4.0_dp)
        print *, "x_est ", x_est

        Cf = 0.0576_dp / (Re * x_est)**(1.0_dp / 5.0_dp)
        tau_w = 0.5_dp * Cf * df%rho_e * df%U_e**2
        u_tau = sqrt(tau_w / df%rhoy(1))

        print *, "Cf ", Cf
        print *, "tau_w ", tau_w
        print *, "u_tau ", u_tau

        ! Scale u'_rms / u* to inflow u*
        do i = 1, df%vel_file_N_values
            u_Mor = sqrt(tau_w / df%rhoy(i))
            u_rms(i) = urms_us(i) * u_Mor
            v_rms(i) = vrms_us(i) * u_Mor
            w_rms(i) = wrms_us(i) * u_Mor
            uv_rms(i) = uvrms_us(i) * u_Mor**2
        end do

        ! Set Reynolds stress terms
        do i = 1, df%vel_file_N_values
            df%R11(i) = u_rms(i)**2
            df%R21(i) = -uv_rms(i)
            df%R22(i) = v_rms(i)**2
            df%R33(i) = w_rms(i)**2
        end do

        df%d_v = 0.0002 * df%d_i    ! Will be changed
        
    end subroutine get_RST


    !  This function writes the velocity fluctuations to a Tecplot file for visualization. It writes the z, y, u_fluc, 
    !  v_fluc, and w_fluc data. @param filename is the name of the file to write to.

    subroutine write_tecplot(df, filename)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        character(len=*), intent(in) :: filename
        integer :: j, k, idx 
        integer :: unit

        inquire(iolength = unit)
        unit = 20

        open(unit=unit, file=filename, status='replace', action='write', iostat=idx)
        if (idx /= 0) then
            print *, 'Error opening file: ', filename
            return
        end if

        write(unit, '(A)') 'VARIABLES = "z", "y", "u_fluc", "v_fluc", "w_fluc"'
        write(unit, '(A,I0,A,I0,A)') 'ZONE T="Flow Field", I=', df%Nz + 1, ', J=', df%Ny + 1, ', F=BLOCK'
        write(unit, '(A)') 'VARLOCATION=([3-5]=CELLCENTERED)'

        ! Write z coordinates
        do j = 1, df%Ny + 1
            do k = 1, df%Nz + 1
                idx = (j - 1) * (df%Nz + 1) + k
                write(unit, '(F12.6)') df%z(idx)
            end do
        end do

        ! Write y coordinates
        do j = 1, df%Ny + 1
            do k = 1, df%Nz + 1
                idx = (j - 1) * (df%Nz + 1) + k
                write(unit, '(F12.6)') df%y(idx)
            end do
        end do

        ! Write u_fluc
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                write(unit, '(F12.6)') df%u_fluc(idx)
            end do
        end do

        ! Write v_fluc
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                write(unit, '(F12.6)') df%v_fluc(idx)
            end do
        end do

        ! Write w_fluc
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                write(unit, '(F12.6)') df%w_fluc(idx)
            end do
        end do

        close(unit)
        print *, 'Finished plotting.'
    end subroutine write_tecplot


    subroutine write_csv(df, filename)
        implicit none
        type(digital_filter_type), intent(in) :: df
        character(len=*), intent(in) :: filename
        integer :: j, k, idx, iidx
        integer :: unit

        ! Assign a free file unit
        unit = 20

        ! Open file for writing
        open(unit=unit, file=filename, status='replace', action='write', iostat=idx)
        if (idx /= 0) then
            print *, 'Error opening file: ', filename
            return
        end if

        ! Write CSV header
        write(unit, '(A)') 'z, y, u_fluc, v_fluc, w_fluc'

        ! Write data row by row
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                iidx = (j - 1) * (df%Nz + 1) + k
                write(unit, '(F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6)') &
                    df%z(iidx), df%y(iidx), df%u_fluc(idx), df%v_fluc(idx), df%w_fluc(idx)
            end do
        end do

        close(unit)
        print *, 'CSV file written to ', filename
    end subroutine write_csv


end module DIGITAL_FILTERING








