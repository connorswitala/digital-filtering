! digital filter_module

module DIGITAL_FILTERING
    implicit none
    private
    public :: digital_filter_type, create_digital_filter, filter

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

contains

    ! Constructor subroutine
    function create_digital_filter(d_i, rho_e, U_e, mu_e, grid_file, vel_fluc_file, vel_file_offset, vel_file_N_values) result(obj)
        implicit none
        real, intent(in) :: d_i, rho_e, U_e 
        character(len=*), intent(in) :: grid_file, vel_fluc_file
        integer, intent(in) :: vel_file_offset, vel_file_N_values
        type(digital_filter_type) :: df

        df%d_i = d_i
        df%rho_e = rho_e
        df%U_e = U_e
        df%mu_e = mu_e
        df%grid_file = grid_file
        df%vel_fluc_file = vel_fluc_file
        df%vel_file_offset = vel_file_offset
        df%vel_file_N_values = vel_file_N_values

        call read_grid(df)
        call rhoy_test(df)

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


    end function create_digital_filter


!   Somehow, this function reads the grid from a file and populates the y, yc, z, dy, dz vectors. 
!   I do not know how to go about this part, but I do know that we need to get the vectors y, yc, z,
!   dy, dz and intgeres Nz, Ny populated with the data from the file. 
 
    subroutine read_grid(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        real(kind=dp) :: y_max, eta, a

        ! set sizes
        df%Nz = 450
        df%Ny = df%vel_file_N_values
        df%n_cells = df%Nz * df%Ny

        !allocate arrays
        allocate(df%y( (df%Ny + 1) * (df%Nz + 1)))
        allocate(df%z( (df%Ny + 1) * (df%Nz + 1)))
        allocate(df%yc( (df%n_cells)))
        allocate(df%dy( (df%n_cells)))
        allocate(df%dz( (df%n_cells)))

        real(kind=dp) :: y_max = 0.01
        real(kind=dp) :: eta = 0.0
        real(kind=dp) :: a = 3.0
        
        integer :: j, k, idx

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

!   This function is used to test the rho_y vector.

    subroutine rhoy_test(df)
        implicit none
        type(digital_filter_type), intent(inout) :: df
        allocate(df%rhoy(df%Ny))

        do j = 1, dy%Ny
            rhoy(j) = df%rho_e * (0.7 * df%yc(j * df%Nz) / df%d_i + 0.5)
        end do
    end subroutine rhoy_test

 !  This function generates white noise for the filtering. It uses the PCG random number generator 
 !  to generate normally distributed random numbers.

    subroutine generate_white_noise(df)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        type(digital_filter_type), intent(inout) :: df
 
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

        call fill_with_noise(df%ru_ys)
        call fill_with_noise(df%ru_zs)
        call fill_with_noise(df%rv_ys)
        call fill_with_noise(df%rv_zs)
        call fill_with_noise(df%rw_ys)
        call fill_with_noise(df%rw_zs)

    end subroutine generate_white_noise
    
    subroutine calculate_filter_properties(df)
        type(digital_filter_type), intent(inout) :: df
        implicit none

        integer Ny_max, Nz_max, j, k, idx, kk
        real(kind=dp) n_int, val


        !=============================================================================================
        ! Find filter half-width and convolution coefficients for u' when filtering in the z-direction
        !=============================================================================================

        ! Define integral length scales
        real(kind=dp) Iz_out = 0.4_dp * df%d_i
        real(kind=dp) Iz_inn = 150_dp * df%d_v

        ! Holdering for maximum filter widths. Used for creating ghost cells.
        integer Ny_max = 0, Nz_max = 0 

        !! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, Ny
            do k = 1, Nz

                idx = (j - 1) * Nz + k
                df%Iz(idx) = Iz_inn + (Iz_out - Iz_inn) * 0.5_dp * (1.0_dp + tanh((df%yc(idx) / df%d_i - 0.2_dp) / 0.03_dp))

                n_int = max(1.0_dp, df%Iz(idx) / df%dx(idx))
                df%Nz_us(idx) = 2 * int(n_int)

                if (n_int > Nz_max) Nz_max = 2 * int(n_int)

            end do  
        end do

        ! allocate space for random data and filter coefficient arrays
        N_holder(1) = Nz_max
        allocate(df%ru_zs(df%Nz + 2 * Nz_max) * df%Ny) 
        allocate(df%bu_z(2 * Nz_max + 1))              

        ! Loop that calculates sum square of filter coefficients for u' in z-direction
        real(kind=dp) sum = 0.0_dp
        do i = 1, 2 * Nz_max + 1
            kk = i - Nz_max
            val = exp(-2.0_dp * pi * abs(kk) / Nz_max)
            sum = sum + val * val
        end do
        sum = sqrt(sum)

        ! Loop that calculates individual filter coefficients
        do i = 1, 2 * Nz_max + 1
            kk = i - Nz_max
            df%bu_z(i) = exp(-2 * pi * abs(kk) / Nz_max) / sum
        end do

        !=============================================================================================
        ! Find filter half-width and convolution coefficients for u' when filtering in the y-direction
        !=============================================================================================

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, Ny
            for k = 1, Nz
                idx = (j - 1) * Nz + k
                df%Iy(idx) = 0.67_dp * df%Iz(idx)
                n_int = max(1.0_dp, df%(Iy(idx)) / df%dy(idx))
                df%Nu_ys(idx) = 2 * int(n_int)
                if (n_int > Ny_max) Ny_max = 2 * int(n_int)
            end do
        end do

        ! allocate space for random data and filter coefficient arrays
        df%N_holder(2) = Ny_max
        allocate(df%ru_ys(df%Nz * (2 * Ny_max + df%Ny)))
        allocate(df%bu_y)(2 * Ny_max + 1)


        ! Loop that calculates sum square of filter coefficients
        sum = 0.0_dp
        do i = 1, 2 * Ny_max + 1
            kk = i - Ny_max
            val = exp(-2.0_dp * pi * abs(kk) / Ny_max)
            sum = sum + val * val
        end do
        sum = sqrt(sum)

        ! Loop that calculates individual filter coefficients
        do i = 1, 2 * Ny_max + 1
            kk = i - Ny_max
            df%bu_y(i) = exp(-2.0_dp * pi * abs(kk) / Ny_max) / sum
        end do  

        !=============================================================================================
        ! Find filter half-width and convolution coefficients for v' when filtering in the z-direction
        !=============================================================================================

        ! Reset max filter widths to 0
        Nz_max = 0
        Ny_max = 0

        ! Define integral length scales
        Iz_out = 0.3_dp * df%d_i
        Iz_inn = 75.0_dp * df%d_v

        ! Loop that calculates integral length scales per cell and filter half-widths
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k 
                df%Iz(idx) = Iz_inn + (Iz_out - Iz_inn) * 0.5_dp * (1.0_dp + tanh((df%yc(idx) / df%d_i - 0.2_dp) / 0.03_dp))
                n_int = max(1.0_dp, df%Iz(idx) / df%dz(idx))
                df%Nv_zs(idx) = 2 * int(n_int)
                if (n_int > Nz_max) Nz_max = 2 * int(n_int)
            end do  
        end do

        ! allocate space for random data and filter coefficient arrays
        df%N_holder(3) = Nz_max
        allocate(df%rv_zs(df%Nz + 2 * Nz_max) * df%Ny)
        allocate(df%bv_z(2 * Nz_max + 1))

        ! Loop that calculates sum square of filter coefficients
        sum = 0.0
        do i = 1, 2 * Nz_max + 1
            kk = i - Nz_max
            val = exp(-2.0_dp * abs(kk) / Nz_max)
            sum = sum + val *  val
        end do
        sum = sqrt(sum)

        ! Loop that calculates individual filter coefficients
        do i = 1, 2 * Nz_max + 1
            kk = i - Nz_max
            df%bv_z(i) = exp(-2.0_dp * pi * abs(kk) / Nz_max) / sum
        end do

        !=============================================================================================
        ! Find filter half-width and convolution coefficients for v' when filtering in the y-direction
        !=============================================================================================

        ! Loop to find integral length scales per cell and get half-widths. 
        do j = 1, df%Ny
            do k = 1, df%Nz
                idx = (j - 1) * df%Nz + k
                df%Iy(idx) = 0.67 * df%Iz(idx)
                n_int = max(1.0_dp, df%Iy(idx) / df%dy(idx))
                df%Nv_ys(idx) = 2 * int(n_int)
                if (n_int > Ny_max) Ny_max = 2 * int(n_int)
            end do
        end do

        ! allocate space for random data and filter coefficient arrays
        df%N_holder(4) = Ny_max
        allocate(df%rv_ys(df%Nz * (2 * Ny_max + df%Ny)))
        allocate(df%bv_y(2 * Ny_max + 1))

        ! Loop that calculates sum square of filter coefficients
        sum = 0.0_dp



        ! Reset max filter widths to 0
        Nz_max = 0
        Ny_max = 0


                    

    end subroutine

    subroutine correlate_fields(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine correlate_fields_ts1(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine filtering_sweeps(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine filter(df, dt_input)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine get_rho_T_fluc(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine set_old(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine test(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine display_data(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine find_mean_variance(df, v)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine get_RST(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine write_tecplot(df, filename)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

end module








