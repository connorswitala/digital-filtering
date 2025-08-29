! digital filter_module

module DIGITAL_FILTERING
    implicit none
    private
    public :: digital_filter_type, create_digital_filter, filter, DFConfig

    INTEGER,PARAMETER :: dp = selected_real_kind(15)
    real(kind=dp), PARAMETER :: pi = acos(-1.0_dp)

    ! Private variables

    type(FilterField) :: u, v, w

    integer :: Ny, Nz, n_cells
    integer :: N_in
    integer :: rms_counter

    real(kind=dp) :: rand1, rand2
    real(kind=dp) :: d_i, d_v
    real(kind=dp) :: rho_e, U_e, mu_e
    real(kind=dp) :: U_w, rho_w, T_w, P, mu, gcon
    character(len=*) :: line_file
    real(kind=dp), allocatable :: Us, rhos, Ts, Ps, Ms
    real(kind=dp) :: u_tau, tau_w

    real(kind=dp), allocatable :: rho_fluc(:), T_fluc(:)
    real(kind=dp), allocatable :: rhoy(:), Uy(:), My(:), Ty(:)
    real(kind=dp), allocatable :: R11(:), R21(:), R22(:), R33(:)
    real(kind=dp), allocatable :: R11_in(:), R21_in(:), R22_in(:), R33_in(:)

    real(kind=dp) :: dt

    real(kind=dp), allocatable :: y(:), & ! y vertices for entire inflow size (Ny * Nz)
    yc(:),      & ! y cell-center for entire inflow
    z(:),       & ! cell heights for entire inflow
    dy(:),      & ! cell widths for entire inflow
    dz(:),      & ! y cell centers normalized by d_i
    yc_d(:),    & ! y locations normalized by delta (from Duan data)
    yin_d(:),   & ! y locations (from Duan data)
    ydline(:),  & ! Assuming y stretching is constant across z, this is all y locations normalized by our delta for DF
    yline(:),   & ! Same as above but not normalized;


    type :: FilterField

        real(kind=dp), allocatable :: fluc(:), filt(:), filt_old(:)
        real(kind=dp), allocatable :: r_ys(:), r_zs(:)
        real(kind=dp), allocatable :: by(:), bz(:)
        real(kind=dp), allocatable :: rms(:), rms_added(:)
        integer, allocatable :: N_ys(:), N_zs(:)
        integer, allocatable :: bz_offsets(:), by_offsets(:)

        integer :: Ny_max, Nz_max
        real(kind=dp) :: Iz_inn, Iz_out, Lt

    end type FilterField

    type :: DFConfig
        real(kind=dp) :: d_i, rho_e, U_e, mu_e
        integer :: vel_file_offset, vel_file_N_values
        character(len=256) :: grid_file, vel_fluc_file
    end type DFConfig

    type :: geometry
        real(kind=dp), allocatable :: x(:), y(:), z(:)
        real(kind=dp), allocatable :: xc(:), yc(:), zc(:)
        real(kind=dp), allocatable :: dx(:), dy(:), dz(:)
    end type geometry
contains

    !   Constructor subroutine

    function create_digital_filter(config) result(DF)
        implicit none
        type(DFConfig), intent(in) :: config
        character(len=50) :: filename
        type(digital_filter_type) :: DF
        
        DF%d_i = config%d_i
        DF%rho_e = config%rho_e
        DF%U_e = config%U_e
        DF%mu_e = config%mu_e
        DF%grid_file = config%grid_file
        DF%vel_fluc_file = config%vel_fluc_file
        DF%vel_file_offset = config%vel_file_offset
        DF%vel_file_N_values = config%vel_file_N_values

        DF%rms_counter = 0

        call read_grid(DF)
        call rhoy_test(DF)
        call get_RST_in(DF)
        call lerp_RST(DF)
        call plot_RST_lerp(DF)

        ! Allocate data structures for u, v, w
        call allocate_data_structures(DF, DF%u)
        call allocate_data_structures(DF, DF%v)
        call allocate_data_structures(DF, DF%w)

        ! Set integral length scales
        DF%u%Iz_out = 0.4_dp * DF%d_i
        DF%u%Iz_inn = 150.0_dp * DF%d_v
        DF%u%Lt = 0.8_dp * DF%d_i / DF%U_e

        DF%v%Iz_out = 0.3_dp * DF%d_i
        DF%v%Iz_inn = 75.0_dp * DF%d_v
        DF%v%Lt = 0.3_dp * DF%d_i / DF%U_e

        DF%w%Iz_out = 0.4_dp * DF%d_i
        DF%w%Iz_inn = 150.0_dp * DF%d_v
        DF%w%Lt = 0.3_dp * DF%d_i / DF%U_e


        allocate(DF%My(DF%Ny))
        allocate(DF%Uy(DF%Ny))
        allocate(DF%Ty(DF%Ny))
        allocate(DF%rho_fluc(DF%n_cells))
        allocate(DF%T_fluc(DF%n_cells))
     
        call calculate_filter_properties(DF, DF%u)
        call calculate_filter_properties(DF, DF%v)
        call calculate_filter_properties(DF, DF%w)

        ! First timestep
        call generate_white_noise(DF)

        call filtering_sweeps(DF, DF%u)
        call filtering_sweeps(DF, DF%v)
        call filtering_sweeps(DF, DF%w)

        call apply_RST_scaling(DF)

        filename = "../files/f_vel_fluc_init.dat"
        call write_tecplot(DF, filename)

    end function create_digital_filter


    ! ========: Subroutines for Digital Filter :========

    subroutine read_grid(DF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        real(kind=dp) :: y_max, eta, a
        integer :: j, jj, k, idx

        ! set sizes
        DF%Nz = 400
        DF%Ny = DF%vel_file_N_values
        DF%n_cells = DF%Nz * DF%Ny

        !allocate arrays
        allocate(DF%y( (DF%Ny + 1) * (DF%Nz + 1)))
        allocate(DF%z( (DF%Ny + 1) * (DF%Nz + 1)))
        allocate(DF%yc( (DF%n_cells)))
        allocate(DF%dy( (DF%n_cells)))
        allocate(DF%dz( (DF%n_cells)))

        y_max = 2 * DF%d_i
        eta = 0.0_dp
        a = 2.0_dp
        

        !--- First loop: y and z
        do j = DF%Ny + 1, 1, -1        ! j decreases: Ny+1 → 1
            jj = DF%Ny + 2 - j          ! jj increases: 1 → Ny+1
            do k = 1, DF%Nz + 1
                eta = real(j - 1, dp) / real(DF%Ny + 1, dp)

                DF%y((jj - 1) * (DF%Nz + 1) + k) = y_max * (1.0_dp - tanh(a * eta) / tanh(a))
                DF%z((jj - 1) * (DF%Nz + 1) + k) = real(k,dp) * 0.000133_dp
            end do
        end do


        !--- Second loop: dy, dz, yc
        do j = 1, DF%Ny
            do k = 1, DF%Nz
                idx = (j - 1) * DF%Nz + k
                DF%dy(idx) = DF%y((j + 0) * (DF%Nz + 1) + k) - DF%y((j - 1) * (DF%Nz + 1) + k)
                DF%dz(idx) = 0.000133_dp
                DF%yc(idx) = 0.25_dp * ( DF%y((j - 1)*(DF%Nz + 1) + k) + DF%y(j*(DF%Nz + 1) + k) &
                                    + DF%y((j - 1)*(DF%Nz + 1) + k + 1) + DF%y(j*(DF%Nz + 1) + k + 1) )
            end do
        end do
    end subroutine read_grid

    subroutine allocate_data_structures(DF, FF) 
            implicit none
            type(digital_filter_type), intent(inout) :: DF
            type(FilterField), intent(inout) :: FF

            allocate(FF%fluc(DF%n_cells))
            allocate(FF%filt(DF%n_cells))
            allocate(FF%filt_old(DF%n_cells))
            allocate(FF%N_ys(DF%n_cells))
            allocate(FF%N_zs(DF%n_cells))
            allocate(FF%by_offsets(DF%n_cells))
            allocate(FF%bz_offsets(DF%n_cells))

    end subroutine allocate_data_structures

    subroutine calculate_filter_properties(DF, FF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        type(FilterField), intent(inout) :: FF
        integer :: i, idx, N, n_val, b_size, offset
        real(kind=dp) :: n_int, sum, val, Iy 
        real(kind=dp), allocatable :: Iz(:)

        allocate(Iz(DF%n_cells))

        !===================================================================================================
        ! Find filter half-width and convolution coefficients for velocity when filtering in the z-direction
        !===================================================================================================

        FF%Ny_max = 0
        FF%Nz_max = 0
        b_size = 0

        do idx = 1, DF%n_cells

            Iz(idx) = FF%Iz_inn + (FF%Iz_out - FF%Iz_inn) * 0.5_dp * (1.0_dp + tanh((DF%yc(idx) / DF%d_i - 0.2_dp) / 0.03_dp))
            n_int = max(1.0_dp, Iz(idx) / DF%dz(idx))
            n_val = 2 * int(n_int)
            FF%N_zs(idx) = n_val

            b_size = b_size + 2 * n_val + 1
            FF%bz_offsets(idx) = b_size - n_val
            if (n_val > FF%Nz_max) FF%Nz_max = n_val
        end do

        allocate(FF%r_zs((DF%Nz + 2 * FF%Nz_max) * DF%Ny))
        allocate(FF%bz(b_size))

        do idx = 1, DF%n_cells

            N = FF%N_zs(idx)
            offset = FF%bz_offsets(idx)

            sum = 0.0_dp

            do i = -N, N
                val = exp(-2.0_dp * pi * real(abs(i), dp) / real(N, dp) )
                sum = sum + val * val
            end do
            sum = sqrt(sum)

            do i = -N, N
                FF%bz(offset + i) = exp(-2.0_dp * pi * abs(i) / real(N, dp)) / sum
            end do
        end do

        !===================================================================================================
        ! Find filter half-width and convolution coefficients for velocity when filtering in the y-direction
        !===================================================================================================

        b_size = 0

        do idx = 1, DF%n_cells

            Iy = 0.67_dp * Iz(idx)
            n_int = max(1.0_dp, Iy / DF%dy(idx))
            n_val = 2 * int(n_int)

            b_size = b_size + 2 * n_val + 1
            FF%by_offsets(idx) = b_size - n_val
            FF%N_ys(idx) = n_val            
            if (n_val > FF%Ny_max) FF%Ny_max = n_val
        end do

        allocate(FF%r_ys(DF%Nz * (2 * FF%Ny_max + DF%Ny)))
        allocate(FF%by(b_size))

        do idx = 1, DF%n_cells

            N = FF%N_ys(idx)
            offset = FF%by_offsets(idx)

            sum = 0.0_dp

            do i = -N, N
                val = exp(-2.0_dp * pi * real(abs(i), dp) / real(N, dp) )
                sum = sum + val * val
            end do
            sum = sqrt(sum)

            do i = -N, N
                FF%by(offset + i) = exp(-2.0_dp * pi * abs(i) / real(N, dp)) / sum
            end do
        end do

    end subroutine calculate_filter_properties

    subroutine get_RST_in(DF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        integer :: i, j, count, ios
        character(len=500) :: line
        real(kind=dp), allocatable :: urms_us(:), vrms_us(:), wrms_us(:), uvrms_us(:), u_rms(:), v_rms(:), w_rms(:), uv_rms(:)
        real(kind=dp) :: values(28)
        real(kind=dp) :: val, x_est, Re, Cf, tau_w, u_tau

        ! Allocate vectors for arrays
        allocate(DF%y_d(DF%vel_file_N_values))
        allocate(urms_us(DF%vel_file_N_values))
        allocate(vrms_us(DF%vel_file_N_values))
        allocate(wrms_us(DF%vel_file_N_values))
        allocate(uvrms_us(DF%vel_file_N_values))
        allocate(u_rms(DF%vel_file_N_values))
        allocate(v_rms(DF%vel_file_N_values))
        allocate(w_rms(DF%vel_file_N_values))
        allocate(uv_rms(DF%vel_file_N_values))
        allocate(DF%R11(DF%vel_file_N_values))
        allocate(DF%R21(DF%vel_file_N_values))
        allocate(DF%R22(DF%vel_file_N_values))
        allocate(DF%R33(DF%vel_file_N_values))
        allocate(DF%R11_in(DF%vel_file_N_values))
        allocate(DF%R21_in(DF%vel_file_N_values))
        allocate(DF%R22_in(DF%vel_file_N_values))
        allocate(DF%R33_in(DF%vel_file_N_values))


        open(unit = 10, file=DF%vel_fluc_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening input file'
            stop 
        end if

        ! Skip header lines
        do i = 1, DF%vel_file_offset
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

            if (count <= DF%vel_file_N_values) then
                DF%y_d(count) = values(2)
                urms_us(count) = values(9)
                vrms_us(count) = values(11)
                wrms_us(count) = values(10)
                uvrms_us(count) = values(16)
            end if
            count = count + 1
        end do

        close(10)

        ! Turbulent boundary layer estimates
        Re = DF%rho_e * DF%U_e / DF%mu_e
        val = DF%d_i / 0.37_dp * Re**(1.0_dp/5.0_dp)
        x_est = val**(5.0_dp / 4.0_dp)

        Cf = 0.0576_dp / (Re * x_est)**(1.0_dp / 5.0_dp)
        tau_w = 33.6_dp
        u_tau = sqrt(tau_w / 0.0264_dp)

        ! Scale u'_rms / u* to inflow u*
        do i = 1, DF%vel_file_N_values
            u_rms(i) = urms_us(i) * u_tau
            v_rms(i) = vrms_us(i) * u_tau
            w_rms(i) = wrms_us(i) * u_tau
            uv_rms(i) = uvrms_us(i) * u_tau**2
        end do

        ! Set Reynolds stress terms
        do i = 1, DF%vel_file_N_values
            DF%R11_in(i) = u_rms(i)**2
            DF%R21_in(i) = uv_rms(i)
            DF%R22_in(i) = v_rms(i)**2
            DF%R33_in(i) = w_rms(i)**2
        end do

        DF%d_v = 0.0002 * DF%d_i    ! Will be changed
        
    end subroutine get_RST_in

    subroutine lerp_RST(DF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        integer :: j, idx, i0, i1
        real(kind=dp) :: yq, x0, x1, dx, t

        do j = 1, DF%Ny
            idx = (j - 1) * DF%Nz + 1

            yq = DF%yc(idx) / DF%d_i

            if (yq <= DF%y_d(1)) then
                DF%R11(j) = DF%R11_in(1)
                DF%R21(j) = DF%R21_in(1)
                DF%R22(j) = DF%R22_in(1)
                DF%R33(j) = DF%R33_in(1)
                cycle
            end if

            if (yq >= DF%y_d(DF%Ny)) then
                DF%R11(j) = DF%R11_in(DF%Ny)
                DF%R21(j) = DF%R21_in(DF%Ny)
                DF%R22(j) = DF%R22_in(DF%Ny)
                DF%R33(j) = DF%R33_in(DF%Ny)
                cycle
            end if

            i1 = lower_bound_bin(DF%y_d, yq)   


            if (i1 <= 1) then
                i0 = 1
                i1 = 2
            else if (i1 > DF%Ny) then
                i1 = DF%Ny
                i0 = DF%Ny-1
            else
                i0 = i1 - 1
            end if

            x0 = DF%y_d(i0)
            x1 = DF%y_d(i1)
            dx = x1 - x0
            if (dx == 0.0_dp) then
                t = 0.0_dp
            else
                t = (yq - x0) / dx
            end if

            DF%R11(j) = (1.0_dp - t) * DF%R11_in(i0) + t * DF%R11_in(i1)
            DF%R21(j) = (1.0_dp - t) * DF%R21_in(i0) + t * DF%R21_in(i1)
            DF%R22(j) = (1.0_dp - t) * DF%R22_in(i0) + t * DF%R22_in(i1)
            DF%R33(j) = (1.0_dp - t) * DF%R33_in(i0) + t * DF%R33_in(i1)

        end do     

    end subroutine lerp_RST 

    function lower_bound_bin(arr, value) result(idx)
        implicit none
        real(kind=dp), intent(in) :: arr(:)     ! sorted ascending
        real(kind=dp), intent(in) :: value
        integer :: left, right, mid, idx

        left  = 1
        right = size(arr)
        idx   = right + 1       ! default if value > all elements

        do while (left <= right)
            mid = (left + right) / 2
            if (arr(mid) >= value) then
                idx   = mid
                right = mid - 1
            else
                left  = mid + 1
            end if
        end do
    end function lower_bound_bin

    subroutine generate_white_noise(DF)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        type(digital_filter_type), intent(inout) :: DF

        call fill_with_noise(DF%u%r_ys)
        call fill_with_noise(DF%u%r_zs)
        call fill_with_noise(DF%v%r_ys)
        call fill_with_noise(DF%v%r_zs)
        call fill_with_noise(DF%w%r_ys)
        call fill_with_noise(DF%w%r_zs)
 
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
    
    subroutine filtering_sweeps(DF, FF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        type(FilterField), intent(inout) :: FF
        integer :: r_idy, r_idz, idx, N, offset, Ny_pad, Nz_pad, i, j, k
        real(kind=dp) :: sum

        ! ========================================
        ! ------: Filter in y-direction :------
        ! ========================================
        
        Nz_pad = FF%Nz_max
        Ny_pad = FF%Ny_max

        do j = 1, DF%Ny

            r_idy = (j - 1 + Ny_pad) * DF%Nz
            r_idz = (j - 1) * (DF%Nz + 2 * Nz_pad) + Nz_pad
            idx = (j - 1) * DF%Nz

            do k = 1, DF%Nz

                r_idy = r_idy + 1
                r_idz = r_idz + 1
                idx = idx + 1

                N = FF%N_ys(idx)
                offset = FF%by_offsets(idx)

                sum = 0.0
                do i = -N, N
                    sum = sum + FF%by(offset + i) * FF%r_ys(r_idy + i * DF%Nz)
                end do

                FF%r_zs(r_idz) = sum

            end do
        end do


        ! ========================================
        ! ------: Filter in z-direction :------
        ! ========================================

        do j = 1, DF%Ny

            r_idz = (j - 1) * (DF%Nz + 2 * Nz_pad) + Nz_pad
            idx = (j - 1) * DF%Nz

            do k = 1, DF%Nz

                r_idz = r_idz + 1
                idx = idx + 1

                N = FF%N_zs(idx)
                offset = FF%bz_offsets(idx)

                sum = 0.0
                do i = -N, N
                    sum = sum + FF%bz(offset + i) * FF%r_zs(r_idz + i)
                end do

                FF%filt(idx) = sum

            end do
        end do

    end subroutine filtering_sweeps

    subroutine correlate_fields(DF, FF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        type(FilterField), intent(inout) :: FF
        integer :: idx


        do idx = 1, DF%n_cells
                FF%filt(idx) = FF%filt_old(idx) * exp(-pi * DF%dt / (2.0_dp * FF%Lt)) &
                                + FF%filt(idx) * sqrt(1.0_dp - exp(-pi * DF%dt / FF%Lt))
        end do

    end subroutine correlate_fields

    subroutine apply_RST_scaling(DF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        integer :: j, k, idx
        real(kind=dp) :: b

        do j = 1, DF%Ny

            idx = (j - 1) * DF%Nz
            if (DF%R11(j) < 1.0e-10_dp) then
                b = 0.0_dp
            else 
                b = DF%R21(j) / sqrt(DF%R11(j))
            end if

            do k = 1, DF%Nz
                idx = idx + 1

                DF%u%fluc(idx) = sqrt(DF%R11(j)) * DF%u%filt(idx)
                DF%v%fluc(idx) = b * DF%u%filt(idx) + sqrt(DF%R22(j) - b**2) * DF%v%filt(idx)
                DF%w%fluc(idx) = sqrt(DF%R33(j)) * DF%w%filt(idx)

                DF%u%filt_old(idx) = DF%u%filt(idx)
                DF%v%filt_old(idx) = DF%v%filt(idx)
                DF%w%filt_old(idx) = DF%w%filt(idx)

            end do
        end do

    end subroutine apply_RST_scaling

    subroutine filter(DF, dt_input)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        real(kind=dp), intent(in) :: dt_input
        real(kind=dp) :: start, end, elapsed
        character(len=50) :: filename

        DF%dt = dt_input

        call cpu_time(start)

        call generate_white_noise(DF)  
 
        call filtering_sweeps(DF, DF%u)
        call filtering_sweeps(DF, DF%v)
        call filtering_sweeps(DF, DF%w)

        call correlate_fields(DF, DF%u)
        call correlate_fields(DF, DF%v)
        call correlate_fields(DF, DF%w)

        call apply_RST_scaling(DF)

        call cpu_time(end)
        elapsed = end - start
        print *, "Time for filtering step (s): ", elapsed

        filename = "../files/f_vel_fluc.csv"
        call write_csv(DF, filename)

    end subroutine filter

    subroutine get_rho_T_fluc(DF)
        type(digital_filter_type), intent(inout) :: DF
        real(kind=dp) :: val
        integer :: j, k, idx

        do j = 1, DF%Ny
            do k = 1, DF%Nz
                idx = (j - 1) * DF%Nz + k
                val = (-1.4_dp - 1.0_dp) * DF%My(j)**2 * DF%Uy(j) * DF%u%fluc(idx)
                DF%T_fluc(idx) = val
                DF%rho_fluc(idx) = -val / DF%Ty(j) * DF%rhoy(j)
            end do
        end do
        
    end subroutine get_rho_T_fluc


    ! ========: Debugging Subroutines :========

    ! subroutine test(DF)
    !     implicit none
    !     type(digital_filter_type), intent(inout) :: DF
    !     character(len=50) :: filename
    !     filename = "../../files/fluc.csv"
    !     call generate_white_noise(DF)
    !     call filtering_sweeps(DF)
    !     call correlate_fields_ts1(DF)
    !     call write_csv(DF, filename)    
    ! end subroutine test

    subroutine display_data(array)
            implicit none
            real(kind=dp), intent(in) :: array(:)       ! the array to display
            integer :: i

            do i = 1, size(array)
                print *, array(i)
            end do
    end subroutine display_data

    subroutine find_mean_variance(DF, v)
        implicit none
        type(digital_filter_type), intent(inout) :: DF

        real(kind=dp), intent(in) :: v(DF%n_cells)
        integer i  
        real(kind=dp) :: mean_sum, var_sum, mean, variance
    
        mean_sum = 0.0_dp
        var_sum = 0.0_dp

        do i = 1, DF%n_cells
            mean_sum = mean_sum + v(i)
        end do

        mean = mean_sum / real(DF%n_cells, dp)

        do i = 1, DF%n_cells
            var_sum = var_sum + (v(i) - mean)**2
        end do
        
        variance = var_sum / real(DF%n_cells, dp)

        print *, "Mean:", mean, " Variance:", variance

    end subroutine find_mean_variance

    subroutine rhoy_test(DF)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        integer j

        allocate(DF%rhoy(DF%Ny))       

        do j = 1, DF%Ny
            DF%rhoy(j) = DF%rho_e * (0.7_dp * DF%yc((j-1) * DF%Nz) / DF%d_i + 0.6_dp)
        end do
    end subroutine rhoy_test


    ! ========: Plotting Subroutines :========
    
    subroutine write_tecplot(DF, filename)
        implicit none
        type(digital_filter_type), intent(inout) :: DF
        character(len=*), intent(in) :: filename
        integer :: j, k, idx 
        integer :: unit, ios

        unit = 20
        open(unit=unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', filename
            return
        end if

        write(unit, '(A)') 'VARIABLES = "z", "y", "u_fluc", "v_fluc", "w_fluc"'
        write(unit, '(A,I0,A,I0,A)') 'ZONE T="Flow Field", I=', DF%Nz + 1, ', J=', DF%Ny + 1, ', F=BLOCK'
        write(unit, '(A)') 'VARLOCATION=([3-5]=CELLCENTERED)'

        ! z coordinates (node-centered)
        do j = 1, DF%Ny + 1
            do k = 1, DF% Nz + 1
                idx = (j - 1) * (DF%Nz + 1) + k
                write(unit,'(F12.6)') DF%z(idx)
            end do
        end do

        ! y coordinates (node-centered)
        do j = 1, DF%Ny + 1
            do k = 1, DF%Nz + 1
                idx = (j - 1) * (DF%Nz + 1) + k
                write(unit,'(F12.6)') DF%y(idx)
            end do
        end do

        ! u, v, w fluctuations (cell-centered)
        do j = 1, DF%Ny
            do k = 1, DF% Nz
                idx = (j - 1) * DF%Nz + k
                write(unit,'(F12.6)') DF%u%fluc(idx)
            end do
        end do

        do j = 1, DF%Ny
            do k = 1, DF% Nz
                idx = (j - 1) * DF%Nz + k
                write(unit,'(F12.6)') DF%v%fluc(idx)
            end do
        end do

        do j = 1, DF%Ny
            do k = 1, DF% Nz
                idx = (j - 1) * DF%Nz + k
                write(unit,'(F12.6)') DF%w%fluc(idx)
            end do
        end do

        close(unit)
        print *, 'Finished plotting: ', filename
    end subroutine write_tecplot

    subroutine write_csv(DF, filename)
        implicit none
        type(digital_filter_type), intent(in) :: DF
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
        do j = 1, DF%Ny
            do k = 1, DF%Nz
                idx = (j - 1) * DF%Nz + k
                iidx = (j - 1) * (DF%Nz + 1) + k
                write(unit, '(F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6)') &
                    DF%z(iidx), DF%y(iidx), DF%u%fluc(idx), DF%v%fluc(idx), DF%w%fluc(idx)
            end do
        end do

        close(unit)
        print *, 'CSV file written to ', filename
    end subroutine write_csv

    subroutine plot_RST_lerp(DF)
        implicit none
        type(digital_filter_type), intent(in) :: DF
        character(len=50) :: filename
        integer :: j, idx
        integer :: unit


        filename = "../files/f_RST.csv"
        ! Assign a free file unit
        unit = 20

        ! Open file for writing
        open(unit=unit, file=filename, status='replace', action='write', iostat=idx)
        if (idx /= 0) then
            print *, 'Error opening file: ', filename
            return
        end if

        ! Write CSV header
        write(unit, '(A)') 'y_d, y, R11_in, R11, R22_in, R22, R33_in, R33, R21_in, R21'

        ! Write data row by row
        do j = 1, DF%Ny
            write(unit, '(F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6)') &
                DF%y_d(j) * 0.0036, DF%y(j * (DF%Nz + 1)), DF%R11_in(j), DF%R11(j), DF%R22_in(j), DF%R22(j), DF%R33_in(j), DF%R33(j), DF%R21_in(j), DF%R21(j)
        end do

        close(unit)
        print *, 'CSV file written to ', filename
    end subroutine plot_RST_lerp

end module DIGITAL_FILTERING