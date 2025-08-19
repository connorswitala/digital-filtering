! digital filter_module

module DIGITAL_FILTERING
    implicit none
    private
    public :: digital_filter_type, create_digital_filter, filter

    INTEGER,PARAMETER :: dp = selected_real_kind(15)

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

        real(kind=dp), allocatable :: rhow(:), Uy(:), My(:), Ty(:)
        real(kind=dp), allocatable :: ru_ys(:), rv_ys(:), rw_ys(:)
        real(kind=dp), allocatable :: ru_zs(:), rv_zs(:), rw_zs(:)
        integer, allocatable :: Nu_ys(:), Nv_ys(:), Nw_ys(:)
        integer, allocatable :: Nu_zs(:), Nv_zs(:), Nw_zs(:)
        real(kind=dp), allocatable :: bu_y(:), bv_y(:), bw_y(:)
        real(kind=dp), allocatable :: bu_z(:), bv_z(:), bw_z(:)

        real(kind=dp), allocatable :: rho_fluc(:), T_fluc(:)
        real(kind=dp), allocatable :: u_fluc(:), v_fluc(:), w_fluc(:)
        real(kind=dp), allocatable :: u_filt(:), v_filt(:), w_filt(:)
        real(kind=dp), allocatable :: u_filt_old(:),, v_filt_old(:), w_filt_old(:)
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

    end function create_digital_filter


    subroutine read_grid(df)
        type(digital_filter_type), intent(inout) :: df
        integer :: j,kind
        real(kind=dp) :: y_max, eta, a

        ! set sizes
        df%Nz = 450
        df%Ny = df%vel_file_N_values
        df%n_cells = df%Nz * df%Ny

        !allocate arrays
        allocate(df%)


    end subroutine



    subroutine generate_white_noise(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine
    
    subroutine calculate_filter_properties(df)
        type(digital_filter_type), intent(inout) :: df
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

    subroutine rhoy_test(df)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

    subroutine write_tecplot(df, filename)
        type(digital_filter_type), intent(inout) :: df
    end subroutine

end module








