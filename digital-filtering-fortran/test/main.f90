! test/main.f90
program test_df
    use DIGITAL_FILTERING   ! refer to the module in df/
    implicit none
    INTEGER,PARAMETER :: dp = selected_real_kind(15)
    type(digital_filter_type) :: df

    real(kind=dp) :: d_i, rho_e, U_e, mu_e
    character(len=100) :: grid_file, vel_fluc_file
    integer :: offset, N_vals
    offset = 142
    N_vals = 330
    d_i = 0.0013_dp
    U_e = 869.0_dp
    rho_e = 0.044_dp
    mu_e = 1.8e-5_dp
    grid_file = 'grid.dat'
    vel_fluc_file = 'fluc.dat'


    ! Call the constructor
    df = create_digital_filter(d_i, rho_e, U_e, mu_e, grid_file, vel_fluc_file, offset, N_vals)

    print *, 'Constructor test complete.'
    print *, 'Number of cells:', df%n_cells
    print *, 'First y-value:', df%y(1)
end program test_df