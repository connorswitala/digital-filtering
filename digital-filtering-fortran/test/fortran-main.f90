! test/main.f90
program test_df
    use DIGITAL_FILTERING   ! refer to the module in df/
    implicit none
    INTEGER,PARAMETER :: dp = selected_real_kind(15)
    type(digital_filter_type) :: df

    real(kind=dp) :: d_i, rho_e, U_e, mu_e, dt
    character(len=100) :: grid_file, vel_fluc_file
    integer :: offset, N_vals
    offset = 142
    N_vals = 330
    d_i = 0.0013_dp
    U_e = 869.0_dp
    rho_e = 0.044_dp
    mu_e = 1.8e-5_dp
    dt = 1e-8_dp
    grid_file = 'grid.dat'
    vel_fluc_file = 'fluc.dat'

    vel_fluc_file = "../files/M6Tw025_Stat.dat"


    ! Call the constructor
    df = create_digital_filter(d_i, rho_e, U_e, mu_e, grid_file, vel_fluc_file, offset, N_vals)
    print *, 'Constructor completed.'
    call test(df)
    print *, "Test completed"

end program test_df