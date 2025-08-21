! test/main.f90
program test_df
    use DIGITAL_FILTERING   ! refer to the digital filtering module in 'df' directory
    implicit none
    INTEGER,PARAMETER :: dp = selected_real_kind(15)
    real(kind=dp) :: dt

    type(digital_filter_type) :: df ! Create the digital filter
    type(df_config) :: config        ! Create the configuration for input to construtor function

    config%d_i = 0.0013_dp
    config%rho_e = 0.044_dp
    config%U_e = 869.0_dp
    config%mu_e = 1.8e-5_dp
    config%grid_file = 'grid.dat'
    config%vel_fluc_file = "../files/M6Tw025_Stat.dat"
    config%vel_file_offset = 142
    config%vel_file_N_values = 330

    ! Call the constructor
    df = create_digital_filter(config)

    ! Call the filter with a timestep
    dt = 1e-8_dp
    call filter(df, dt)


end program test_df