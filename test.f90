program test
    use split_fluxes_mod
    use octant_fluxes_mod
    use wall_flux_dGxneg_mod
    use wall_flux_dGyneg_mod
    use wall_flux_dGxpos_mod
    use wall_flux_dGypos_mod
    use wall_flux_dGzneg_mod
    use q_derivatives_mod
    use q_variables_mod
    use point_preprocessor_mod
    use compute_conserved_vector_mod
    use timestep_delt_mod
    use implicit_aliasing_mod
    use generate_connectivity_mod
    use initial_conditions_mod
    use interior_flux_dGxneg_mod
    use flux_residual_mod
    implicit none

    double precision G(5)
    integer i,k
    ! open(1,file='test_data.dat',status='old')
    ! open(2,file='fortran_out.dat',status='old')
    G=0
    call allocate_size()
    call read_input_point_data()
    call initial_conditions()
    call generate_split_stencils()
    call aliasing()
    call compute_conserved_vector()
    call eval_q_variables()
    call eval_q_derivatives()
    call timestep_delt()
    call eval_flux_residual()
end program test