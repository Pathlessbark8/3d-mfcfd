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
    ! print *,point%nbhs(1)
    call eval_q_variables()
    !  print *,point%dq(1,:,1)
    call eval_q_derivatives()
    !  print *,point%dq(1,:,1)
    call timestep_delt()
    call eval_flux_residual()
    
   
    ! OPEN(UNIT=80,FILE="flux_res_var_fortran",STATUS="OLD")
    ! do k=1,wall_points
    !     i=wall_points_index(k)
    !     write(80) point%flux_res(:,i)
    ! end do
    ! ! print *,wall_points
    ! close(80)


    ! do i=1,wall_points
    !     k=wall_points_index(i)
    !     print *,k, point%dq(1,:,k)
    !     ! ! print *,k
    !     ! call wall_dGx_neg(G,k)
    !     ! call wall_dGx_pos(G,k)
    !     ! call wall_dGy_neg(G,k)
    !     ! call wall_dGy_pos(G,k)
    !     ! call wall_dGz_neg(G,k)
    !         ! print *, G
    ! end do
    !  do i=1,interior_points
    !     k=interior_points_index(i)
    !     ! ! print *,k
    !     call interior_dGx_neg(G,k)
    !         print *, G
    ! end do
end program test