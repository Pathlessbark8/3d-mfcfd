!
module fpi_solver_mod
!
!	First written on 06.02.2020
!
!
	use data_structure_mod
	use compute_conserved_vector_mod
	use q_variables_mod
	use q_derivatives_mod
	use timestep_delt_mod
	use flux_residual_mod
	use implicit_forward_and_reverse_sweeps_mod
	use state_update_mod
!
!
contains
!
!
		subroutine fpi_solver(t)
!
			implicit none
!
				integer t,i		
!
!
				call compute_conserved_vector()
				call eval_q_variables()
				call eval_q_derivatives()
!					do i = 1, max_points
!						point%dq(:,:,i)	= 0.00
!					enddo	
				call timestep_delt()
				call eval_flux_residual()
				call implicit_forward_backward_sweeps()
				call state_update()
!
				if(t .le. 2) then
					res_old = res_new
					residue = 0.00
				else 
					residue = dlog10(res_new/res_old)					
				endif	
!
		end subroutine
!
!
end module fpi_solver_mod							
