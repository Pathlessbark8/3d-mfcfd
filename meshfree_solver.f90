!
!	First written on 02.02.2021
!
!
program meshfree_solver
!
!
	use data_structure_mod
	use point_preprocessor_mod
	use initial_conditions_mod
	use q_lskum_mod
	use post_processor_mod
!
!
		implicit none
!	
!
!	Reading the input data ..
!
		print*, max_points
		call allocate_size()
		call read_input_point_data()
!
		print*, wall_points, outer_points, interior_points
		print*, supersonic_inlet_points, supersonic_outlet_points
!
!	Assign the initial conditions for the primitive variables ..	
!
		call initial_conditions()
!
!	Primal fixed point iterative solver ..
!
		call q_lskum()
!
!	Printing the output (post-processing) ..
!
	call print_primal_output()
!
!	
end program meshfree_solver
