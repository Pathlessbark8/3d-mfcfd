module q_lskum_mod
!
!	First written on 05.02.2021
!	
!
	use data_structure_mod
	use generate_connectivity_mod
	use fpi_solver_mod
	use implicit_aliasing_mod
!	
!
contains
!
!
	subroutine q_lskum()
!
			implicit none
!		
			integer t
!
!
			OPEN(UNIT=104,FILE="residue.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")         
!
!
			call generate_split_stencils()
			call aliasing()
!
			do t = 1, max_iters
				call fpi_solver(t)
				print*, t, res_new, residue
				write(104,*) t, res_new, residue
			enddo	
			
!
!
	end subroutine
!
!	
!
end module q_lskum_mod
