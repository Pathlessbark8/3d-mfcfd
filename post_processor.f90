module post_processor_mod
!
!	First written on 05.02.2021
!
	
!
	use data_structure_mod
!	
!
contains
!
!
	subroutine print_primal_output()
	!
	!
		implicit none
!			
		integer :: k, r
		
		OPEN(UNIT=101,FILE="primal-solution.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		OPEN(UNIT=102,FILE="normals.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")		
		OPEN(UNIT=103,FILE="flux-residual.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")		
		OPEN(UNIT=104,FILE="qderivatives.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")		
		
!		write(101,*) 'TITLE="Transonic flow over NACA0012 airfoil"'
!		write(101,*) 'Variables="x","y","rho","u1","u2","pr","entropy","vorticity","vorticity_sqr"'
!		write(101,*) 'Zone I = 160, J = 60, F = POINT'
!			do k=1, max_points
!				write(101,*) point(k)%x, point(k)%y, point(k)%rho, point(k)%u1, point(k)%u2, &
!				point(k)%pr, point(k)%entropy, point(k)%vorticity, point(k)%vorticity_sqr
!			enddo
!
		do k = 1, max_points
			write(101,*) point%x(k), point%y(k), point%z(k), point%prim(1,k),point%prim(2,k),point%prim(3,k),point%prim(4,k),point%prim(5,k)
			write(104,*) k, point%status(k),point%dq(3,1,k), point%dq(3,2,k), point%dq(3,3,k), point%dq(3,4,k), point%dq(3,5,k)
				write(103,*) k, point%status(k), (point%flux_res(r,k), r = 1, 5)
!				write(101,*) k, (point%prim(r,k), r = 1, 5)				
!				write(101,*) k, point%delt(k)
		enddo	
!
!
	!				
		CLOSE(UNIT=101)
	!		
	end subroutine		
!
!
!
!
end module post_processor_mod
