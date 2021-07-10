module initial_conditions_mod
!
	use data_structure_mod
!
contains
!
!
		subroutine initial_conditions()
!
			implicit none
!
			integer k, i 
			real*8 U(5),x1,y1,z1
!
		OPEN(UNIT=101,FILE="/home/dhruv/primal-solution.dat",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
!
			do k=1, max_points
				! point%prim(1,k) = rho_inf
				! point%prim(2,k) = u1_inf
				! point%prim(3,k) = u2_inf
				! point%prim(4,k) = u3_inf
				! point%prim(5,k) = pr_inf
				! point%delUp(:,k) = 0.00
				! point%delUn(:,k) = 0.00
!				read(101,*) point%prim(1,k), point%prim(2,k), point%prim(3,k), point%prim(4,k), point%prim(5,k)
				read(101,*) x1,y1,z1,point%prim(1,k), point%prim(2,k), point%prim(3,k), point%prim(4,k), point%prim(5,k)
			enddo
!
		end subroutine
!
!
end module initial_conditions_mod
