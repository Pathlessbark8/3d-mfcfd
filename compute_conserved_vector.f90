module compute_conserved_vector_mod
!
!	First written on 10.04.2021.
!	Revised on 17.04.2021.
!
!
	use data_structure_mod
!
!
	contains
!
!
		subroutine compute_conserved_vector()
!
			implicit none
!
			integer :: k
			real*8 :: temp1, temp2, temp3
!
			do k = 1, max_points
				point%U(1,k) = point%prim(1,k)
				temp1 = point%prim(1,k)*point%prim(2,k)
				temp2 = point%prim(1,k)*point%prim(3,k)
				temp3 = point%prim(1,k)*point%prim(4,k)
!
				point%U(2,k) = temp1*point%tan1(1,k) + temp2*point%tan1(2,k) + temp3*point%tan1(3,k)
				point%U(3,k) = temp1*point%tan2(1,k) + temp2*point%tan2(2,k) + temp3*point%tan2(3,k)
				point%U(4,k) = temp1*point%nor(1,k) + temp2*point%nor(2,k) + temp3*point%nor(3,k)
!
				temp1 = point%prim(2,k)*point%prim(2,k) + point%prim(3,k)*point%prim(3,k) + point%prim(4,k)*point%prim(4,k)								
				point%U(5,k) = 2.5d0*point%prim(5,k) + 0.5d0*point%U(1,k)*temp1
			enddo
!
		end subroutine 
!
end module compute_conserved_vector_mod													
