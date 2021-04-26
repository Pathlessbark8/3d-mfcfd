module conserved_to_primitive_vector_mod
!
!
	use data_structure_mod
!
!
	contains
!
!
		subroutine conserved_to_primitive(k, U) 
!
			implicit none
!		
			real*8 :: U(5), temp
			real*8 :: U2_rot, U3_rot, U4_rot
			integer :: k
!		
!			
			point%prim(1,k) = U(1)
			temp = 1.0d0/U(1)	
!
!
			U2_rot = U(2)
			U3_rot = U(3)
			U4_rot = U(4)		
!		
!
			U(2) = point%tan1(1,k)*U2_rot + point%tan2(1,k)*U3_rot + point%nor(1,k)*U4_rot
			U(3) = point%tan1(2,k)*U2_rot + point%tan2(2,k)*U3_rot + point%nor(2,k)*U4_rot
			U(4) = point%tan1(3,k)*U2_rot + point%tan2(3,k)*U3_rot + point%nor(3,k)*U4_rot
!
!			
			point%prim(2,k) = temp*U(2)
			point%prim(3,k) = temp*U(3)
			point%prim(4,k) = temp*U(4)
!		
			temp = 	point%prim(2,k)*point%prim(2,k) + point%prim(3,k)*point%prim(3,k) + point%prim(4,k)*point%prim(4,k)

			point%prim(5,k) = 0.4d0*(U(5) - 0.5d0*point%prim(1,k)*temp)
!
!	
		end subroutine
!		
!		
end module conserved_to_primitive_vector_mod
!	
