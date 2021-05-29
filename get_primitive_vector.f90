module get_primitive_vector_mod
!
!
	use data_structure_mod
!
!
	contains
!
!
		subroutine get_primitive_vector(k, U, prim) 
!
			implicit none
!		
			real*8 :: U(5), prim(5), temp
			real*8 :: U2_rot, U3_rot, U4_rot
			integer :: k
!		
!			
			prim(1) = U(1)
			temp = 1.00/U(1)	
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
			prim(2) = temp*U(2)
			prim(3) = temp*U(3)
			prim(4) = temp*U(4)
!		
			temp = prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4)

			prim(5) = 0.40*(U(5) - 0.50*prim(1)*temp)
!	
		end subroutine
!		
!		
end module get_primitive_vector_mod
