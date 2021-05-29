module primitive_to_conserved_vector_mod
!
!
	use data_structure_mod
!
!
	contains
!
!
		subroutine primitive_to_conserved(k, U) 
!
			implicit none
!		
			real*8 :: U(5)
			real*8 :: temp1, temp2, temp3
			integer :: k
!		
!			
			U(1) = point%prim(1,k)
			temp1 = point%prim(1,k)*point%prim(2,k)
			temp2 = point%prim(1,k)*point%prim(3,k)
			temp3 = point%prim(1,k)*point%prim(4,k)
!			U(5) = 2.50*point%prim(5,k) + 0.50*U(1)*(temp1*temp1 + temp2*temp2 + temp3*temp3)
!
!			The momentum components in the rotational frame ..
!
			U(2) = temp1*point%tan1(1,k) + temp2*point%tan1(2,k) + temp3*point%tan1(3,k)
			U(3) = temp1*point%tan2(1,k) + temp2*point%tan2(2,k) + temp3*point%tan2(3,k)
			U(4) = temp1*point%nor(1,k) + temp2*point%nor(2,k) + temp3*point%nor(3,k)
!
			temp1 = point%prim(2,k)*point%prim(2,k) + point%prim(3,k)*point%prim(3,k) + point%prim(4,k)*point%prim(4,k)
			U(5) = 2.50*point%prim(5,k) + 0.50*U(1)*temp1
!	
		end subroutine
!		
!		
end module primitive_to_conserved_vector_mod
!	
