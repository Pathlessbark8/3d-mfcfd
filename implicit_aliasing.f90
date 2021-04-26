module implicit_aliasing_mod
!
!	First written on 23.04.2021.
!
	use data_structure_mod
!
!
contains
!
!	This subroutine sets the alias indexing for all points ..
!
!	
		subroutine aliasing()
!
!
			implicit none
!
			integer:: i, j, p, k, count
!
!
			do i = 1, max_points
				point%alias(i) = 0
				point%point_with_alias(i) = 0
			enddo
!
!	
			point%alias(1) = 1
			point%point_with_alias(1) = 1
			count = 0
!
!
			do i = 1, max_points
				p = point%point_with_alias(i)
				do j = 1, point%nbhs(p)
					k = point%conn(j,p)
					if(point%alias(k) == 0) then 
						count = count + 1
						point%alias(k) = count
						point%point_with_alias(count) = k
					endif
				enddo
			enddo			
!
		end subroutine
!
end module implicit_aliasing_mod						
					
