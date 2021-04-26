module timestep_delt_mod
!
!
	use data_structure_mod
!
	contains	

	subroutine timestep_delt()
!
!
		implicit none
!
        integer :: i, k, r
		real*8 :: rho, u1, u2, u3, pr
		real*8 :: delx, dely, delz, dist		
		real*8 :: min_delt
		real*8 :: mod_u, delt
!
!
		do i = 1, max_points
!
			min_delt = 1.0d0						
!			
			do r = 1, point%nbhs(i)
!			
                k = point%conn(r,i)
!
                rho = point%prim(1,k)
                u1 = point%prim(2,k)
                u2 = point%prim(3,k)
				u3 = point%prim(4,k)
                pr = point%prim(5,k)
!
				delx = point%x(k) - point%x(i)
				dely = point%y(k) - point%y(i)
				delz = point%z(k) - point%z(i)
!
                dist = dsqrt(delx*delx + dely*dely + delz*delz)
!
                mod_u = dsqrt(u1*u1 + u2*u2 + u3*u3)

                delt = cfl*dist/(mod_u + 3.0d0*dsqrt(pr/rho))
!
                if(delt < min_delt) then
					min_delt = delt
                endif
!
			enddo
!
			point%delt(i) = min_delt						    
!
		end do
!		
	end subroutine
!
!
end module timestep_delt_mod
