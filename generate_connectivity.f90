!
!	First written on 03.02.2021
!
!
module generate_connectivity_mod
    
	use data_structure_mod
!
!
	contains
!
		subroutine generate_split_stencils()
!
			implicit none
!
			integer :: i, k

			do k = 1, interior_points
					i = interior_points_index(k)
					call get_interior_point_neighbours(i)
			enddo
!
			do k = 1, wall_points
					i = wall_points_index(k)
					call get_wall_point_neighbours(i)
			enddo
!
			do k = 1, outer_points
					i = outer_points_index(k)
					call get_outer_point_neighbours(i)
			enddo
!

    	end subroutine 
!
!
!
		subroutine get_interior_point_neighbours(i)
!
			implicit none
!
			integer i,k,r,count
			real*8 xi,yi,zi,xk,yk,zk
			real*8 delx,dely,delz 	
			real*8 dels,delt,deln
!
			xi = point%x(i)
			yi = point%y(i)
			zi = point%z(i)
!
			point%xpos_nbhs(i) = 0
			point%xneg_nbhs(i) = 0
			point%ypos_nbhs(i) = 0
			point%yneg_nbhs(i) = 0
			point%zpos_nbhs(i) = 0
			point%zneg_nbhs(i) = 0
!
			do r = 1, point%nbhs(i)
				k = point%conn(r,i)
				xk = point%x(k)
				yk = point%y(k)
				zk = point%z(k)
!
				delx = xk - xi
				dely = yk - yi
				delz = zk - zi
!
				dels = delx*point%tan1(1,i) + dely*point%tan1(2,i) + delz*point%tan1(3,i)
				delt = delx*point%tan2(1,i) + dely*point%tan2(2,i) + delz*point%tan2(3,i)
				deln = delx*point%nor(1,i) + dely*point%nor(2,i) + delz*point%nor(3,i)
!				
				if(dels .le. 0.00) then
					point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1
					count = point%xpos_nbhs(i)
					point%xpos_conn(count,i) = k
				endif
!
				if(dels .gt. 0.00) then
					point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1
					count = point%xneg_nbhs(i)
					point%xneg_conn(count,i) = k
				endif
!
				if(delt .le. 0.00) then
					point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1
					count = point%ypos_nbhs(i)
					point%ypos_conn(count,i) = k
				endif
!
				if(delt .gt. 0.00) then
					point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1
					count = point%yneg_nbhs(i)
					point%yneg_conn(count,i) = k
				endif
!
				if(deln .le. 0.00) then
					point%zpos_nbhs(i) = point%zpos_nbhs(i) + 1
					count = point%zpos_nbhs(i)
					point%zpos_conn(count,i) = k
				endif
!
				if(deln .gt. 0.00) then
					point%zneg_nbhs(i) = point%zneg_nbhs(i) + 1
					count = point%zneg_nbhs(i)
					point%zneg_conn(count,i) = k
				endif
!
        	enddo
!
!        
		end subroutine
!
!
!
		subroutine get_wall_point_neighbours(i)
!    
			implicit none
!    
			integer i,k,r,count
			real*8 xi,yi,zi,xk,yk,zk
			real*8 delx,dely,delz 	
			real*8 dels,delt,deln
!     
			xi = point%x(i)
			yi = point%y(i)
			zi = point%z(i)
!
			point%xpos_nbhs(i) = 0
			point%xneg_nbhs(i) = 0
			point%ypos_nbhs(i) = 0
			point%yneg_nbhs(i) = 0
			point%zneg_nbhs(i) = 0
!    
			do r = 1, point%nbhs(i)
				k = point%conn(r,i)
				xk = point%x(k)
				yk = point%y(k)
				zk = point%z(k)
!
				delx = xk - xi
				dely = yk - yi
				delz = zk - zi
!
				dels = delx*point%tan1(1,i) + dely*point%tan1(2,i) + delz*point%tan1(3,i)
				delt = delx*point%tan2(1,i) + dely*point%tan2(2,i) + delz*point%tan2(3,i)
				deln = delx*point%nor(1,i) + dely*point%nor(2,i) + delz*point%nor(3,i)
!
				if(dels .le. 0.00) then
					point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1
					count = point%xpos_nbhs(i)
					point%xpos_conn(count,i) = k
				endif
!
				if(dels .ge. 0.00) then
					point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1
					count = point%xneg_nbhs(i)
					point%xneg_conn(count,i) = k
				endif
!
				if(delt .le. 0.00) then
					point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1
					count = point%ypos_nbhs(i)
					point%ypos_conn(count,i) = k
				endif
!
				if(delt .ge. 0.00) then
					point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1
					count = point%yneg_nbhs(i)
					point%yneg_conn(count,i) = k
				endif
!
				point%zneg_nbhs(i) = point%zneg_nbhs(i) + 1
				count = point%zneg_nbhs(i)
				point%zneg_conn(count,i) = k
			enddo
!
!
		end subroutine
!
!
!
		subroutine get_outer_point_neighbours(i)
!
!    
			implicit none
!
			integer i,k,r,count
			real*8 xi,yi,zi,xk,yk,zk
			real*8 delx,dely,delz 	
			real*8 dels,delt,deln
!
			xi = point%x(i)
			yi = point%y(i)
			zi = point%z(i)
!
			point%xpos_nbhs(i) = 0
			point%xneg_nbhs(i) = 0
			point%ypos_nbhs(i) = 0
			point%yneg_nbhs(i) = 0
			point%zneg_nbhs(i) = 0
!    
			do r = 1, point%nbhs(i)
				k = point%conn(r,i)
				xk = point%x(k)
				yk = point%y(k)
				zk = point%z(k)
!
				delx = xk - xi
				dely = yk - yi
				delz = zk - zi
!
				dels = delx*point%tan1(1,i) + dely*point%tan1(2,i) + delz*point%tan1(3,i)
				delt = delx*point%tan2(1,i) + dely*point%tan2(2,i) + delz*point%tan2(3,i)
				deln = delx*point%nor(1,i) + dely*point%nor(2,i) + delz*point%nor(3,i)
!
				if(dels .le. 0.00) then
					point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1
					count = point%xpos_nbhs(i)
					point%xpos_conn(count,i) = k
				endif
!
				if(dels .ge. 0.00) then
					point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1
					count = point%xneg_nbhs(i)
					point%xneg_conn(count,i) = k
				endif
!
				if(delt .le. 0.00) then
					point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1
					count = point%ypos_nbhs(i)
					point%ypos_conn(count,i) = k
				endif
!
				if(delt .ge. 0.00) then
					point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1
					count = point%yneg_nbhs(i)
					point%yneg_conn(count,i) = k
				endif
!
				point%zpos_nbhs(i) = point%zpos_nbhs(i) + 1
				count = point%zpos_nbhs(i)
				point%zpos_conn(count,i) = k
			enddo
!
!
		end subroutine
!    
!
end module generate_connectivity_mod
