module implicit_forward_and_reverse_sweeps_mod
!
	use data_structure_mod
	use implicit_D_inv_wall_mod
	use implicit_D_inv_outer_mod
	use implicit_D_inv_interior_mod
	use implicit_LdelU_forward_sweep_wall_mod
	use implicit_LdelU_forward_sweep_outer_mod
	use implicit_LdelU_forward_sweep_interior_mod
	use implicit_UdelU_backward_sweep_wall_mod
	use implicit_UdelU_backward_sweep_outer_mod
	use implicit_UdelU_backward_sweep_interior_mod
!
	contains
!
	subroutine implicit_forward_backward_sweeps()
!
		implicit none
!
		integer i, k, r, p, q, nbh, point_status
		real*8 :: dx, dy, dz, ds, min_dist
		real*8 :: D_inv, LdelU(5), UdelU(5)
!
!	We first perform the forward sweep ..
!
		do i = 1, max_points
			p = point%point_with_alias(i)
			point_status = point%status(p)
!
			if(point_status == 1) then
				call implicit_D_inv_wall(p, D_inv)
				call implicit_LdelU_forward_sweep_wall(p, LdelU)
				point%delUp(:,p) = - D_inv*(point%flux_res(:,p) + LdelU)
			else if(point_status == 0) then
				call implicit_D_inv_interior(p, D_inv)
				call implicit_LdelU_forward_sweep_interior(p, LdelU)
				point%delUp(:,p) = - D_inv*(point%flux_res(:,p) + LdelU)
			else if(point_status == 5 .OR. point_status == 6) then ! supersonic inlet/outlet points ..
				min_dist = 1000000.0d0
				do r = 1, point%nbhs(p)
					nbh = point%conn(r,p)
					dx = point%x(nbh) - point%x(p)
					dy = point%y(nbh) - point%y(p)
					dz = point%z(nbh) - point%z(p)
					ds = dsqrt(dx*dx + dy*dy + dz*dz)
					if(ds < min_dist) then
						min_dist = ds
						q = nbh
					endif
				enddo
				point%delUp(:,p) = point%delUp(:,q)
			endif	
		enddo	
!
!
!	The following steps perform the reverse sweep ..	  	
! 
		do i = max_points, 1, -1
			p = point%point_with_alias(i)
			point_status = point%status(p)
			if(point_status == 1) then
				call implicit_D_inv_wall(p, D_inv)
				call implicit_UdelU_backward_sweep_wall(p, UdelU)
				point%delUn(:,p) = point%delUp(:,p) - D_inv*UdelU
			else if(point_status == 0) then
				call implicit_D_inv_interior(p, D_inv)
				call implicit_UdelU_backward_sweep_interior(p, UdelU)
				point%delUn(:,p) = point%delUp(:,p) - D_inv*UdelU
			else if(point_status == 5 .OR. point_status == 6) then ! supersonic inlet/outlet points ..
				min_dist = 1000000.0d0
				do r = 1, point%nbhs(p)
					nbh = point%conn(r,p)
					dx = point%x(nbh) - point%x(p)
					dy = point%y(nbh) - point%y(p)
					dz = point%z(nbh) - point%z(p)
					ds = dsqrt(dx*dx + dy*dy + dz*dz)
					if(ds < min_dist) then
						min_dist = ds
						q = nbh
					endif
				enddo
				point%delUn(:,p) = point%delUn(:,q)
			endif		
		enddo	
!
		do i = 1, max_points
			point%delUp(:,i) = point%delUn(:,i)
		enddo	
	end subroutine 
end module implicit_forward_and_reverse_sweeps_mod

