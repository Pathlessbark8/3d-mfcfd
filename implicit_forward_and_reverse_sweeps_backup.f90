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
		integer i, k, r, p, nbh
		real*8 :: dx, dy, dz, ds, min_dist
		real*8 :: D_inv, LdelU(5), UdelU(5)
!
!	We first perform the forward sweep ..
!
		do i = 1, wall_points
			k = wall_points_index(i)
			call implicit_D_inv_wall(k, D_inv)
			call implicit_LdelU_forward_sweep_wall(k, LdelU)
			point%delUp(:,k) = - D_inv*(point%flux_res(:,k) + LdelU)
		enddo	
!
		do i = 1, interior_points
			k = interior_points_index(i)
			call implicit_D_inv_interior(k, D_inv)
			call implicit_LdelU_forward_sweep_interior(k, LdelU)
			point%delUp(:,k) = - D_inv*(point%flux_res(:,k) + LdelU)
		enddo	
!
		do i = 1, outer_points
			k = outer_points_index(i)
			call implicit_D_inv_outer(k, D_inv)
			call implicit_LdelU_forward_sweep_outer(k, LdelU)
			point%delUp(:,k) = - D_inv*(point%flux_res(:,k) + LdelU)
		enddo	
!
		do i = 1, supersonic_inlet_points
			min_dist = 100000.0d0
			k = supersonic_inlet_points_index(i)
			do r = 1, point%nbhs(k)
				nbh = point%conn(r,k)
				dx = point%x(nbh) - point%x(k)
				dy = point%y(nbh) - point%y(k)
				dz = point%z(nbh) - point%z(k)
				ds = dsqrt(dx*dx + dy*dy + dz*dz)
				if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
					min_dist = ds
					p = nbh
				endif
			enddo
			point%delUp(:,k) = point%delUp(:,p)
		enddo
!
		do i = 1, supersonic_outlet_points
			min_dist = 100000.0d0
			k = supersonic_outlet_points_index(i)
			do r = 1, point%nbhs(k)
				nbh = point%conn(r,k)
				dx = point%x(nbh) - point%x(k)
				dy = point%y(nbh) - point%y(k)
				dz = point%z(nbh) - point%z(k)
				ds = dsqrt(dx*dx + dy*dy + dz*dz)
				if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
					min_dist = ds
					p = nbh
				endif
			enddo
			point%delUp(:,k) = point%delUp(:,p)
		enddo
!
!

!	The following steps perform the reverse sweep ..	  	
!
		do i = wall_points, 1, -1
			k = wall_points_index(i)
			call implicit_D_inv_wall(k, D_inv)
			call implicit_UdelU_backward_sweep_wall(k, UdelU)
			point%delUn(:,k) = point%delUp(:,k) - D_inv*UdelU
		enddo	
!
		do i = interior_points, 1, -1
			k = interior_points_index(i)
			call implicit_D_inv_interior(k, D_inv)
			call implicit_UdelU_backward_sweep_interior(k, UdelU)
			point%delUn(:,k) = point%delUp(:,k) - D_inv*UdelU
!			print*, k, UdelU 
		enddo	
!
		do i = outer_points, 1, -1
			k = outer_points_index(i)
			call implicit_D_inv_outer(k, D_inv)
			call implicit_UdelU_backward_sweep_outer(k, UdelU)
			point%delUn(:,k) = point%delUp(:,k) - D_inv*UdelU
		enddo	
!
		do i = supersonic_inlet_points, 1, -1
			min_dist = 100000.0d0
			k = supersonic_inlet_points_index(i)
			do r = 1, point%nbhs(k)
				nbh = point%conn(r,k)
				dx = point%x(nbh) - point%x(k)
				dy = point%y(nbh) - point%y(k)
				dz = point%z(nbh) - point%z(k)
				ds = dsqrt(dx*dx + dy*dy + dz*dz)
				if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
					min_dist = ds
					p = nbh
				endif
			enddo
			point%delUn(:,k) = point%delUn(:,p)
		enddo
!
		do i = supersonic_outlet_points, 1, -1
			min_dist = 100000.0d0
			k = supersonic_outlet_points_index(i)
			do r = 1, point%nbhs(k)
				nbh = point%conn(r,k)
				dx = point%x(nbh) - point%x(k)
				dy = point%y(nbh) - point%y(k)
				dz = point%z(nbh) - point%z(k)
				ds = dsqrt(dx*dx + dy*dy + dz*dz)
				if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
					min_dist = ds
					p = nbh
				endif
			enddo
			point%delUn(:,k) = point%delUn(:,p)
		enddo
!
	end subroutine 
end module implicit_forward_and_reverse_sweeps_mod

