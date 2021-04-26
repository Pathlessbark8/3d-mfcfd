!
!	First written on 02.02.2021 ..
!
!
module point_preprocessor_mod
!
!
	use data_structure_mod
!
!
contains
!
	subroutine read_input_point_data()
!
!
		implicit none		
!
		integer k,r,counter
		integer interior_count, wall_count, outer_count
		integer supersonic_inlet_count, supersonic_outlet_count
!		
!
!
!		OPEN(UNIT=101,FILE="hemisphere-testcase/3d_input_data",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
		OPEN(UNIT=101,FILE="3d_input_data",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
!
!
!		read(101,*) max_points
!
		do k = 1, max_points
			read(101,*) counter, point%status(k), point%x(k), point%y(k), point%z(k), &
			& point%tan1(1,k), point%tan1(2,k), point%tan1(3,k), &
			& point%tan2(1,k), point%tan2(2,k), point%tan2(3,k), &
			& point%nor(1,k), point%nor(2,k), point%nor(3,k), &
			& point%nbhs(k), (point%conn(r,k),r=1,point%nbhs(k)), point%min_dist(k)
		enddo	
!
!		Finding the number of interior, wall, outer and other boundary points ..
!
!
		interior_points = 0
		wall_points = 0
		outer_points = 0
		supersonic_outlet_points = 0
		supersonic_inlet_points = 0
!
		do k = 1, max_points
			if(point%status(k) == 0) then
				interior_points = interior_points + 1
			else if(point%status(k) == 1) then 	
				wall_points = wall_points + 1
			else if(point%status(k) == 2) then
				outer_points = outer_points + 1
			else if(point%status(k) == 6) then
				supersonic_outlet_points = supersonic_outlet_points + 1
			else if(point%status(k) == 5) then
				supersonic_inlet_points = supersonic_inlet_points + 1
			endif
		enddo	
!
!		Allocating the size of the respective points ..
!
		allocate(interior_points_index(interior_points))
		allocate(wall_points_index(wall_points))
		allocate(outer_points_index(outer_points))
		allocate(supersonic_outlet_points_index(supersonic_outlet_points))
		allocate(supersonic_inlet_points_index(supersonic_inlet_points))
!		
!		
!		Finding the indices of the interior, wall, outer and other boundary points ..
!
		interior_count = 0
		wall_count = 0
		outer_count = 0
		supersonic_inlet_count = 0
		supersonic_outlet_count = 0
!		
		do k = 1, max_points
			if(point%status(k) == 0) then
				interior_count = interior_count + 1
				interior_points_index(interior_count) = k
			else if(point%status(k) == 1) then
				wall_count = wall_count + 1
				wall_points_index(wall_count) = k
			else if(point%status(k) == 2) then
				outer_count = outer_count + 1
				outer_points_index(outer_count) = k
			else if(point%status(k) == 6) then
				supersonic_outlet_count = supersonic_outlet_count + 1
				supersonic_outlet_points_index(supersonic_outlet_count) = k
			else if(point%status(k) == 5) then
				supersonic_inlet_count = supersonic_inlet_count + 1
				supersonic_inlet_points_index(supersonic_inlet_count) = k
			endif
		enddo
!							
!				
	end subroutine
!
!
end module point_preprocessor_mod 						
					 							
			
