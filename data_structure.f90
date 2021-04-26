module data_structure_mod
!
	use parameter_mod
!
	implicit none
!
	type :: points
!
!		geometry based attributes ..
!
		real*8, dimension(:), allocatable :: x,y,z
!		real*8, dimension(:,:), allocatable :: nor_neigh
		real*8, dimension(:,:), allocatable :: tan1, tan2, nor
		integer, dimension(:), allocatable :: status
		integer, dimension(:), allocatable :: nbhs
		integer, dimension(:,:), allocatable :: conn
		real*8, dimension(:), allocatable :: min_dist
!
		integer, dimension(:), allocatable :: xpos_nbhs,xneg_nbhs
		integer, dimension(:), allocatable :: ypos_nbhs,yneg_nbhs
		integer, dimension(:), allocatable :: zpos_nbhs,zneg_nbhs
		integer, dimension(:,:), allocatable :: xpos_conn,xneg_conn
		integer, dimension(:,:), allocatable :: ypos_conn,yneg_conn
		integer, dimension(:,:), allocatable :: zpos_conn,zneg_conn
!
!		flow field based attributes ..
!
		real*8, dimension(:,:), allocatable :: U
		real*8, dimension(:,:), allocatable :: delUp
		real*8, dimension(:,:), allocatable :: delUn
		real*8, dimension(:,:), allocatable :: prim
		real*8, dimension(:,:), allocatable :: flux_res
		real*8, dimension(:,:), allocatable :: q
		real*8, dimension(:,:,:), allocatable :: qm
		real*8, dimension(:,:,:), allocatable :: dq
		real*8, dimension(:,:,:), allocatable :: temp
		real*8, dimension(:), allocatable :: delt
!
!		for the implicit solver ..
!
		integer, dimension(:), allocatable :: alias, point_with_alias
!
!
	end type points
!
	type(points) :: point
!
!
	integer interior_points, wall_points, outer_points
	integer supersonic_inlet_points, supersonic_outlet_points
!
	integer, dimension(:), allocatable :: interior_points_index
	integer, dimension(:), allocatable :: wall_points_index
	integer, dimension(:), allocatable :: outer_points_index
	integer, dimension(:), allocatable :: supersonic_outlet_points_index
	integer, dimension(:), allocatable :: supersonic_inlet_points_index
!
	real*8 res_old, res_new, max_res, residue
	integer max_res_point
!	 
	contains
!
	subroutine allocate_size()
!	        
			implicit none
!
			allocate(point%x(max_points))
			allocate(point%y(max_points))
			allocate(point%z(max_points))
			allocate(point%tan1(3,max_points))
			allocate(point%tan2(3,max_points))
			allocate(point%nor(3,max_points))              
			allocate(point%status(max_points))
			allocate(point%nbhs(max_points))
			allocate(point%conn(27,max_points))
			allocate(point%min_dist(max_points))
			allocate(point%xpos_nbhs(max_points))
			allocate(point%xneg_nbhs(max_points))
			allocate(point%ypos_nbhs(max_points))
			allocate(point%yneg_nbhs(max_points))
			allocate(point%zpos_nbhs(max_points))
			allocate(point%zneg_nbhs(max_points))
			allocate(point%xpos_conn(22,max_points))
			allocate(point%xneg_conn(22,max_points))
			allocate(point%ypos_conn(22,max_points))
			allocate(point%yneg_conn(22,max_points))
			allocate(point%zpos_conn(22,max_points))
			allocate(point%zneg_conn(22,max_points))
			allocate(point%U(5,max_points))
			allocate(point%delUp(5,max_points))
			allocate(point%delUn(5,max_points))
			allocate(point%prim(5,max_points))
			allocate(point%flux_res(5,max_points))
			allocate(point%delt(max_points))
			allocate(point%q(5,max_points))
			allocate(point%qm(2,5,max_points))
			allocate(point%dq(3,5,max_points))
			allocate(point%temp(3,5,max_points))
!
			allocate(point%alias(max_points))
			allocate(point%point_with_alias(max_points))

!
	end subroutine
		

end module data_structure_mod
