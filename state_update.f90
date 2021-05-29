module state_update_mod
!
!
!
!
	use data_structure_mod
	use primitive_to_conserved_vector_mod
	use conserved_to_primitive_vector_mod
	use conserved_vector_Ubar_mod
!
!
	contains
!
!
		subroutine state_update()
		!
		!	
			implicit none			
			!
				integer :: i, k, nbh, p, r
				real*8 :: U(5), temp, min_dist
				real*8 :: res_sqr, sum_res_sqr
				real*8 :: dx,dy,dz,ds
!			
				sum_res_sqr = 0.00
				max_res = 0.00	
!
!
				do i =	1, wall_points
!
					k = wall_points_index(i)
!					call primitive_to_conserved(k, U)
!					temp = U(1)
					temp = point%U(1,k)
!					U = U - point%delt(k)*point%flux_res(:,k)
				        U = point%U(:,k) + point%delUn(:,k)
					U(4) = 0.00
!
					res_sqr = (U(1) - temp)*(U(1) - temp)					
					if(res_sqr .gt. max_res) then 
						max_res = res_sqr
						max_res_point = k
					endif							
					sum_res_sqr = sum_res_sqr + res_sqr
!
					call conserved_to_primitive(k, U)			
					point%delUp(:,k) = point%delUn(:,k)
!
				enddo
!
				do i =	1, outer_points
!
					k = outer_points_index(i)
					call conserved_vector_Ubar(k, U) 
					temp = U(1)
					U = U - point%delt(k)*point%flux_res(:,k)
!
					call conserved_to_primitive(k, U)			
!
				enddo
!
				do i = 1, interior_points
					k = interior_points_index(i)	
!					call primitive_to_conserved(k, U)
					temp = point%U(1,k)
			        U = point%U(:,k) + point%delUn(:,k)					
!					U = U - point%delt(k)*point%flux_res(:,k)
!
					res_sqr = (U(1) - temp)*(U(1) - temp)
					if(res_sqr .gt. max_res) then 
						max_res = res_sqr
						max_res_point = k
					endif
					sum_res_sqr = sum_res_sqr + res_sqr
!
					call conserved_to_primitive(k, U)			
					point%delUp(:,k) = point%delUn(:,k)					
				enddo
!				
				do i = 1, supersonic_outlet_points
					min_dist = 100000.00
					k = supersonic_outlet_points_index(i)
						do r = 1, point%nbhs(k)
							nbh = point%conn(r,k)
							dx = point%x(nbh) - point%x(k)
							dy = point%y(nbh) - point%y(k)
							dz = point%z(nbh) - point%z(k)
							ds = sqrt(dx*dx + dy*dy + dz*dz)
!							if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
							if(ds < min_dist) then
								min_dist = ds
								p = nbh
							endif
						enddo
					point%prim(:,k) = point%prim(:,p)
				enddo	
!
				do i = 1, supersonic_inlet_points
					min_dist = 100000.00
					k = supersonic_inlet_points_index(i)
						do r = 1, point%nbhs(k)
							nbh = point%conn(r,k)
							dx = point%x(nbh) - point%x(k)
							dy = point%y(nbh) - point%y(k)
							dz = point%z(nbh) - point%z(k)
							ds = sqrt(dx*dx + dy*dy + dz*dz)
!							if(ds < min_dist .AND. point%status(nbh) .ne. 1) then
							if(ds < min_dist) then
								min_dist = ds
								p = nbh
							endif
						enddo
					point%prim(:,k) = point%prim(:,p)
				enddo	
				res_new = sqrt(sum_res_sqr)/max_points
!
!
		end subroutine 							
!
!
		!	
!
end module state_update_mod	
!	
