module implicit_LdelU_forward_sweep_wall_mod

!	First written on 07.03.2021.
!	Revised on 16.04.2021.
!

        use data_structure_mod
        use get_primitive_vector_mod
        use octant_fluxes_mod
        use split_fluxes_mod    
!
!
contains
!
!
!	This subroutine evaluates the LdelU of a wall point during the forward sweep ..


        subroutine implicit_LdelU_forward_sweep_wall(i, LdelU)
!
!
		implicit none

		integer :: i, j, k
		real*8 :: prim(5)
		real*8 :: x_i, y_i, z_i, x_k, y_k, z_k
		real*8 :: tan1(3), tan2(3), nor(3)
		real*8 :: delx, dely, delz, det
		real*8 :: dels, delt, deln
!		
		real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delz_sqr
		real*8 :: sum_delx_dely, sum_dely_delz, sum_delz_delx
		real*8 :: sum_delx_delf(5), sum_dely_delf(5), sum_delz_delf(5)
		real*8 :: dist, weights
		real*8 :: temp(5), U(5), Gnew(5), Gold(5)
		real*8 :: dGxp(5), dGxn(5), dGyp(5), dGyn(5), dGzp(5), dGzn(5)
		real*8 :: LdelU(5)
		real*8 :: dels_weights, delt_weights, deln_weights
		integer:: alias_of_i, alias_of_k
!
!
		x_i = point%x(i)
		y_i = point%y(i)
		z_i = point%z(i)
!
		tan1 = point%tan1(:,i)
		tan2 = point%tan2(:,i)
		nor = point%nor(:,i)
!
!	Compute dGxpos contributions to LdelU ..
!
		sum_delx_sqr = 0.00
		sum_dely_sqr = 0.00
		sum_delz_sqr = 0.00
!
		sum_delx_dely = 0.00
		sum_dely_delz = 0.00
		sum_delz_delx = 0.00
!
		sum_delx_delf = 0.00
		sum_dely_delf = 0.00
		sum_delz_delf = 0.00
!
		do j = 1, point%xpos_nbhs(i)
			k = point%xpos_conn(j,i)
!
			x_k = point%x(k)
			y_k = point%y(k)
			z_k = point%z(k)
!
			delx = x_k - x_i
			dely = y_k - y_i
			delz = z_k - z_i
!
			dels = delx*tan1(1) + dely*tan1(2) + delz*tan1(3)
			delt = delx*tan2(1) + dely*tan2(2) + delz*tan2(3)
			deln = delx*nor(1) + dely*nor(2) + delz*nor(3)
!
			dist = sqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.00/(dist**power)
!
			dels_weights = dels*weights
			delt_weights = delt*weights
			deln_weights = deln*weights
!
			sum_delx_sqr = sum_delx_sqr + dels*dels_weights
			sum_dely_sqr = sum_dely_sqr + delt*delt_weights
			sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
			sum_delx_dely = sum_delx_dely + dels*delt_weights
			sum_dely_delz = sum_dely_delz + delt*deln_weights
			sum_delz_delx = sum_delz_delx + deln*dels_weights
!
			alias_of_i = point%alias(i)
			alias_of_k = point%alias(k)
			if(alias_of_k < alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gwxp(Gold, tan1, tan2, nor, prim)
				U = point%delUp(:,k) + point%U(:,k) 
				call get_primitive_vector(k, U, prim)
				call flux_Gwxp(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
	    enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!   	                     
		temp = sum_delx_delf*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_dely_delf*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
				+ sum_delz_delf*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
		dGxp = temp/det
!
!
!	Compute dGxneg contributions to LdelU ..
!
		sum_delx_sqr = 0.00
		sum_dely_sqr = 0.00
		sum_delz_sqr = 0.00
!
		sum_delx_dely = 0.00
		sum_dely_delz = 0.00
		sum_delz_delx = 0.00
!
		sum_delx_delf = 0.00
		sum_dely_delf = 0.00
		sum_delz_delf = 0.00
!
		do j = 1, point%xneg_nbhs(i)
			k = point%xneg_conn(j,i)
!								
			x_k = point%x(k)
			y_k = point%y(k)
			z_k = point%z(k)
!
			delx = x_k - x_i
			dely = y_k - y_i
			delz = z_k - z_i
!
			dels = delx*tan1(1) + dely*tan1(2) + delz*tan1(3)
			delt = delx*tan2(1) + dely*tan2(2) + delz*tan2(3)
			deln = delx*nor(1) + dely*nor(2) + delz*nor(3)
!
			dist = sqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.00/(dist**power)
!
			dels_weights = dels*weights
			delt_weights = delt*weights
			deln_weights = deln*weights
!
			sum_delx_sqr = sum_delx_sqr + dels*dels_weights
			sum_dely_sqr = sum_dely_sqr + delt*delt_weights
			sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
			sum_delx_dely = sum_delx_dely + dels*delt_weights
			sum_dely_delz = sum_dely_delz + delt*deln_weights
			sum_delz_delx = sum_delz_delx + deln*dels_weights
!
			alias_of_i = point%alias(i)
			alias_of_k = point%alias(k)
			if(alias_of_k < alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gwxn(Gold, tan1, tan2, nor, prim)				
				U = point%delUp(:,k) + point%U(:,k) 
				call get_primitive_vector(k, U, prim)
				call flux_Gwxn(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!                        
		temp = sum_delx_delf*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_dely_delf*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
				+ sum_delz_delf*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
		dGxn = temp/det

!
!	Compute dGypos contributions to LdelU ..
!
		sum_delx_sqr = 0.00
		sum_dely_sqr = 0.00
		sum_delz_sqr = 0.00
!
		sum_delx_dely = 0.00
		sum_dely_delz = 0.00
		sum_delz_delx = 0.00
!
		sum_delx_delf = 0.00
		sum_dely_delf = 0.00
		sum_delz_delf = 0.00
!
		do j = 1, point%ypos_nbhs(i)
			k = point%ypos_conn(j,i)
				
			x_k = point%x(k)
			y_k = point%y(k)
			z_k = point%z(k)
!
			delx = x_k - x_i
			dely = y_k - y_i
			delz = z_k - z_i
!
			dels = delx*tan1(1) + dely*tan1(2) + delz*tan1(3)
			delt = delx*tan2(1) + dely*tan2(2) + delz*tan2(3)
			deln = delx*nor(1) + dely*nor(2) + delz*nor(3)
!
			dist = sqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.00/(dist**power)
!
			dels_weights = dels*weights
			delt_weights = delt*weights
			deln_weights = deln*weights
!
			sum_delx_sqr = sum_delx_sqr + dels*dels_weights
			sum_dely_sqr = sum_dely_sqr + delt*delt_weights
			sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
			sum_delx_dely = sum_delx_dely + dels*delt_weights
			sum_dely_delz = sum_dely_delz + delt*deln_weights
			sum_delz_delx = sum_delz_delx + deln*dels_weights
!
			alias_of_i = point%alias(i)
			alias_of_k = point%alias(k)
			if(alias_of_k < alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gwyp(Gold, tan1, tan2, nor, prim)
				U = point%delUp(:,k) + point%U(:,k) 
				call get_primitive_vector(k, U, prim)
				call flux_Gwyp(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!
   	    temp = sum_delx_sqr*(sum_dely_delf*sum_delz_sqr - sum_dely_delz*sum_delz_delf) &
				- sum_delx_dely*(sum_delx_delf*sum_delz_sqr - sum_delz_delx*sum_delz_delf) &
				+ sum_delz_delx*(sum_delx_delf*sum_dely_delz - sum_delz_delx*sum_dely_delf)
!
		dGyp = temp/det
!
!
!	Compute dGyneg contributions to LdelU ..
!
		sum_delx_sqr = 0.00
		sum_dely_sqr = 0.00
		sum_delz_sqr = 0.00
!
		sum_delx_dely = 0.00
		sum_dely_delz = 0.00
		sum_delz_delx = 0.00
!
		sum_delx_delf = 0.00
		sum_dely_delf = 0.00
		sum_delz_delf = 0.00
!
		do j = 1, point%yneg_nbhs(i)
			k = point%yneg_conn(j,i)
!
				x_k = point%x(k)
				y_k = point%y(k)
				z_k = point%z(k)
!
				delx = x_k - x_i
				dely = y_k - y_i
				delz = z_k - z_i
!
				dels = delx*tan1(1) + dely*tan1(2) + delz*tan1(3)
				delt = delx*tan2(1) + dely*tan2(2) + delz*tan2(3)
				deln = delx*nor(1) + dely*nor(2) + delz*nor(3)
!
				dist = sqrt(dels*dels + delt*delt + deln*deln)
				weights = 1.00/(dist**power)
!
				dels_weights = dels*weights
				delt_weights = delt*weights
				deln_weights = deln*weights
!
				sum_delx_sqr = sum_delx_sqr + dels*dels_weights
				sum_dely_sqr = sum_dely_sqr + delt*delt_weights
				sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
				sum_delx_dely = sum_delx_dely + dels*delt_weights
				sum_dely_delz = sum_dely_delz + delt*deln_weights
				sum_delz_delx = sum_delz_delx + deln*dels_weights
!
				alias_of_i = point%alias(i)
				alias_of_k = point%alias(k)
				if(alias_of_k < alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gwyn(Gold, tan1, tan2, nor, prim)
				U = point%delUp(:,k) + point%U(:,k) 
				call get_primitive_vector(k, U, prim)
				call flux_Gwyn(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!  
   	    temp = sum_delx_sqr*(sum_dely_delf*sum_delz_sqr - sum_dely_delz*sum_delz_delf) &
				- sum_delx_dely*(sum_delx_delf*sum_delz_sqr - sum_delz_delx*sum_delz_delf) &
				+ sum_delz_delx*(sum_delx_delf*sum_dely_delz - sum_delz_delx*sum_dely_delf)
!
		dGyn = temp/det
!
!
!	Compute dGzneg contributions to LdelU ..
!
		sum_delx_sqr = 0.00
		sum_dely_sqr = 0.00
		sum_delz_sqr = 0.00
!
		sum_delx_dely = 0.00
		sum_dely_delz = 0.00
		sum_delz_delx = 0.00
!
		sum_delx_delf = 0.00
		sum_dely_delf = 0.00
		sum_delz_delf = 0.00
!
		do j = 1, point%nbhs(i)
			k = point%conn(j,i)
!
			x_k = point%x(k)
			y_k = point%y(k)
			z_k = point%z(k)
!
			delx = x_k - x_i
			dely = y_k - y_i
			delz = z_k - z_i
!
			dels = delx*tan1(1) + dely*tan1(2) + delz*tan1(3)
			delt = delx*tan2(1) + dely*tan2(2) + delz*tan2(3)
			deln = delx*nor(1) + dely*nor(2) + delz*nor(3)
!
			dist = sqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.00/(dist**power)
!
			dels_weights = dels*weights
			delt_weights = delt*weights
			deln_weights = deln*weights
!
			sum_delx_sqr = sum_delx_sqr + dels*dels_weights
			sum_dely_sqr = sum_dely_sqr + delt*delt_weights
			sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
			sum_delx_dely = sum_delx_dely + dels*delt_weights
			sum_dely_delz = sum_dely_delz + delt*deln_weights
			sum_delz_delx = sum_delz_delx + deln*dels_weights
!
			alias_of_i = point%alias(i)
			alias_of_k = point%alias(k)
			if(alias_of_k < alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gzn(Gold, tan1, tan2, nor, prim)
				U = point%delUp(:,k) + point%U(:,k) 
				call get_primitive_vector(k, U, prim)
				call flux_Gzn(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!   
		temp = sum_delx_sqr*(sum_dely_sqr*sum_delz_delf - sum_dely_delf*sum_dely_delz) &
 				- sum_delx_dely*(sum_delx_dely*sum_delz_delf - sum_delx_delf*sum_dely_delz) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delf - sum_delx_delf*sum_dely_sqr)
!
		dGzn = temp/det
!
!
		LdelU = 2.00*(dGxp + dGxn + dGyp + dGyn + dGzn)
!
!			
        end subroutine
!
!
end module implicit_LdelU_forward_sweep_wall_mod
