module implicit_UdelU_backward_sweep_interior_mod

!	First written on 07.03.2021.
!	Revised on 16.04.2021.
!	

        use data_structure_mod
        use get_primitive_vector_mod
        use split_fluxes_mod    
!
!
contains
!
!
!	This subroutine evaluates the UdelU of an interior point during the backward sweep ..


        subroutine implicit_UdelU_backward_sweep_interior(i, UdelU)
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
		real*8 :: UdelU(5)
		real*8 :: dels_weights, delt_weights, deln_weights
		integer :: alias_of_i, alias_of_k
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
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gxp(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
				call get_primitive_vector(k, U, prim)
				call flux_Gxp(Gnew, tan1, tan2, nor, prim)
				temp = Gnew - Gold
				if(i == 24476) then
!				print*, i, k, point%delUn(:,k), point%U(:,k)
				endif
!
				sum_delx_delf = sum_delx_delf + temp*dels_weights
				sum_dely_delf = sum_dely_delf + temp*delt_weights
				sum_delz_delf = sum_delz_delf + temp*deln_weights
			endif	
!
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
		if(i == 846) then
!			print*, i, det, temp
		endif
!
!	Compute dGxneg contributions to LdelU ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gxn(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
				call get_primitive_vector(k, U, prim)
				call flux_Gxn(Gnew, tan1, tan2, nor, prim)
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
!
!	Compute dGypos contributions to LdelU ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
!
		do j = 1, point%ypos_nbhs(i)
			k = point%ypos_conn(j,i)
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gyp(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
				call get_primitive_vector(k, U, prim)
				call flux_Gyp(Gnew, tan1, tan2, nor, prim)
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
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gyn(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
				call get_primitive_vector(k, U, prim)
				call flux_Gyn(Gnew, tan1, tan2, nor, prim)
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
!	Compute dGzpos contributions to LdelU ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
!
		do j = 1, point%zpos_nbhs(i)
			k = point%zpos_conn(j,i)
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
!
			dels_weights = dels*weights
			delt_weights = delt*weights
			deln_weights = deln*weights

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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gzp(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
				call get_primitive_vector(k, U, prim)
				call flux_Gzp(Gnew, tan1, tan2, nor, prim)
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
		dGzp = temp/det
!
!
!	Compute dGzneg contributions to LdelU ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx_delf = 0.0d0
		sum_dely_delf = 0.0d0
		sum_delz_delf = 0.0d0
!
		do j = 1, point%zneg_nbhs(i)
			k = point%zneg_conn(j,i)
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
			dist = dsqrt(dels*dels + delt*delt + deln*deln)
			weights = 1.0d0/(dist**power)
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
			if(alias_of_k > alias_of_i)	then
				prim = point%prim(:,k)
				call flux_Gzn(Gold, tan1, tan2, nor, prim)
				U = point%delUn(:,k) + point%U(:,k)
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
		UdelU = dGxp + dGxn + dGyp + dGyn + dGzp + dGzn
!
!			
        end subroutine
!
!
end module implicit_UdelU_backward_sweep_interior_mod
