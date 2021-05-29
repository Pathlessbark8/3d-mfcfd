module outer_flux_dGxpos_mod

!	First written on 01.01.2021.
!	

        use data_structure_mod
        use octant_fluxes_mod    
        use q_variables_mod
        use q_derivatives_mod
        use qvariables_to_primitive_variables_mod
        use limiters_mod


contains


!	This subroutine evaluates the outer flux derivative dGx_pos


        subroutine outer_dGx_pos(G, i)
!
!
		implicit none

		integer :: i, j, k, r
		real*8 :: prim(5)
		real*8 :: x_i, y_i, z_i, x_k, y_k, z_k
		real*8 :: tan1(3), tan2(3), nor(3)
		real*8 :: G_i(5), G_k(5), G(5)
		real*8 :: delx, dely, delz, det, one_by_det
		real*8 :: dels, delt, deln
!		
		real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delz_sqr
		real*8 :: sum_delx_dely, sum_dely_delz, sum_delz_delx
		real*8 :: sum_delx_delf(5), sum_dely_delf(5), sum_delz_delf(5)
		real*8 :: dist, weights
		real*8 :: temp(5), qtilde(5), phi(5)
		real*8 :: dels_weights, delt_weights, deln_weights
!
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
		x_i = point%x(i)
		y_i = point%y(i)
		z_i = point%z(i)
!
		tan1 = point%tan1(:,i)
		tan2 = point%tan2(:,i)
		nor = point%nor(:,i)
!
		do j = 1, point%xpos_nbhs(i)
!
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
			temp = delx*point%dq(1,:,i) + dely*point%dq(2,:,i) + delz*point%dq(3,:,i)  
			qtilde = point%q(:,i) - 0.50*temp
			call venkat_limiter(qtilde, phi, i)
			qtilde = point%q(:,i) - 0.50*phi*temp
			call qtilde_to_primitive(qtilde, prim)
			call flux_Goxp(G_i, tan1, tan2, nor, prim)
!
            temp = delx*point%dq(1,:,k) + dely*point%dq(2,:,k) + delz*point%dq(3,:,k)  
			qtilde = point%q(:,k) - 0.50*temp
			call venkat_limiter(qtilde, phi, k)
			qtilde = point%q(:,k) - 0.50*phi*temp
			call qtilde_to_primitive(qtilde, prim)
			call flux_Goxp(G_k, tan1, tan2, nor, prim)
!       
			temp = G_k - G_i
!
			sum_delx_delf = sum_delx_delf + temp*dels_weights
			sum_dely_delf = sum_dely_delf + temp*delt_weights
			sum_delz_delf = sum_delz_delf + temp*deln_weights
!
        enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
			- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
			+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!                        
		one_by_det = 1.0/det
!
		temp = sum_delx_delf*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
			- sum_dely_delf*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
			+ sum_delz_delf*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
		G = temp*one_by_det
!
        end subroutine
!
!
end module outer_flux_dGxpos_mod
