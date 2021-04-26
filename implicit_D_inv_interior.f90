module implicit_D_inv_interior_mod

!	First written on 08.03.2021.
!	Revised on 07.04.2021.
!	

        use data_structure_mod
!
!
contains
!
!
!	This subroutine evaluates the D_inv of an interior point during
!	the forward and backward sweeps ..


        subroutine implicit_D_inv_interior(i, D_inv)
!
!
		implicit none

		integer :: i, j, k
		real*8 :: x_i, y_i, z_i, x_k, y_k, z_k
		real*8 :: tan1(3), tan2(3), nor(3)
		real*8 :: delx, dely, delz, det
		real*8 :: dels, delt, deln
!		
		real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delz_sqr
		real*8 :: sum_delx_dely, sum_dely_delz, sum_delz_delx
		real*8 :: sum_delx, sum_dely, sum_delz
		real*8 :: dist, weights
		real*8 :: temp, as
		real*8 :: dGxp, dGxn, dGyp, dGyn, dGzp, dGzn
		real*8 :: delta_t, D_inv, temp1
		real*8 :: dels_weights, delt_weights, deln_weights
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
		delta_t = point%delt(i)
		as = dsqrt(1.40d0*point%prim(5,i)/point%prim(1,i))
!
!	Compute dGxpos contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
!
	    enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!                        
		temp = sum_delx*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_dely*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
				+ sum_delz*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
		temp1 = point%prim(2,i) + as
		dGxp = 0.50d0*dabs(temp1 + dabs(temp1))*temp/det
!
!
!	Compute dGxneg contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!                        
		temp = sum_delx*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_dely*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
				+ sum_delz*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
		temp1 = point%prim(2,i) - as
		dGxn = -0.50d0*dabs(temp1 - dabs(temp1))*temp/det
!
!
!	Compute dGypos contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!
        temp = sum_delx_sqr*(sum_dely*sum_delz_sqr - sum_dely_delz*sum_delz) &
				- sum_delx_dely*(sum_delx*sum_delz_sqr - sum_delz_delx*sum_delz) &
				+ sum_delz_delx*(sum_delx*sum_dely_delz - sum_delz_delx*sum_dely)
!
		temp1 = point%prim(3,i) + as
		dGyp = 0.50d0*dabs(temp1 + dabs(temp1))*temp/det
!
!
!	Compute dGyneg contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!  
        temp = sum_delx_sqr*(sum_dely*sum_delz_sqr - sum_dely_delz*sum_delz) &
				- sum_delx_dely*(sum_delx*sum_delz_sqr - sum_delz_delx*sum_delz) &
				+ sum_delz_delx*(sum_delx*sum_dely_delz - sum_delz_delx*sum_dely)
!
		temp1 = point%prim(3,i) - as
		dGyn = -0.50d0*dabs(temp1 - dabs(temp1))*temp/det
!
!
!	Compute dGzpos contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
!
				sum_delx_sqr = sum_delx_sqr + dels*dels_weights
				sum_dely_sqr = sum_dely_sqr + delt*delt_weights
				sum_delz_sqr = sum_delz_sqr + deln*deln_weights
!
				sum_delx_dely = sum_delx_dely + dels*delt_weights
				sum_dely_delz = sum_dely_delz + delt*deln_weights
				sum_delz_delx = sum_delz_delx + deln*dels_weights
!
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!
		temp = sum_delx_sqr*(sum_dely_sqr*sum_delz - sum_dely*sum_dely_delz) &
 				- sum_delx_dely*(sum_delx_dely*sum_delz - sum_delx*sum_dely_delz) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely - sum_delx*sum_dely_sqr)
!
		temp1 = point%prim(4,i) + as
		dGzp = 0.50d0*dabs(temp1 + dabs(temp1))*temp/det
!
!
!	Compute dGzneg contributions to D_inv ..
!
		sum_delx_sqr = 0.0d0
		sum_dely_sqr = 0.0d0
		sum_delz_sqr = 0.0d0
!
		sum_delx_dely = 0.0d0
		sum_dely_delz = 0.0d0
		sum_delz_delx = 0.0d0
!
		sum_delx = 0.0d0
		sum_dely = 0.0d0
		sum_delz = 0.0d0
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
				sum_delx = sum_delx + dels_weights
				sum_dely = sum_dely + delt_weights
				sum_delz = sum_delz + deln_weights
		enddo
!
		det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
				- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)
!   
		temp = sum_delx_sqr*(sum_dely_sqr*sum_delz - sum_dely*sum_dely_delz) &
 				- sum_delx_dely*(sum_delx_dely*sum_delz - sum_delx*sum_dely_delz) &
				+ sum_delz_delx*(sum_delx_dely*sum_dely - sum_delx*sum_dely_sqr)
!
		temp1 = point%prim(4,i) - as
		dGzn = -0.50d0*dabs(temp1 - dabs(temp1))*temp/det
!
!
		temp = (1.0d0/delta_t) - (dGxp + dGxn + dGyp + dGyn + dGzp + dGzn)
		D_inv = 1.0d0/temp
!			
        end subroutine
!
!
end module implicit_D_inv_interior_mod
