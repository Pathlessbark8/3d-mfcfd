module q_derivatives_mod

       
        use data_structure_mod
        use q_variables_mod

contains
!
!
		subroutine eval_q_derivatives()
!
!
			implicit none
!
!
			integer :: i, j, k, r, nbh
			real*8 :: x_i, y_i, z_i, x_k, y_k, z_k
			real*8 :: delx, dely, delz, dist, weights
			real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delz_sqr
			real*8 :: sum_delx_dely, sum_dely_delz, sum_delz_delx
			real*8 :: sum_delx_delq(5), sum_dely_delq(5), sum_delz_delq(5)
			real*8 :: det, one_by_det
			real*8 :: temp(5), qtilde_i(5), qtilde_nbh(5)
!
			do i = 1, max_points
!
				x_i = point%x(i)
				y_i = point%y(i)
				z_i = point%z(i)
!			
				sum_delx_sqr = 0.0
        	    sum_dely_sqr = 0.0
            	sum_delz_sqr = 0.0
!            
	            sum_delx_dely = 0.0
    	        sum_dely_delz = 0.0
        	    sum_delz_delx = 0.0           
!
	 			sum_delx_delq = 0.0
    	        sum_dely_delq = 0.0
        	    sum_delz_delq = 0.0
!
	 			point%qm(1, :, i) = point%q(:, i)	! q_maximum ..
    	        point%qm(2, :, i) = point%q(:, i)	! q_minimum ..
!
	            do k = 1, point%nbhs(i)
!
					nbh = point%conn(k,i)
!
					do r = 1, 5
						if(point%q(r,nbh) > point%qm(1,r,i)) then
							point%qm(1,r,i) = point%q(r,nbh)
						end if
						if(point%q(r,nbh) < point%qm(2,r,i)) then
							point%qm(2,r,i) = point%q(r,nbh)
						end if
					end do
!
					x_k = point%x(nbh)
					y_k = point%y(nbh)
					z_k = point%z(nbh)
!			
					delx = x_k - x_i
					dely = y_k - y_i
					delz = z_k - z_i
!			
					dist = sqrt(delx*delx + dely*dely + delz*delz)
					weights = 1.00/(dist**power)
!
					sum_delx_sqr = sum_delx_sqr + delx*delx*weights
					sum_dely_sqr = sum_dely_sqr + dely*dely*weights
					sum_delz_sqr = sum_delz_sqr + delz*delz*weights
!			
					sum_delx_dely = sum_delx_dely + delx*dely*weights
					sum_dely_delz = sum_dely_delz + dely*delz*weights
					sum_delz_delx = sum_delz_delx + delz*delx*weights
!			
					temp = (point%q(:,nbh) - point%q(:,i))
					sum_delx_delq = sum_delx_delq + weights*delx*temp
					sum_dely_delq = sum_dely_delq + weights*dely*temp
					sum_delz_delq = sum_delz_delq + weights*delz*temp
!
				enddo
!
				det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
					- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
					+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)					
!				
				one_by_det = 1.00/det
!
				temp = sum_delx_delq*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
					 - sum_dely_delq*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
					 + sum_delz_delq*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
				point%dq(1,:,i) = temp*one_by_det
!
				temp = sum_delx_sqr*(sum_dely_delq*sum_delz_sqr - sum_dely_delz*sum_delz_delq) &
						- sum_delx_dely*(sum_delx_delq*sum_delz_sqr - sum_delz_delx*sum_delz_delq) &
						+ sum_delz_delx*(sum_delx_delq*sum_dely_delz - sum_delz_delx*sum_dely_delq)				
!
				point%dq(2,:,i) = temp*one_by_det
!
				temp = sum_delx_sqr*(sum_dely_sqr*sum_delz_delq - sum_dely_delq*sum_dely_delz) &
						- sum_delx_dely*(sum_delx_dely*sum_delz_delq - sum_delx_delq*sum_dely_delz) &
						+ sum_delz_delx*(sum_delx_dely*sum_dely_delq - sum_delx_delq*sum_dely_sqr)
!
				point%dq(3,:,i) = temp*one_by_det
!
!
			enddo
!
!
!	Inner iterations to compute second order accurate q-derivatives ..
!
!
			do j = 1, inner_iterations
!			
				do i = 1, max_points
!
					x_i = point%x(i)
					y_i = point%y(i)
					z_i = point%z(i)
!			
					sum_delx_sqr = 0.0
					sum_dely_sqr = 0.0
					sum_delz_sqr = 0.0
!            
					sum_delx_dely = 0.0
					sum_dely_delz = 0.0
					sum_delz_delx = 0.0           
!
					sum_delx_delq = 0.0
					sum_dely_delq = 0.0
					sum_delz_delq = 0.0
!
					do k = 1, point%nbhs(i)
!
						nbh = point%conn(k,i)
!
						x_k = point%x(nbh)
						y_k = point%y(nbh)
						z_k = point%z(nbh)
!			
						delx = x_k - x_i
						dely = y_k - y_i
						delz = z_k - z_i
!		
						dist = sqrt(delx*delx + dely*dely + delz*delz)
						weights = 1.00/(dist**power)
!
						sum_delx_sqr = sum_delx_sqr + delx*delx*weights
						sum_dely_sqr = sum_dely_sqr + dely*dely*weights
						sum_delz_sqr = sum_delz_sqr + delz*delz*weights
!			
						sum_delx_dely = sum_delx_dely + delx*dely*weights
						sum_dely_delz = sum_dely_delz + dely*delz*weights
						sum_delz_delx = sum_delz_delx + delz*delx*weights
!			
						qtilde_i = point%q(:,i) - 0.50*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i) + delz*point%dq(3,:,i))
						qtilde_nbh = point%q(:,nbh) - 0.50*(delx*point%dq(1,:,nbh) + dely*point%dq(2,:,nbh) + delz*point%dq(3,:,nbh)) 						
						temp = qtilde_nbh - qtilde_i
!						
						sum_delx_delq = sum_delx_delq + weights*delx*temp
						sum_dely_delq = sum_dely_delq + weights*dely*temp
						sum_delz_delq = sum_delz_delq + weights*delz*temp
!
					enddo
!
					det = sum_delx_sqr*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
						- sum_delx_dely*(sum_delx_dely*sum_delz_sqr - sum_dely_delz*sum_delz_delx) &
						+ sum_delz_delx*(sum_delx_dely*sum_dely_delz - sum_dely_sqr*sum_delz_delx)					
!				
					one_by_det = 1.00/det
!
					temp = sum_delx_delq*(sum_dely_sqr*sum_delz_sqr - sum_dely_delz*sum_dely_delz) &
						- sum_dely_delq*(sum_delx_dely*sum_delz_sqr - sum_delz_delx*sum_dely_delz) &
						+ sum_delz_delq*(sum_delx_dely*sum_dely_delz - sum_delz_delx*sum_dely_sqr)
!
					point%temp(1,:,i) = temp*one_by_det
!
					temp = sum_delx_sqr*(sum_dely_delq*sum_delz_sqr - sum_dely_delz*sum_delz_delq) &
						- sum_delx_dely*(sum_delx_delq*sum_delz_sqr - sum_delz_delx*sum_delz_delq) &
						+ sum_delz_delx*(sum_delx_delq*sum_dely_delz - sum_delz_delx*sum_dely_delq)				
!
					point%temp(2,:,i) = temp*one_by_det
!
					temp = sum_delx_sqr*(sum_dely_sqr*sum_delz_delq - sum_dely_delq*sum_dely_delz) &
						- sum_delx_dely*(sum_delx_dely*sum_delz_delq - sum_delx_delq*sum_dely_delz) &
						+ sum_delz_delx*(sum_delx_dely*sum_dely_delq - sum_delx_delq*sum_dely_sqr)
!
					point%temp(3,:,i) = temp*one_by_det
!
				enddo	!	i loop (for points)..
!
				do i = 1, max_points
					point%dq(1,:,i) = point%temp(1,:,i)
					point%dq(2,:,i) = point%temp(2,:,i)
					point%dq(3,:,i) = point%temp(3,:,i)
				enddo					
!				
			enddo	!	j loop (for inner iterations)..			
!
!
		end subroutine 


end module q_derivatives_mod
