module q_variables_mod

       
        use data_structure_mod

contains



		subroutine eval_q_variables()

			implicit none

			integer :: k
			real*8 :: rho, u1, u2, u3, pr, beta
			real*8 :: two_times_beta
!
!
			do k = 1, max_points
!				
				rho = point%prim(1,k)
				u1 = point%prim(2,k)
				u2 = point%prim(3,k)
				u3 = point%prim(4,k)                
				pr = point%prim(5,k)
!
				beta = 0.5d0*rho/pr
!
				point%q(1,k) = dlog(rho) + (dlog(beta)*2.5d0) - beta*(u1*u1 + u2*u2 + u3*u3)
!
				two_times_beta = 2.0d0*beta
!
				point%q(2,k) = two_times_beta*u1
				point%q(3,k) = two_times_beta*u2
				point%q(4,k) = two_times_beta*u3				
				point%q(5,k) = -two_times_beta
!
            enddo
!
		end subroutine 
!
!
end module q_variables_mod
