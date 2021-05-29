module conserved_vector_Ubar_mod
!
!
	use data_structure_mod
!
	contains
!
		subroutine conserved_vector_Ubar(k, Ubar)
!
		implicit none
		!
		!
			real*8	::	u1_inf_rot, u2_inf_rot, u3_inf_rot, rho_e_inf
			real*8	::	u1, u2, u3, pr, rho, u1_rot, u2_rot, u3_rot, rho_e
			real*8	::	beta, S3, B3_inf, A3n_inf
			real*8	::	B3, A3p, temp1, temp2
			real*8	:: Ubar(5)
			integer	:: k
!			
!
!			
			u1_inf_rot = u1_inf*point%tan1(1,k) + u2_inf*point%tan1(2,k) + u3_inf*point%tan1(3,k)
			u2_inf_rot = u1_inf*point%tan2(1,k) + u2_inf*point%tan2(2,k) + u3_inf*point%tan2(3,k)
			u3_inf_rot = u1_inf*point%nor(1,k) + u2_inf*point%nor(2,k) + u3_inf*point%nor(3,k)
!			
			rho_e_inf = pr_inf*2.50 + 0.5*rho_inf*(u1_inf_rot*u1_inf_rot + u2_inf_rot*u2_inf_rot + u3_inf_rot*u3_inf_rot)
!
			beta = (0.5*rho_inf)/pr_inf
			S3 = u3_inf_rot*sqrt(beta)
			B3_inf = exp(-S3*S3)/(2*sqrt(pi*beta))
			A3n_inf = 0.5*(1-derf(S3))
!	
			rho = point%prim(1,k)
			u1 = point%prim(2,k)
			u2 = point%prim(3,k)
			u3 = point%prim(4,k)
			pr = point%prim(5,k)
!			
			u1_rot = u1*point%tan1(1,k) + u2*point%tan1(2,k) + u3*point%tan1(3,k)
			u2_rot = u1*point%tan2(1,k) + u2*point%tan2(2,k) + u3*point%tan2(3,k)
			u3_rot = u1*point%nor(1,k) + u2*point%nor(2,k) + u3*point%nor(3,k)
!
			rho_e = pr*2.50 + 0.5*rho*(u1_rot*u1_rot + u2_rot*u2_rot + u3_rot*u3_rot)
!
			beta = (0.5*rho)/pr
			S3 = u3_rot*sqrt(beta)
			B3 = exp(-S3*S3)/(2*sqrt(pi*beta))
			A3p = 0.5*(1+derf(S3))
!
			Ubar(1) = (rho_inf*A3n_inf) + (rho*A3p)
			Ubar(2) = (rho_inf*u1_inf_rot*A3n_inf) + (rho*u1_rot*A3p)
			Ubar(3) = (rho_inf*u2_inf_rot*A3n_inf) + (rho*u2_rot*A3p)
!
			temp1 = rho_inf*(u3_inf_rot*A3n_inf - B3_inf)
			temp2 = rho*(u3_rot*A3p + B3)
			Ubar(4) = temp1 + temp2
!
			temp1 = (rho_e_inf*A3n_inf - 0.5*rho_inf*u3_inf_rot*B3_inf)
			temp2 = (rho_e*A3p + 0.5*rho*u3_rot*B3)
!
			Ubar(5) = temp1 + temp2
!
!
		end subroutine
		!
		!	
!
end module conserved_vector_Ubar_mod
!	
