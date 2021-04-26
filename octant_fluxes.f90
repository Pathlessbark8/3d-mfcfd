module octant_fluxes_mod
!
!
	use parameter_mod
!
contains
!
!
	subroutine flux_Gwxn(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S1, B1, S3, B3, A1neg, A3neg
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
			
!
!
		beta = 0.5*rho/pr
		S1 = ut1*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A1neg = 0.5*(1 - derf(S1))     
		A3neg = 0.5*(1 - derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut1*A1neg - B1)  
		G(1) = rho*temp1*A3neg 
!      
		G(2) = (pr*A1neg + rho*ut1*temp1)*A3neg
!
		G(3) = rho*ut2*temp1*A3neg
!
		G(4) = rho*temp1*(un*A3neg - B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut1*A1neg - (temp2 + 0.5*pr)*B1
		
		G(5) =  temp3*A3neg - 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Gwxp(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S1, B1, S3, B3, A1pos, A3neg
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
			
!
!
		beta = 0.5*rho/pr
		S1 = ut1*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A1pos = 0.5*(1 + derf(S1))     
		A3neg = 0.5*(1 - derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut1*A1pos + B1)  
		G(1) = rho*temp1*A3neg 
!      
		G(2) = (pr*A1pos + rho*ut1*temp1)*A3neg
!
		G(3) = rho*ut2*temp1*A3neg
!
		G(4) = rho*temp1*(un*A3neg - B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut1*A1pos + (temp2 + 0.5*pr)*B1
		
		G(5) =  temp3*A3neg - 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Goxp(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S1, B1, S3, B3, A1pos, A3pos
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S1 = ut1*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A1pos = 0.5*(1 + derf(S1))     
		A3pos = 0.5*(1 + derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut1*A1pos + B1)  
		G(1) = rho*temp1*A3pos
!      
		G(2) = (pr*A1pos + rho*ut1*temp1)*A3pos
!
		G(3) = rho*ut2*temp1*A3pos
!
		G(4) = rho*temp1*(un*A3pos + B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut1*A1pos + (temp2 + 0.5*pr)*B1
		
		G(5) =  temp3*A3pos + 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Goxn(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S1, B1, S3, B3, A1neg, A3pos
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S1 = ut1*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A1neg = 0.5*(1 - derf(S1))     
		A3pos = 0.5*(1 + derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut1*A1neg - B1)  
		G(1) = rho*temp1*A3pos
!      
		G(2) = (pr*A1neg + rho*ut1*temp1)*A3pos
!
		G(3) = rho*ut2*temp1*A3pos
!
		G(4) = rho*temp1*(un*A3pos + B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut1*A1neg - (temp2 + 0.5*pr)*B1
		
		G(5) =  temp3*A3pos + 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Gwyn(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S2, B2, S3, B3, A2neg, A3neg
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S2 = ut2*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A2neg = 0.5*(1 - derf(S2))     
		A3neg = 0.5*(1 - derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut2*A2neg - B2)  
		G(1) = rho*temp1*A3neg 
!      
		G(2) = rho*ut1*temp1*A3neg

		G(3) = (pr*A2neg + rho*ut2*temp1)*A3neg
!
		G(4) = rho*temp1*(un*A3neg - B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut2*A2neg - (temp2 + 0.5*pr)*B2
		
		G(5) =  temp3*A3neg - 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Gwyp(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S2, B2, S3, B3, A2pos, A3neg
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S2 = ut2*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A2pos = 0.5*(1 + derf(S2))     
		A3neg = 0.5*(1 - derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut2*A2pos + B2)  
		G(1) = rho*temp1*A3neg 
!      
		G(2) = rho*ut1*temp1*A3neg

		G(3) = (pr*A2pos + rho*ut2*temp1)*A3neg
!
		G(4) = rho*temp1*(un*A3neg - B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut2*A2pos + (temp2 + 0.5*pr)*B2
		
		G(5) =  temp3*A3neg - 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Goyp(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S2, B2, S3, B3, A2pos, A3pos
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S2 = ut2*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A2pos = 0.5*(1 + derf(S2))     
		A3pos = 0.5*(1 + derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut2*A2pos + B2)  
		G(1) = rho*temp1*A3pos
!      
		G(2) = rho*ut1*temp1*A3pos

		G(3) = (pr*A2pos + rho*ut2*temp1)*A3pos
!
		G(4) = rho*temp1*(un*A3pos + B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut2*A2pos + (temp2 + 0.5*pr)*B2
		
		G(5) =  temp3*A3pos + 0.5*rho*un*B3*temp1
!
!
      end
!
!
	subroutine flux_Goyn(G, t1, t2, n, prim)
!
!
		implicit none
!
		double precision G(5), prim(5)
		double precision u1, u2, u3, rho, pr
		double precision t1(3), t2(3), n(3)
		double precision ut1, ut2, un
		double precision beta, S2, B2, S3, B3, A2neg, A3pos
		double precision temp1, temp2, temp3, u_sqr
! 
!
		rho = prim(1)
		u1 = prim(2)
		u2 = prim(3)
		u3 = prim(4)
		pr = prim(5)
!
		ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
		ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
		un = n(1)*u1 + n(2)*u2 + n(3)*u3
!
		beta = 0.5*rho/pr
		S2 = ut2*dsqrt(beta) 
		S3 = un*dsqrt(beta) 		
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
		A2neg = 0.5*(1 - derf(S2))     
		A3pos = 0.5*(1 + derf(S3))     		
!
		u_sqr = ut1*ut1 + ut2*ut2 + un*un
!
!     Expressions for the split fluxes ..	
!
		temp1 = (ut2*A2neg - B2)  
		G(1) = rho*temp1*A3pos
!      
		G(2) = rho*ut1*temp1*A3pos

		G(3) = (pr*A2neg + rho*ut2*temp1)*A3pos
!
		G(4) = rho*temp1*(un*A3pos + B3)
!
		temp2 = 2.5d0*pr + 0.5*rho*u_sqr
		
		temp3 = (temp2 + pr)*ut2*A2neg - (temp2 + 0.5*pr)*B2
		
		G(5) =  temp3*A3pos + 0.5*rho*un*B3*temp1
!
!
      end
!
!
end module octant_fluxes_mod
