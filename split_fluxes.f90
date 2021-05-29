module split_fluxes_mod
    !
    !
        use parameter_mod
    !
    !
    !
    contains
    !
    !
        subroutine flux_Gxp(Gxp, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gxp(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S1, B1, A1pos
            double precision temp1, temp2, u_sqr
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
            beta = 0.5d0*rho/pr
            S1 = ut1*dsqrt(beta) 
            B1 = 0.5d0*dexp(-S1*S1)/dsqrt(pi*beta)
            A1pos = 0.5d0*(1.0d0 + derf(S1))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (ut1*A1pos + B1)  
            Gxp(1) = rho*temp1 
    !      

            Gxp(2) = pr*A1pos + rho*ut1*temp1
    !
            Gxp(3) = rho*ut2*temp1
    !
            Gxp(4) = rho*un*temp1
    !
            temp1 = 2.5d0*pr + 0.5d0*rho*u_sqr
    !		temp1 = 2.5d0*pr 
            Gxp(5) =  (temp1 + pr)*ut1*A1pos + (temp1 + 0.5d0*pr)*B1
    !
    !		Gxp(5) = (3.0d0*pr + 0.5d0*rho*u_sqr)*B1 
    !
          end
    !
    !
    !
        subroutine flux_Gxn(Gxn, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gxn(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S1, B1, A1neg
            double precision temp1, temp2, u_sqr
    ! 
    !
            rho = prim(1)
            u1 = prim(2)
            u2 = prim(3)
            u3 = prim(4)
            pr = prim(5)
    !
    !
            ut1 = t1(1)*u1 + t1(2)*u2 + t1(3)*u3
            ut2 = t2(1)*u1 + t2(2)*u2 + t2(3)*u3
            un = n(1)*u1 + n(2)*u2 + n(3)*u3
    !
            beta = 0.5*rho/pr
            S1 = ut1*dsqrt(beta) 
            B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
            A1neg = 0.5*(1 - derf(S1))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (ut1*A1neg - B1)  
            Gxn(1) = rho*temp1 
    !      
            Gxn(2) = pr*A1neg + rho*ut1*temp1
    !
            Gxn(3) = rho*ut2*temp1
    !
            Gxn(4) = rho*un*temp1
    !
            temp1 = 2.5d0*pr + 0.5*rho*u_sqr
            Gxn(5) =  (temp1 + pr)*ut1*A1neg - (temp1 + 0.5*pr)*B1
    !
    !
          end
    !
    !
        subroutine flux_Gyp(Gyp, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gyp(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S2, B2, A2pos
            double precision temp1, temp2, u_sqr
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
            S2 = ut2*dsqrt(beta) 
            B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
            A2pos = 0.5*(1 + derf(S2))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (ut2*A2pos + B2)  
            Gyp(1) = rho*temp1 
    !      
            Gyp(2) = rho*ut1*temp1
    !		
            Gyp(3) = pr*A2pos + rho*ut2*temp1
    !
            Gyp(4) = rho*un*temp1
    !
            temp1 = 2.5d0*pr + 0.5*rho*u_sqr
            Gyp(5) =  (temp1 + pr)*ut2*A2pos + (temp1 + 0.5*pr)*B2
    !
    !
          end
    !
    !
        subroutine flux_Gyn(Gyn, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gyn(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S2, B2, A2neg
            double precision temp1, temp2, u_sqr
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
            S2 = ut2*dsqrt(beta) 
            B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
            A2neg = 0.5*(1 - derf(S2))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (ut2*A2neg - B2)  
            Gyn(1) = rho*temp1 
    !      
            Gyn(2) = rho*ut1*temp1
    !		
            Gyn(3) = pr*A2neg + rho*ut2*temp1
    !
            Gyn(4) = rho*un*temp1
    !
            temp1 = 2.5d0*pr + 0.5*rho*u_sqr
            Gyn(5) =  (temp1 + pr)*ut2*A2neg - (temp1 + 0.5*pr)*B2
    !
    !
          end
    !
    !
        subroutine flux_Gzp(Gzp, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gzp(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S3, B3, A3pos
            double precision temp1, temp2, u_sqr
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
            S3 = un*dsqrt(beta) 
            B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
            A3pos = 0.5*(1 + derf(S3))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (un*A3pos + B3)  
            Gzp(1) = rho*temp1 
    !      
            Gzp(2) = rho*ut1*temp1
    !
            Gzp(3) = rho*ut2*temp1		
    !		
            Gzp(4) = pr*A3pos + rho*un*temp1
    !
    !
            temp1 = 2.5d0*pr + 0.5*rho*u_sqr
            Gzp(5) =  (temp1 + pr)*un*A3pos + (temp1 + 0.5*pr)*B3
    !
    !
          end
    !
    !
        subroutine flux_Gzn(Gzn, t1, t2, n, prim)
    !
    !
            implicit none
    !
            double precision Gzn(5), prim(5)
            double precision u1, u2, u3, rho, pr
            double precision t1(3), t2(3), n(3)
            double precision ut1, ut2, un
            double precision beta, S3, B3, A3neg
            double precision temp1, temp2, u_sqr
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
            S3 = un*dsqrt(beta) 
            B3 = 0.5*dexp(-S3*S3)/dsqrt(pi*beta)
            A3neg = 0.5*(1 - derf(S3))     
    !
            u_sqr = ut1*ut1 + ut2*ut2 + un*un
    !
    !     Expressions for the split fluxes ..	
    !
            temp1 = (un*A3neg - B3)  
            Gzn(1) = rho*temp1 
    !      
            Gzn(2) = rho*ut1*temp1
    !
            Gzn(3) = rho*ut2*temp1		
    !		
            Gzn(4) = pr*A3neg + rho*un*temp1
    !
    !
            temp1 = 2.5d0*pr + 0.5*rho*u_sqr
            Gzn(5) =  (temp1 + pr)*un*A3neg - (temp1 + 0.5*pr)*B3
    !
    !
          end
    !
    !
    !
    end module split_fluxes_mod