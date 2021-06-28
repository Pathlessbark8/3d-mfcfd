!
module flux_residual_mod
!
!
    use data_structure_mod
!    
    use interior_flux_dGxpos_mod
    use interior_flux_dGxneg_mod
    use interior_flux_dGypos_mod
    use interior_flux_dGyneg_mod
    use interior_flux_dGzpos_mod
    use interior_flux_dGzneg_mod
!
	use wall_flux_dGxpos_mod
	use wall_flux_dGxneg_mod
	use wall_flux_dGypos_mod
	use wall_flux_dGyneg_mod
	use wall_flux_dGzneg_mod
!
	use outer_flux_dGxpos_mod
	use outer_flux_dGxneg_mod
	use outer_flux_dGypos_mod
	use outer_flux_dGyneg_mod
	use outer_flux_dGzpos_mod	   
!	
!
contains
!
!
    subroutine eval_flux_residual()
!
!
        implicit none
!
		integer :: i, k
		real*8 :: Gxp(5), Gyp(5), Gzp(5)
		real*8 :: Gxn(5), Gyn(5), Gzn(5)
!
!
		print *,Gxn
		do i = 1, wall_points
!
			k = wall_points_index(i)
!	
			
            call wall_dGx_pos(Gxp, k) 
			
            ! call wall_dGx_neg(Gxn, k) 
            ! call wall_dGy_pos(Gyp, k)
            ! call wall_dGy_neg(Gyn, k)
            ! call wall_dGz_neg(Gzn, k)            
!
			point%flux_res(:,k) = 2.00*(Gxp + Gxn + Gyp + Gyn + Gzn)
!			
			! print *,point%flux_res(:,k)
		enddo
!
! 		do i = 1, outer_points
! !
! 			k = outer_points_index(i)
! !
! 			call outer_dGx_pos(Gxp, k)
! 			call outer_dGx_neg(Gxn, k) 
! 			call outer_dGy_pos(Gyp, k) 
! 			call outer_dGy_neg(Gyn, k) 
! 			call outer_dGz_pos(Gzp, k) 			
! !
! 			point%flux_res(:,k) = (Gxp + Gxn + Gyp + Gyn + Gzp)
! !
! 		enddo
! !
! 		do i = 1, interior_points
! !
! 			k = interior_points_index(i)
! !
! 			call interior_dGx_pos(Gxp, k) 
! 			call interior_dGx_neg(Gxn, k) 
! 			call interior_dGy_pos(Gyp, k) 
! 			call interior_dGy_neg(Gyn, k) 
! 			call interior_dGz_pos(Gzp, k) 
! 			call interior_dGz_neg(Gzn, k) 			
! !
! 			point%flux_res(:,k) = (Gxp + Gxn + Gyp + Gyn + Gzp + Gzn)
! !
! 		enddo
!
!
	end subroutine
!
!
end module flux_residual_mod
