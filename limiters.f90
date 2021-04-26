module limiters_mod


    use data_structure_mod

contains


!	The following subroutine computes the Venkatakrishnan limiter ..
!
!
	subroutine venkat_limiter(qtilde, phi, k)
!
!
		implicit none
!
		integer :: r, k
		real*8 :: qtilde(5), phi(5)
		real*8 :: q, del_neg, del_pos
		real*8 :: max_q, min_q, ds, epsi, num, den, temp 
!
		do r = 1, 5
			q = point%q(r, k) 
			del_neg = qtilde(r) - q
!
			if(dabs(del_neg) .le. 10e-6) then
				phi(r)=1.0d0
			else if(dabs(del_neg) .gt. 10e-6) then              
				if(del_neg .gt. 0.0d0) then 
					del_pos = point%qm(1,r,k) - q
				else if(del_neg .lt. 0.0d0) then
                    del_pos = point%qm(2,r,k) - q
                endif
!
				epsi = VL_CONST*point%min_dist(k)
				epsi = epsi**3.0d0
				num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
				num = num*del_neg + 2.0d0*del_neg*del_neg*del_pos
				den = del_pos*del_pos + 2.0d0*del_neg*del_neg ! Denominator ..
				den = den + del_neg*del_pos + epsi*epsi
				den = den*del_neg
				temp = num/den
!
				if(temp .lt. 1.d0) then
					phi(r) = temp
				else 
					phi(r) = 1.d0
				endif
			endif
        enddo
!
	end subroutine 
!
!
end module limiters_mod
