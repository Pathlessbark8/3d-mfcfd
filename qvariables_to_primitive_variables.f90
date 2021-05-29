module qvariables_to_primitive_variables_mod
!
!       
        use data_structure_mod
!
contains
!
!
		subroutine qtilde_to_primitive(q, prim)
!
        implicit none
!
        real*8 :: q(5), prim(5)
        real*8 :: beta, temp, temp1, temp2
!
!
        beta = -q(5)*0.50
!
        temp = 0.50/beta
!
        prim(2) = q(2)*temp
        prim(3) = q(3)*temp
	prim(4) = q(4)*temp
!
        temp1 = q(1) + beta*(prim(2)*prim(2) + prim(3)*prim(3) + prim(4)*prim(4))
        temp2 = temp1 - (dlog(beta)*2.50)
!
		prim(1) = exp(temp2)
		prim(5) = prim(1)*temp

    end subroutine
!
!
!
end module qvariables_to_primitive_variables_mod
