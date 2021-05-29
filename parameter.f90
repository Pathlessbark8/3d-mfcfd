module parameter_mod
    !
    !    
            implicit none
    !
    !		integer, parameter :: max_points = 136161
    !		integer, parameter :: max_points = 7623
            integer, parameter :: max_points = 26133
            integer, parameter :: max_iters = 1000
    !
    !		Flow conditions ..
    !
            real*8, parameter :: Mach = 2.0d0
            real*8, parameter :: aoa = 0.0d0
    !
    !		Other parameters ..
    !
            real*8, parameter :: gamma = 1.4d0
            real*8, parameter :: pi=4.0d0*datan(1.0d0)
            real*8, parameter :: theta = aoa*pi/180.0d0
            real*8, parameter :: power = 0.0d0
            real*8, parameter :: VL_CONST = 2.0d0
            real*8, parameter :: CFL = 5.0d0
            integer, parameter :: inner_iterations = 0
    !
    !
    !		Freestream values of the primitive variables ..
    !
            real*8, parameter :: rho_inf = 1.0d0
            real*8, parameter :: u1_inf = Mach*dcos(theta)
            real*8, parameter :: u2_inf = Mach*dsin(theta)
            real*8, parameter :: u3_inf = 0.0d0
            real*8, parameter :: pr_inf = 1.0d0/1.4d0
    !
    !
    end module parameter_mod