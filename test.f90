program test
    use split_fluxes_mod
    use octant_fluxes_mod
    implicit none
    double precision G(5), prim(5)
    double precision t1(3), t2(3), n(3)
    integer i
    open(1,file='test_data.dat',status='old')
    open(2,file='fortran_out.dat',status='old')
    do i=1,1000
            read(1,*) G,t1,t2,n,prim
        call flux_Gxp(G,t1,t2,n,prim)
        call flux_Gxn(G,t1,t2,n,prim)
        call flux_Gyn(G,t1,t2,n,prim)
        call flux_Gyp(G,t1,t2,n,prim)
        call flux_Gzn(G,t1,t2,n,prim)
        call flux_Gzp(G,t1,t2,n,prim)
        call flux_Gwxn(G,t1,t2,n,prim)
        call flux_Gwxp(G,t1,t2,n,prim)
        call flux_Gwyn(G,t1,t2,n,prim)
        call flux_Gwyp(G,t1,t2,n,prim)
        call flux_Goxn(G,t1,t2,n,prim)
        call flux_Goxp(G,t1,t2,n,prim)
        call flux_Goyn(G,t1,t2,n,prim)
        call flux_Goyp(G,t1,t2,n,prim)
            write(2,*) G,t1,t2,n,prim
    end do
end program test