!##############################################################################
! include amrvacusrpar - nul


INTEGER,PARAMETER:: grav0_=neqpar, grav^D_=grav0_+^D, BB1_=grav0_+^ND+1,BB2_=BB1_+1, BB3_=BB2_+1, nspecialpar=^ND+4

CHARACTER*20,PARAMETER:: specialparname='g1 g2 bb1 bb2 bb3'
INTEGER, PARAMETER:: jmax=500
COMMON, double precision :: pa(jmax),Temper(jmax),rhoa(jmax),ya(jmax),mua(jmax)
COMMON, double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit, &
                            k_B,miu0,vunit,tunit,punit,heatunit
COMMON, double precision :: dr,Ti,J_sp,J_d
COMMON, double precision :: rho1,rho2,rho3,rho4,rho5,rho6,p1,p2,p3,p4,p5,p6


! end include amrvacusrpar - nul
!##############################################################################
