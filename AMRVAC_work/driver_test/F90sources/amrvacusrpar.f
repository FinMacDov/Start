!##############################################################################
! include amrvacusrpar - nul


INTEGER,PARAMETER:: grav0_=neqpar, grav1_=grav0_+1,grav2_=grav0_&
   +2, BB1_=grav0_+2+1,BB2_=BB1_+1, BB3_=BB2_+1, nspecialpar=2+4

CHARACTER*20,PARAMETER:: specialparname='g1 g2 bb1 bb2 bb3'
INTEGER, PARAMETER:: jmax=500
double precision :: pa(jmax),Temper(jmax),rhoa(jmax),ya(jmax),mua(jmax)
double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit, k_B,miu0,vunit,&
   tunit,punit,heatunit
double precision :: dr,Ti,J_sp,J_d
double precision :: rho1,rho2,rho3,rho4,rho5,rho6,p1,p2,p3,p4,p5,p6


! end include amrvacusrpar - nul
!##############################################################################
COMMON /doub/ pa,Temper,rhoa,ya,mua,Lunit,Teunit,nHunit,runit,Bunit,mHunit,&
    k_B,miu0,vunit,tunit,punit,heatunit,dr,Ti,J_sp,J_d,rho1,rho2,rho3,rho4,&
   rho5,rho6,p1,p2,p3,p4,p5,p6
