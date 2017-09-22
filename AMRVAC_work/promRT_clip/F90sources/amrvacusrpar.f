!##############################################################################
! include amrvacusrpar - promRTideal

INTEGER,PARAMETER:: grav0_=neqpar, grav1_=grav0_+1,grav2_=grav0_&
   +2, eps_=grav0_+2+1, nxmodes_=eps_+1, BB1_=nxmodes_+1,BB2_=BB1_&
   +1,BB3_=BB2_+1,nspecialpar=2+5

CHARACTER*21,PARAMETER:: specialparname='g1 g2 eps nxm bb1 bb3'



INTEGER, PARAMETER:: jmax=80000
double precision :: pa(jmax),ra(jmax),ya(jmax),raext(jmax),paext(jmax)
double precision :: Lunit,Teunit,nHunit,runit,Bunit,mHunit,k_B,miu0,vunit,&
   tunit,punit,heatunit
double precision :: SRadius,dr,gzone,rho0,Tch,Tco,bsca,htra1,htra2,htra3,&
    ybot,ytop,Tpromax,Tpromin,pwidth,Tcoext

DOUBLE PRECISION:: randphase(1000)
! end include amrvacusrpar - promRTideal
!##############################################################################
COMMON /DOUB/ randphase
COMMON /doub/ pa,ra,ya,raext,paext,Lunit,Teunit,nHunit,runit,Bunit,mHunit,k_B,&
   miu0,vunit,tunit,punit,heatunit,SRadius,dr,gzone,rho0,Tch,Tco,bsca,htra1,&
   htra2,htra3, ybot,ytop,Tpromax,Tpromin,pwidth,Tcoext
