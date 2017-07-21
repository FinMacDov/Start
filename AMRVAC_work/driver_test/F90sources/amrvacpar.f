!===============================================================================
! include file amrvacpar.t.mhd

CHARACTER*3,PARAMETER:: typephys='mhd'            ! VACPHYS module name


 

CHARACTER*11,PARAMETER:: eqparname='gamma eta etah etahyper' !Equation parameter names


INTEGER,PARAMETER:: rho_=1,m0_=rho_,m1_=m0_+1,m2_=m0_+2,m3_=m0_+3



INTEGER,PARAMETER:: e_=m3_+1
INTEGER,PARAMETER:: b0_=e_,b1_=b0_+1,b2_=b0_+2,b3_=b0_+3  ! flow variables
INTEGER,PARAMETER:: nwflux=2+2*3+1
INTEGER,PARAMETER:: nwprim=2+2*3


INTEGER,PARAMETER:: Dtr1_=b3_+1


INTEGER,PARAMETER:: nwwave=8


integer,parameter:: nwaux=0
integer,parameter:: nwextra=0
integer,parameter:: nw=nwflux+nwaux+nwextra


integer,parameter:: ee_=e_
integer,parameter:: rhos_=e_
INTEGER, PARAMETER:: p_=e_, pp_=ee_    ! Primitive variables

integer,parameter:: tr1_=Dtr1_

INTEGER,PARAMETER:: v0_=m0_, v1_=m1_,v2_=m2_,v3_=m3_

INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_  ! Polar var. names
INTEGER,PARAMETER:: br_=b0_+r_,bphi_=b0_+phi_,bz_=b0_+z_

integer, parameter :: nvector=2                             ! No. vector vars
integer, dimension(nvector), parameter :: iw_vector=(/ m0_, b0_ /)

INTEGER,PARAMETER:: fastRW_=3,fastLW_=4,slowRW_=5,slowLW_=6 ! Characteristic
INTEGER,PARAMETER:: entroW_=8,diverW_=7,alfvRW_=1,alfvLW_=2 ! waves
INTEGER,PARAMETER:: nworkroe=15



 

INTEGER,PARAMETER:: gamma_=1,eta_=2,etah_=3,etahyper_=4,neqpar=4 !equation params


DOUBLE PRECISION::smalle,minrho,minp
INTEGER,PARAMETER:: nflag_=nw+1
INTEGER:: flags(nflag_)
DOUBLE PRECISION:: wflags(nflag_)

! xprob: problem box; iprob: problem
INTEGER:: iprob
DOUBLE PRECISION:: xprobmin1,xprobmin2,xprobmax1,xprobmax2

! end include file amrvacpar.t.mhd
!===============================================================================
COMMON /INTE/ flags,iprob
COMMON /DOUB/ smalle,minrho,minp,wflags,xprobmin1,xprobmin2,xprobmax1,&
   xprobmax2
