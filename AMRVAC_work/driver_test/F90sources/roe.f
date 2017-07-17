!##############################################################################
! module vacphys.mhdroe - subroutines for Roe-type Riemann solver for MHD
!=============================================================================
subroutine average(wL,wR,x,ixmin1,ixmin2,ixmax1,ixmax2,idim,wroe,workroe)

! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
! Calculate the wroe average of primitive variables in wL and wR, assignment:
! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
!
! wL,wR,wroe are all interface centered quantities

include 'amrvacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim,iw
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nworkroe):: workroe
!-----------------------------------------------------------------------------
call average2(wL,wR,x,ixmin1,ixmin2,ixmax1,ixmax2,idim,wroe,&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   2), workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,3),workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,4),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,5),&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,6), workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,7),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,8))

end subroutine average
!=============================================================================
subroutine average2(wL,wR,x,ixmin1,ixmin2,ixmax1,ixmax2,idim,wroe,cfast,cslow,&
   afast,aslow,csound2,dp, rhodv,tmp)

! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
! Calculate the wroe average of primitive variables in wL and wR, assignment:
! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
!
! wL,wR,wroe are all interface centered quantities

include 'amrvacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim,idir,jdir,iw
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: cfast,cslow,afast,&
   aslow,csound2,dp, rhodv,tmp
!!common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
!-----------------------------------------------------------------------------

if (ndir==1) call mpistop("MHD with d=11 is the same as HD")

oktest=index(teststr,'average')>=1
if(oktest)write(*,*)'Average wL,wR:',wL(ixtest1,ixtest2,iwtest),wR(ixtest1,&
   ixtest2,iwtest)

!Averaging primitive variables
wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))

wroe(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   m1_)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,&
   m1_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
wroe(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   b1_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,b1_))


wroe(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   m2_)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,&
   m2_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
wroe(ixmin1:ixmax1,ixmin2:ixmax2,b2_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   b2_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,b2_))


wroe(ixmin1:ixmax1,ixmin2:ixmax2,v3_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   m3_)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,&
   m3_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
wroe(ixmin1:ixmax1,ixmin2:ixmax2,b3_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,&
   b3_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,b3_))


! Use afast and aslow for pressures pL and pR
call getpthermal(wL,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
   afast)
call getpthermal(wR,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
   aslow)

wroe(ixmin1:ixmax1,ixmin2:ixmax2,pp_)=half*(afast(ixmin1:ixmax1,&
   ixmin2:ixmax2)+aslow(ixmin1:ixmax1,ixmin2:ixmax2))

if(oktest)write(*,*)'Calculate saved variables'



if(useprimitive.or.eqpar(gamma_)<=zero)then
   ! dp=pR-pL
   dp(ixmin1:ixmax1,ixmin2:ixmax2)=aslow(ixmin1:ixmax1,ixmin2:ixmax2)&
      -afast(ixmin1:ixmax1,ixmin2:ixmax2)
else
   !CONSERVATIVE dp=(g-1)*(de-v*dm+0.5*v**2*drho-0.5*d(B**2))
   dp(ixmin1:ixmax1,ixmin2:ixmax2)=(eqpar(gamma_)-one)*(wR(ixmin1:ixmax1,&
      ixmin2:ixmax2,e_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,e_)&
      -(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m1_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,&
         m1_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m1_))+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m2_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,m2_)&
         -wL(ixmin1:ixmax1,ixmin2:ixmax2,m2_))+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m3_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,m3_)&
         -wL(ixmin1:ixmax1,ixmin2:ixmax2,m3_)))&
      +half*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m1_)**2+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m2_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,m3_)&
         **2)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))&
      -half*(wR(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2-wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,b1_)**2+wR(ixmin1:ixmax1,ixmin2:ixmax2,b2_)**2&
         -wL(ixmin1:ixmax1,ixmin2:ixmax2,b2_)**2+wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,b3_)**2-wL(ixmin1:ixmax1,ixmin2:ixmax2,b3_)**2))
endif

!CONSERVATIVE rho*dv_idim=dm_idim-v_idim*drho
rhodv(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
   +idim)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)-wroe(ixmin1:ixmax1,&
   ixmin2:ixmax2,m0_+idim)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
   -wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))

!Calculate csound2,cfast,cslow,alphafast and alphaslow

! get csound**2
call getcsound2prim(wroe,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,&
   ixmax2,csound2)

! aa=B**2/rho+a**2
cfast(ixmin1:ixmax1,ixmin2:ixmax2)=( wroe(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2&
   +wroe(ixmin1:ixmax1,ixmin2:ixmax2,b2_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
   b3_)**2 )/wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+csound2(ixmin1:ixmax1,&
   ixmin2:ixmax2)

! cs**2=0.5*(aa+dsqrt(aa**2-4*a**2*(b_i**2/rho)))
cslow(ixmin1:ixmax1,ixmin2:ixmax2)=half*(cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
   -dsqrt(cfast(ixmin1:ixmax1,ixmin2:ixmax2)**2-4d0*csound2(ixmin1:ixmax1,&
   ixmin2:ixmax2)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)**2&
   /wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)))

! cf**2=aa-cs**2
cfast(ixmin1:ixmax1,ixmin2:ixmax2)=cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
   -cslow(ixmin1:ixmax1,ixmin2:ixmax2)

! alpha_f**2=(a**2-cs**2)/(cf**2-cs**2)
afast(ixmin1:ixmax1,ixmin2:ixmax2)=(csound2(ixmin1:ixmax1,ixmin2:ixmax2)&
   -cslow(ixmin1:ixmax1,ixmin2:ixmax2))/(cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
   -cslow(ixmin1:ixmax1,ixmin2:ixmax2))
afast(ixmin1:ixmax1,ixmin2:ixmax2)=min(one,max(afast(ixmin1:ixmax1,&
   ixmin2:ixmax2),zero))

! alpha_s=dsqrt(1-alpha_f**2)
aslow(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(one-afast(ixmin1:ixmax1,&
   ixmin2:ixmax2))

! alpha_f=dsqrt(alpha_f**2)
afast(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(afast(ixmin1:ixmax1,ixmin2:ixmax2))

! cf=dsqrt(cf**2)
cfast(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(cfast(ixmin1:ixmax1,ixmin2:ixmax2))

! cs=dsqrt(cs**2)
cslow(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(cslow(ixmin1:ixmax1,ixmin2:ixmax2))

if(oktest)write(*,*)'Average:rho,csound2,dp,rhodv',wroe(ixtest1,ixtest2,rho_),&
   csound2(ixtest1,ixtest2),dp(ixtest1,ixtest2),rhodv(ixtest1,ixtest2)
if(oktest)write(*,*)'Average:cf,cs,af,as',cfast(ixtest1,ixtest2),&
   cslow(ixtest1,ixtest2),afast(ixtest1,ixtest2),aslow(ixtest1,ixtest2)

!Replace the primitive variables with more useful quantities:
! rho -> dsqrt(rho)
wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=dsqrt(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
   rho_))

! Avoid sgn(b_idim)==0
where(dabs(wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))<smalldouble)wroe&
   (ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)=smalldouble
! B_idir,jdir -> beta_idir,jdir
idir=idim+1-ndir*(idim/ndir)
if(ndir==2)then
    where(wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)>=zero)
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)=one
    elsewhere
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)=-one
    end where
else
    !beta_j=B_j/dsqrt(B_i**2+B_j**2); beta_i=B_i/dsqrt(B_i**2+B_j**2)
    jdir=idir+1-ndir*(idir/ndir)
    tmp(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
       b0_+idir)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)**2)
    where(tmp(ixmin1:ixmax1,ixmin2:ixmax2)>smalldouble)
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,b0_+idir)/tmp(ixmin1:ixmax1,ixmin2:ixmax2)
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,b0_+jdir)/tmp(ixmin1:ixmax1,ixmin2:ixmax2)
    elsewhere
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)=dsqrt(half)
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)=dsqrt(half)
       !!wroe(ix^S,b0_+idir)=zero
       !!wroe(ix^S,b0_+jdir)=zero
    end where
endif

end subroutine average2
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ixmin1,ixmin2,ixmax1,ixmax2,il,idim,&
   smalla,a,jump,workroe)

! Calculate the il-th characteristic speed and the jump in the il-th
! characteristic variable in the idim direction within ixL.
! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
! variables. However part of the summation is done in advance and saved into
! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
! used in the entropy fix.
!
! All the variables are centered on the cell interface, thus the 
! "*C" notation is omitted for sake of brevity.

include 'amrvacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,il,idim
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: smalla,a,jump
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nworkroe) :: workroe
!!common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
!!!save bdv,bdb
!!!save cs2L,cs2R,cs2ca2L,cs2ca2R
!-----------------------------------------------------------------------------
call geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmax1,ixmax2,il,idim,smalla,a,&
   jump, workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1),workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,2), workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,3),&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,4),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   5),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,6), workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,7),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,8),&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,9),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   10), workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,11),workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,12),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,13))

end subroutine geteigenjump
!=============================================================================
subroutine geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmax1,ixmax2,il,idim,&
   smalla,a,jump, cfast,cslow,afast,aslow,csound2,dp,rhodv,bdv,bdb,cs2L,cs2R,&
   cs2ca2L,cs2ca2R)

! Calculate the il-th characteristic speed and the jump in the il-th
! characteristic variable in the idim direction within ixL.
! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
! variables. However part of the summation is done in advance and saved into
! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
! used in the entropy fix.
!
! All the variables are centered on the cell interface, thus the 
! "*C" notation is omitted for sake of brevity.

include 'amrvacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2,il,idim,idir,jdir
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw):: wL,wR,wroe
double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: smalla,a,jump
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   ::cfast,cslow,&
   afast,aslow,csound2,dp,rhodv
!!common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: bdv,bdb
!!!save bdv,bdb
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)   :: aL,aR,cs2L,cs2R,&
   cs2ca2L,cs2ca2R
!!!save cs2L,cs2R,cs2ca2L,cs2ca2R
!-----------------------------------------------------------------------------

oktest=index(teststr,'geteigenjump')>=1

idir=idim+1-ndir*(idim/ndir)
jdir=idir+1-ndir*(idir/ndir)

if(il==fastRW_)then
   !Fast and slow waves use bdv=sqrho**2*sign(bx)*(betay*dvy+betaz*dvz)
   !                        bdb=sqrho*a*          (betay*dBy+betaz*dBz)
   bdv(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
      +idir)* (wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)/wR(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +idir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
   if(ndir==3)bdv(ixmin1:ixmax1,ixmin2:ixmax2)=bdv(ixmin1:ixmax1,&
      ixmin2:ixmax2)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)&
      * (wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+jdir)/wR(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +jdir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
   bdv(ixmin1:ixmax1,ixmin2:ixmax2)=bdv(ixmin1:ixmax1,ixmin2:ixmax2)&
      *sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**2,wroe(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+idim))

   bdb(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
      +idir)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)-wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+idir))
   if(ndir==3)bdb(ixmin1:ixmax1,ixmin2:ixmax2)=bdb(ixmin1:ixmax1,&
      ixmin2:ixmax2)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)&
      *(wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)-wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+jdir))
   bdb(ixmin1:ixmax1,ixmin2:ixmax2)=bdb(ixmin1:ixmax1,ixmin2:ixmax2)&
      *dsqrt(csound2(ixmin1:ixmax1,ixmin2:ixmax2))*wroe(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)
   if(oktest)write(*,*)'rhobetadv,sqrhoabetadb:',bdv(ixtest1,ixtest2),&
      bdb(ixtest1,ixtest2)
endif

if(il==alfvRW_)then
   !Alfven waves use      bdv=0.5*sqrho**2*      (betaz*dvy-betay*dvz)
   !                      bdb=0.5*sqrho*sign(bx)*(betaz*dBy-betay*dBz)
   bdv(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
      +jdir)* (wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)/wR(ixmin1:ixmax1,&
      ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +idir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)) -wroe(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+idir)* (wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +jdir)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,m0_+jdir)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
   bdb(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
      +jdir)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)-wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+idir)) -wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
      +idir)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)-wL(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+jdir))
   bdv(ixmin1:ixmax1,ixmin2:ixmax2)=bdv(ixmin1:ixmax1,ixmin2:ixmax2)*half&
      *wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)**2
   bdb(ixmin1:ixmax1,ixmin2:ixmax2)=bdb(ixmin1:ixmax1,ixmin2:ixmax2)*half&
      *sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_),wroe(ixmin1:ixmax1,&
      ixmin2:ixmax2,b0_+idim))
   if(oktest)write(*,*)'rhobetaXdv/2,sqrhobetaXdb/2:',bdv(ixtest1,ixtest2),&
      bdb(ixtest1,ixtest2)
endif

select case(il)
   case(fastRW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)+cfast(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half/csound2(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(afast(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         +cfast(ixmin1:ixmax1,ixmin2:ixmax2)*rhodv(ixmin1:ixmax1,&
         ixmin2:ixmax2)+dp(ixmin1:ixmax1,ixmin2:ixmax2))+aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(-cslow(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bdv(ixmin1:ixmax1,ixmin2:ixmax2)+bdb(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(fastLW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-cfast(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half/csound2(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(afast(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         -cfast(ixmin1:ixmax1,ixmin2:ixmax2)*rhodv(ixmin1:ixmax1,&
         ixmin2:ixmax2)+dp(ixmin1:ixmax1,ixmin2:ixmax2))+aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(+cslow(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bdv(ixmin1:ixmax1,ixmin2:ixmax2)+bdb(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(slowRW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)+cslow(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half/csound2(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(aslow(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         +cslow(ixmin1:ixmax1,ixmin2:ixmax2)*rhodv(ixmin1:ixmax1,&
         ixmin2:ixmax2)+dp(ixmin1:ixmax1,ixmin2:ixmax2))+afast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(+cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bdv(ixmin1:ixmax1,ixmin2:ixmax2)-bdb(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(slowLW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-cslow(ixmin1:ixmax1,ixmin2:ixmax2)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=half/csound2(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(aslow(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         -cslow(ixmin1:ixmax1,ixmin2:ixmax2)*rhodv(ixmin1:ixmax1,&
         ixmin2:ixmax2)+dp(ixmin1:ixmax1,ixmin2:ixmax2))+afast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*(-cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bdv(ixmin1:ixmax1,ixmin2:ixmax2)-bdb(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(entroW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)&
         -wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)-dp(ixmin1:ixmax1,&
         ixmin2:ixmax2)/csound2(ixmin1:ixmax1,ixmin2:ixmax2)
   case(diverW_)
      if(divbwave)then
         a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
            +idim)
         jump(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
            +idim)-wL(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)
      else
         a(ixmin1:ixmax1,ixmin2:ixmax2)=zero
         jump(ixmin1:ixmax1,ixmin2:ixmax2)=zero
      endif
   case(alfvRW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)+dabs(wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))&
         /wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=+bdv(ixmin1:ixmax1,ixmin2:ixmax2)&
         -bdb(ixmin1:ixmax1,ixmin2:ixmax2)
   case(alfvLW_)
      a(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
         +idim)-dabs(wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))&
         /wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
      jump(ixmin1:ixmax1,ixmin2:ixmax2)=-bdv(ixmin1:ixmax1,ixmin2:ixmax2)&
         -bdb(ixmin1:ixmax1,ixmin2:ixmax2)
end select

! Calculate "smalla" or modify "a" based on the "typeentropy" switch

select case(typeentropy(il))
case('yee')
   ! Based on Yee JCP 68,151 eq 3.23
   smalla(ixmin1:ixmax1,ixmin2:ixmax2)=entropycoef(il)
case('harten','powell', 'ratio')
   ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
   ! Initialize left and right eigenvalues by velocities
   aL(ixmin1:ixmax1,ixmin2:ixmax2)= wL(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +idim)/wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
   aR(ixmin1:ixmax1,ixmin2:ixmax2)= wR(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
      +idim)/wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
   ! Calculate the final "aL" and "aR"
   select case(il)
   case(fastRW_)
      ! These quantities will be used for all the fast and slow waves
      ! Calculate soundspeed**2 and cs**2+ca**2.
      call getcsound2(wL,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,&
         ixmax2,cs2L)
      call getcsound2(wR,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,&
         ixmax2,cs2R)
      cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2)=cs2L(ixmin1:ixmax1,ixmin2:ixmax2)&
         +(wL(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2+wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,b2_)**2+wL(ixmin1:ixmax1,ixmin2:ixmax2,b3_)**2)&
         /wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
      cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2)=cs2R(ixmin1:ixmax1,ixmin2:ixmax2)&
         +(wR(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2+wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,b2_)**2+wR(ixmin1:ixmax1,ixmin2:ixmax2,b3_)**2)&
         /wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
      ! Save the discriminants into cs2L and cs2R
      cs2L(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(cs2ca2L(ixmin1:ixmax1,&
         ixmin2:ixmax2)**2-4d0*cs2L(ixmin1:ixmax1,ixmin2:ixmax2)&
         *wL(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)**2/wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))
      cs2R(ixmin1:ixmax1,ixmin2:ixmax2)=dsqrt(cs2ca2R(ixmin1:ixmax1,&
         ixmin2:ixmax2)**2-4d0*cs2R(ixmin1:ixmax1,ixmin2:ixmax2)&
         *wR(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)**2/wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,rho_))

      ! The left and right eigenvalues for the fast wave going to right
      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         + dsqrt(half*(cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2L(ixmin1:ixmax1,ixmin2:ixmax2)))
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         + dsqrt(half*(cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2R(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(fastLW_)
      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         - dsqrt(half*(cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2L(ixmin1:ixmax1,ixmin2:ixmax2)))
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         - dsqrt(half*(cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2R(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(slowRW_)
      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         + dsqrt(half*(cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2L(ixmin1:ixmax1,ixmin2:ixmax2)))
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         + dsqrt(half*(cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2R(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(slowLW_)
      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         - dsqrt(half*(cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2L(ixmin1:ixmax1,ixmin2:ixmax2)))
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         - dsqrt(half*(cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2R(ixmin1:ixmax1,ixmin2:ixmax2)))
   case(entroW_,diverW_)
      ! These propagate by the velocity
   case(alfvRW_)
      ! Store the Alfven speeds into cs2ca2L and cs2ca2R
      cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(wL(ixmin1:ixmax1,&
         ixmin2:ixmax2,b0_+idim))/dsqrt(wL(ixmin1:ixmax1,ixmin2:ixmax2,rho_))
      cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(wR(ixmin1:ixmax1,&
         ixmin2:ixmax2,b0_+idim))/dsqrt(wR(ixmin1:ixmax1,ixmin2:ixmax2,rho_))

      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2)
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         + cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2)
   case(alfvLW_)
      aL(ixmin1:ixmax1,ixmin2:ixmax2)=aL(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2)
      aR(ixmin1:ixmax1,ixmin2:ixmax2)=aR(ixmin1:ixmax1,ixmin2:ixmax2) &
         - cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2)
   end select
end select

call entropyfix(ixmin1,ixmin2,ixmax1,ixmax2,il,aL,aR,a,smalla)

end subroutine geteigenjump2
!=============================================================================
subroutine rtimes(q,wroe,ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim
double precision:: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: q,rq
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nworkroe):: workroe
!!common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
!!save:: bv,v2a2
!-----------------------------------------------------------------------------
call rtimes2(q,wroe,ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,rq,&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   2), workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,3),workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,4),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,5),&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,6), workroe(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2,7),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,14),&
   workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,15))

end subroutine rtimes
!=============================================================================
subroutine rtimes2(q,wroe,ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,rq, cfast,&
   cslow,afast,aslow,csound2,dp,rhodv,bv,v2a2)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

include 'amrvacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2,iw,il,idim,idir,jdir
double precision:: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: q,rq
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: cfast,cslow,afast,&
   aslow,csound2,dp,rhodv
!!common /roe/ cfast,cslow,afast,aslow,csound2,dp,rhodv
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: bv,v2a2
!!save:: bv,v2a2
!-----------------------------------------------------------------------------

oktest=index(teststr,'rtimes')>=1

idir=idim+1-ndir*(idim/ndir)
jdir=idir+1-ndir*(idir/ndir)

select case(iw)
  case(rho_)
    select case(il)
      case(fastRW_,fastLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
           *afast(ixmin1:ixmax1,ixmin2:ixmax2)
      case(slowRW_,slowLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
           *aslow(ixmin1:ixmax1,ixmin2:ixmax2)
      case(entroW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)
      case(diverW_,alfvRW_,alfvLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
    end select

 case(e_)
   if(il==fastRW_)then
     if(eqpar(gamma_)>zero)then
!!!IDEAL GAS
       ! Store 0.5*v**2+(2-gamma)/(gamma-1)*a**2
       v2a2(ixmin1:ixmax1,ixmin2:ixmax2)=half*(wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,m1_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,m2_)**2&
          +wroe(ixmin1:ixmax1,ixmin2:ixmax2,m3_)**2)+ &
         (two-eqpar(gamma_))/(eqpar(gamma_)-one)*csound2(ixmin1:ixmax1,&
            ixmin2:ixmax2)
     else
       call mpistop("Correct rTimes for NONIDEAL gas in src/mhd/roe.t")
         ! Express the partial derivative de/dp using wroe
 !v2a2(ixmin1:ixmax1,ixmin2:ixmax2)=half*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m1_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,m2_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,m3_)**2)+ &
 !(??dedp(ixmin1:ixmax1,ixmin2:ixmax2)??-1)*csound2(ixmin1:ixmax1,ixmin2:ixmax2)
     endif
     ! Store sgn(bx)*(betay*vy+betaz*vz) in bv
     bv(ixmin1:ixmax1,ixmin2:ixmax2)=wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
        +idir)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)
     if(ndir==3)bv(ixmin1:ixmax1,ixmin2:ixmax2)=bv(ixmin1:ixmax1,&
        ixmin2:ixmax2)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
        +jdir)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+jdir)
     bv(ixmin1:ixmax1,ixmin2:ixmax2)=bv(ixmin1:ixmax1,ixmin2:ixmax2)&
        *sign(one,wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))
     if(oktest)write(*,*)'v2/2+(2-g)/(g-1)a2,betav:',&
          v2a2(ixtest1,ixtest2),bv(ixtest1,ixtest2)
   else if(il==alfvRW_)then
     !Store betaz*vy-betay*vz in bv
     bv(ixmin1:ixmax1,ixmin2:ixmax2)=(wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
        +jdir)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idir)-&
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,m0_+jdir))
     if(oktest)write(*,*)'betaXv:',bv(ixtest1,ixtest2)
   endif
   select case(il)
   case(fastRW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         -aslow(ixmin1:ixmax1,ixmin2:ixmax2)*cslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*bv(ixmin1:ixmax1,ixmin2:ixmax2)+afast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*&
           (v2a2(ixmin1:ixmax1,ixmin2:ixmax2)+cfast(ixmin1:ixmax1,&
              ixmin2:ixmax2)*(cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
              +wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))))
   case(fastLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         +aslow(ixmin1:ixmax1,ixmin2:ixmax2)*cslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*bv(ixmin1:ixmax1,ixmin2:ixmax2)+afast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*&
           (v2a2(ixmin1:ixmax1,ixmin2:ixmax2)+cfast(ixmin1:ixmax1,&
              ixmin2:ixmax2)*(cfast(ixmin1:ixmax1,ixmin2:ixmax2)&
              -wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))))
   case(slowRW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         +afast(ixmin1:ixmax1,ixmin2:ixmax2)*cfast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*bv(ixmin1:ixmax1,ixmin2:ixmax2)+aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*&
           (v2a2(ixmin1:ixmax1,ixmin2:ixmax2)+cslow(ixmin1:ixmax1,&
              ixmin2:ixmax2)*(cslow(ixmin1:ixmax1,ixmin2:ixmax2)&
              +wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))))
   case(slowLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)*(&
         -afast(ixmin1:ixmax1,ixmin2:ixmax2)*cfast(ixmin1:ixmax1,&
         ixmin2:ixmax2)*bv(ixmin1:ixmax1,ixmin2:ixmax2)+aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2)*&
              (v2a2(ixmin1:ixmax1,ixmin2:ixmax2)+cslow(ixmin1:ixmax1,&
                 ixmin2:ixmax2)*(cslow(ixmin1:ixmax1,ixmin2:ixmax2)&
                 -wroe(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim))))
   case(entroW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)= q(ixmin1:ixmax1,ixmin2:ixmax2)*half&
         *(wroe(ixmin1:ixmax1,ixmin2:ixmax2,m1_)**2+wroe(ixmin1:ixmax1,&
         ixmin2:ixmax2,m2_)**2+wroe(ixmin1:ixmax1,ixmin2:ixmax2,m3_)**2)
   case(diverW_)
      if(divbwave)then
         rq(ixmin1:ixmax1,ixmin2:ixmax2)= q(ixmin1:ixmax1,ixmin2:ixmax2)&
            *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)
      else
         rq(ixmin1:ixmax1,ixmin2:ixmax2)= zero
      endif
   case(alfvRW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=+q(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bv(ixmin1:ixmax1,ixmin2:ixmax2)
   case(alfvLW_)
      rq(ixmin1:ixmax1,ixmin2:ixmax2)=-q(ixmin1:ixmax1,ixmin2:ixmax2)&
         *bv(ixmin1:ixmax1,ixmin2:ixmax2)
   end select

  case(m1_,m2_,m3_)
    if(iw==m0_+idim)then
      select case(il)
        case(fastRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *afast(ixmin1:ixmax1,ixmin2:ixmax2)*(wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)+cfast(ixmin1:ixmax1,ixmin2:ixmax2))
        case(fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *afast(ixmin1:ixmax1,ixmin2:ixmax2)*(wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)-cfast(ixmin1:ixmax1,ixmin2:ixmax2))
        case(slowRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *aslow(ixmin1:ixmax1,ixmin2:ixmax2)*(wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)+cslow(ixmin1:ixmax1,ixmin2:ixmax2))
        case(slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *aslow(ixmin1:ixmax1,ixmin2:ixmax2)*(wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)-cslow(ixmin1:ixmax1,ixmin2:ixmax2))
        case(entroW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *wroe(ixmin1:ixmax1,ixmin2:ixmax2,iw)
        case(diverW_,alfvLW_,alfvRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
      end select
    else
      select case(il)
        case(fastRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *(afast(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)-aslow(ixmin1:ixmax1,ixmin2:ixmax2)&
             *cslow(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_-m0_+iw)*sign(one,wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_+idim)))
        case(fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *(afast(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)+aslow(ixmin1:ixmax1,ixmin2:ixmax2)&
             *cslow(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_-m0_+iw)*sign(one,wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_+idim)))
        case(slowRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *(aslow(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)+afast(ixmin1:ixmax1,ixmin2:ixmax2)&
             *cfast(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_-m0_+iw)*sign(one,wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_+idim)))
        case(slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *(aslow(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,iw)-afast(ixmin1:ixmax1,ixmin2:ixmax2)&
             *cfast(ixmin1:ixmax1,ixmin2:ixmax2)*wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_-m0_+iw)*sign(one,wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,b0_+idim)))
        case(entroW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *wroe(ixmin1:ixmax1,ixmin2:ixmax2,iw)
        case(diverW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
        case(alfvRW_)
          if(iw==m0_+idir)then
            rq(ixmin1:ixmax1,ixmin2:ixmax2)=+q(ixmin1:ixmax1,ixmin2:ixmax2)&
               *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)
          else
            rq(ixmin1:ixmax1,ixmin2:ixmax2)=-q(ixmin1:ixmax1,ixmin2:ixmax2)&
               *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)
          endif
        case(alfvLW_)
          if(iw==m0_+idir)then
            rq(ixmin1:ixmax1,ixmin2:ixmax2)=-q(ixmin1:ixmax1,ixmin2:ixmax2)&
               *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)
          else
            rq(ixmin1:ixmax1,ixmin2:ixmax2)=+q(ixmin1:ixmax1,ixmin2:ixmax2)&
               *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)
          endif        
      end select
    end if ! iw=m_idir,m_jdir
  case(b1_,b2_,b3_)
    if(iw==b0_+idim)then
      if(il==diverW_ .and. divbwave)then
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)
      else
        rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
      endif
    else
      select case(il)
        case(fastRW_,fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=+q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *aslow(ixmin1:ixmax1,ixmin2:ixmax2)*dsqrt(csound2(ixmin1:ixmax1,&
             ixmin2:ixmax2))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,iw)&
             /wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
        case(slowRW_,slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=-q(ixmin1:ixmax1,ixmin2:ixmax2)&
             *afast(ixmin1:ixmax1,ixmin2:ixmax2)*dsqrt(csound2(ixmin1:ixmax1,&
             ixmin2:ixmax2))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,iw)&
             /wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_)
        case(entroW_,diverW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2)=zero
        case(alfvRW_,alfvLW_)
          if(iw==b0_+idir)then
             rq(ixmin1:ixmax1,ixmin2:ixmax2)=-q(ixmin1:ixmax1,ixmin2:ixmax2)&
                *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+jdir)&
                /sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_),&
                wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))
          else
             rq(ixmin1:ixmax1,ixmin2:ixmax2)=+q(ixmin1:ixmax1,ixmin2:ixmax2)&
                *wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idir)&
                /sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,rho_),&
                wroe(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim))
          end if
      end select
    end if ! iw=b_idir,b_jdir
end select

end subroutine rtimes2
!=============================================================================
! end module vacphys.mhdroe
!##############################################################################
