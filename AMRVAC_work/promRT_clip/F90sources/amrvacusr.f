!#############################################################################
! module amrvacusr - promRTideal

! this setup to simulate Rayleigh-Taylor dynamics in a solar prominence
! It can be used in 2.5D or 3D, and assumes ideal MHD (with external gravity)
! One can use a tracer, and exploit GLM (or not)
! Related literature is in:
!   'Solar prominences: "double, double ... boil and bubble',
!     R. Keppens, X. Cia, & O. Porth, 2015, ApJ Letters 806, L13 (7pp)

!=============================================================================
subroutine specialsource_impl(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,dtnew,dx1,dx2,x)

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
    dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl
!=============================================================================
subroutine fixp_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x)
include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
!----------------------------------------------------------------------------


end subroutine fixp_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,flag)

include 'amrvacdef.f'

integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
!=============================================================================
module usr_bc
implicit none
double precision, allocatable, save :: pbc(:),rbc(:)
end module usr_bc
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
logical, save :: firstinitglobal=.true.

integer,dimension(:),allocatable:: seed
integer::  seed_size,ix
real:: randphaseP1(1:1000)
!-----------------------------------------------------------------------------
! CGS unit
k_B=1.3806d-16    ! erg*K-1,erg*K-1
miu0=4.d0*dpi     ! Gauss2,Gauss2 cm2,cm2 dyne-1,dyne-1
Lunit=1.d9        ! cm = 10 Mm
UNIT_LENGTH=Lunit
Teunit=1.d6       ! K
nHunit=1.d9       ! cm-3,cm-3
mHunit=1.67262d-24 ! g
runit=1.4d0*mHunit*nHunit     ! 2.341668000000000E-015 g*cm-3,g*cm-3
UNIT_DENSITY=runit
punit=2.3d0*nHunit*k_B*Teunit !0.317538000000000 erg*cm-3,erg*cm-3 !where 1 erg*cm-3,erg*cm-3 =0.1 pascals
Bunit=dsqrt(miu0*punit)       ! 1.99757357615242 Gauss
vunit=Bunit/dsqrt(miu0*runit) ! 1.16448846777562E007 cm/s = 116.45 km/s
UNIT_VELOCITY=vunit
tunit=Lunit/vunit     ! 85.8746159942810 s
heatunit=punit/tunit  ! 3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s

! units for convert !!ie this is deciding how we normilise
if(iprob==-1) then
  normvar(0) = one
else
  normvar(0) = UNIT_LENGTH
endif
normvar(rho_) = UNIT_DENSITY
normvar(v1_)   = UNIT_VELOCITY 
normvar(v2_)   = UNIT_VELOCITY 
normvar(v3_)   = UNIT_VELOCITY 
normvar(pp_)     = UNIT_VELOCITY**2 * UNIT_DENSITY !! I think this is a normilised pressure
normvar(b1_)   = dsqrt(4.0d0*dpi*normvar(pp_)) 
normvar(b2_)   = dsqrt(4.0d0*dpi*normvar(pp_)) 
normvar(b3_)   = dsqrt(4.0d0*dpi*normvar(pp_)) !!this has no units so is a normilized mag feild and is eq 1.9975735761524234 on t=0
normt = UNIT_LENGTH/UNIT_VELOCITY ! t=D/S this is nomrisiled time

eqpar(grav1_)=0.d0
eqpar(grav2_)=-2.74d4*Lunit/vunit**2 !solar surface gravity -2.74d2 m*s-2,m*s-2


eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=0.0d0

! gzone gives the distance from the lower boundary xprobmin2
! to the so-called photosphere where density is fixed to rho0 further on
! this has to be wider than the ghost layer zone for the coarsest mesh....
gzone=0.03d0 ! our grid starts above photosphere
dr=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
! Solar radius
SRadius=6.961d10/Lunit
! base density and temperature
rho0=1.d17/nHunit !! density of photosphere
Tch=8.d3/Teunit !!
Tco=1.8d6/Teunit !!
Tpromin=6.0d3/Teunit !!
Tpromax=1.4d4/Teunit !!
htra1=0.2d0 !! tracer for prom
htra2=1.125d0 !! tracer for
htra3=2.0d0 !!  tracer for
ybot=1.75d0 !!
ytop=2.0d0  !!
bsca=1.5d0  !!
pwidth=0.5d0 !! width in Lunit

   eqpar(nxmodes_)=50 !! number of random modes for the velocity pertabtions.
   eqpar(BB1_)=zero!0.1d0 !!
   eqpar(BB2_)=10.d0/Bunit !!
   eqpar(BB3_)=zero!4.d0 !!
   eqpar(eps_)=0.05d0 !!


eqpar(Cr_)=-0.2d0


!! I think on the first run it goes in here to set up the HDstatic conditions
if(firstinitglobal) then
  call inithdstatic
  firstinitglobal=.false.
endif


!! this add random
randphase(1:1000)=zero
if(eqpar(nxmodes_)>1000) call mpistop('too many modes, edit amrvacusrpar')

if(mype==0)then
    call random_seed(SIZE=seed_size)
    allocate(seed(seed_size))
    call random_seed(GET=seed(1:seed_size))
    call random_number(randphaseP1(1:nint(eqpar(nxmodes_))))
    randphase(1:nint(eqpar(nxmodes_)))=-dpi+two*dpi*dble(randphaseP1(1:nint&
       (eqpar(nxmodes_))))
endif
call MPI_BARRIER(icomm,ierrmpi)
if(npe>1)then
     call MPI_BCAST(randphase,1000,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
endif

if(mype==0)then
    print *,'number of modes=',eqpar(nxmodes_)
    open(123,file='phaseinfo',form='formatted')
    write(123,*) nint(eqpar(nxmodes_))
    do ix=1,nint(eqpar(nxmodes_))
        write(123,"(i4,1es12.4)") ix,randphase(ix)
    enddo
    close(123)
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
use usr_bc
include 'amrvacdef.f'

integer :: j,na,ibc,ix
double precision:: Ta(jmax),gg(jmax),Taext(jmax)
double precision:: res,rhob,pb,wtra1,wtra2,wtra3,b3
!----------------------------------------------------------------------------
wtra1=0.01d0 !! tracer
wtra2=0.01d0 !! tracer
wtra3=0.01d0 !! tracer
!! this is where the temp profile is set
do j=1,jmax
   ya(j)=(dble(j)-0.5d0)*dr-gzone+xprobmin2
!! i think this is for the photosphere
   if(ya(j)<=0.4d0) then
     Ta(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
     Taext(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
   endif
!! this is for between the photosphere and prominance
   if(ya(j)>0.4d0.and.ya(j)<ybot) then
     Ta(j)=Tco-0.5d0*(Tco-Tpromin)*(tanh((ya(j)-htra2)/wtra2)+1.d0)
     Taext(j)=Tco
   endif
!! this is for the prom itself
   if(ya(j)>=ybot.and.ya(j)<ytop) then
     !! T = t_min + delt_T*(0 to 1)
     Ta(j)=Tpromin+(Tpromax-Tpromin)*(ya(j)-ybot)/(ytop-ybot)
     Taext(j)=Tco
   endif
!! this is the temp above the prom
   if(ya(j)>ytop) then
     Ta(j)=Tpromax+0.5d0*(Tco-Tpromax)*(tanh((ya(j)-htra3)/wtra3)+1.d0)
     Taext(j)=Tco
   endif
   gg(j)=eqpar(grav2_)*SRadius**2/(SRadius+ya(j))**2
   !!gg(j)=eqpar(grav2_)
enddo
! solution of hydrostatic equation
ra(1)=rho0
raext(1)=rho0 !! external density
pa(1)=rho0*Tch !! p = rho*T
paext(1)=rho0*Tch !! external pressure p = rho*T
do j=2,jmax
  if(ya(j)<htra2) then
    pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)&
       +gg(j-1))/Ta(j)/4.d0)
    paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0)/(one-dr*(gg(j)&
       +gg(j-1))/Taext(j)/4.d0)
  else
    if(ya(j)<ybot) then
    b3=(eqpar(BB3_)*dexp(-(ya(j-1)-htra2)/bsca))**2+(eqpar(BB3_)*dexp(-(ya(j)&
       -htra2)/bsca))**2
    else
    b3=zero
    endif
    pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0+dr/bsca*b3/2.d0)/(one&
       -dr*(gg(j)+gg(j-1))/Ta(j)/4.d0)
    paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0&
       +dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/Taext(j)/4.d0)
  endif
  ra(j)=pa(j)/Ta(j) !! rho in prom rho = p/T
  raext(j)=paext(j)/Taext(j) !! rho external to prom
end do
! initialized rho and p in the fixed bottom boundary
na=floor(gzone/dr+0.5d0) !! i think these are some sort of step size
res=gzone-(dble(na)-0.5d0)*dr !! maybe resolution
rhob=ra(na)+res/dr*(ra(na+1)-ra(na)) !! maybe rho bc
pb=pa(na)+res/dr*(pa(na+1)-pa(na)) !! pressure bc


allocate(rbc(dixB))
allocate(pbc(dixB))
do ibc=dixB,1,-1 !! dixb = 2, this counts backwards ie 2,1
  na=floor((gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0))/dr+0.5d0)
  res=gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
  rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
  pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
end do

!! this creates txt file pruns with the output shown below
if (mype==0) then
 open(123,file='pruns',form='formatted')
 write(123,*) jmax, 'ya(ix)     ','pa(ix)    ','paext(ix)    ','ra(ix)    ',&
    'raext(ix)    ','Ta(ix)    ','Taext    '
 do ix=1,jmax
    write(123,"(i7,7es12.4)") ix,ya(ix)*Lunit,pa(ix)*punit,paext(ix)&
       *punit,ra(ix)*runit,raext(ix)*runit,Ta(ix),Taext(ix)
 enddo
 close(123)
endif

end subroutine inithdstatic
!=============================================================================
subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

double precision:: psi(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tmp(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision:: res,sigma,lxsize,sigma3
integer :: ix1,ix2,na,idims,imode
logical, save:: first=.true.
logical :: patchw(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
!----------------------------------------------------------------------------

if (first) then
   if (mype==0) then
      print *,'2.5D or 3D MHD Rayleigh-Taylor instability in prominence'
      print *, 'Each time step is: ', tunit, ' s'
   end if
   first=.false.
end if
sigma=0.02d0
sigma3=0.1d0

w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=0.d0
w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=0.d0
w(ixmin1:ixmax1,ixmin2:ixmax2,v3_)=0.d0
! now add the incompressible perturbations
psi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=zero
lxsize=(xprobmax1-xprobmin1)
do imode=1,nint(eqpar(nxmodes_))
   psi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=psi(ixGlo1:ixGhi1,ixGlo2:ixGhi2) &
      +eqpar(eps_)*(dcos(two*dpi*dble(imode)*(x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
      1)/lxsize)+randphase(imode)) /dble(imode)) *dexp(-((x(ixGlo1:ixGhi1,&
      ixGlo2:ixGhi2,2)-htra2)/sigma)**2)
enddo
! compute dpsi/dy
idims=2
select case(typegrad)
    case("central")
     call gradient(psi,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
        ixmax2,idims,tmp)
    case("limited")
     call gradientS(psi,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
        ixmax2,idims,tmp)
end select
w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)&
   -tmp(ixmin1:ixmax1,ixmin2:ixmax2)
! compute dpsi/dx
idims=1
select case(typegrad)
     case("central")
      call gradient(psi,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
         ixmax2,idims,tmp)
     case("limited")
      call gradientS(psi,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
         ixmax2,idims,tmp)
end select
w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)&
   +tmp(ixmin1:ixmax1,ixmin2:ixmax2)


do ix2=ixmin2,ixmax2
do ix1=ixmin1,ixmax1
   na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dr+0.5d0)
   res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dr
   
   
   w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
   w(ix1,ix2,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
  
end do
end do

w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)  =eqpar(BB1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)  =eqpar(BB2_)
w(ixmin1:ixmax1,ixmin2:ixmax2,b3_)  =eqpar(BB3_)
where(x(ixmin1:ixmax1,ixmin2:ixmax2,2)>htra2.and.x(ixmin1:ixmax1,&
   ixmin2:ixmax2,2)<ybot)
   w(ixmin1:ixmax1,ixmin2:ixmax2,b3_)  =eqpar(BB3_)*dexp(-(x(ixmin1:ixmax1,&
      ixmin2:ixmax2,2)-htra2)/bsca)
endwhere
where(x(ixmin1:ixmax1,ixmin2:ixmax2,2)>=ybot)
   w(ixmin1:ixmax1,ixmin2:ixmax2,b3_)  =eqpar(BB3_)*dexp(-(ybot-htra2)/bsca)
endwhere


w(ixmin1:ixmax1,ixmin2:ixmax2,tr1_)=zero
where(x(ixmin1:ixmax1,ixmin2:ixmax2,2)>htra2.and.x(ixmin1:ixmax1,&
   ixmin2:ixmax2,2)<htra3)
 w(ixmin1:ixmax1,ixmin2:ixmax2,tr1_)=one
endwhere
where(x(ixmin1:ixmax1,ixmin2:ixmax2,2)<htra1)
 w(ixmin1:ixmax1,ixmin2:ixmax2,tr1_)=-one
endwhere



w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,psi_)=0.d0


patchw(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=.false.
call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,x,&
   patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iw,iB,w,x)
use usr_bc
! special boundary types, user defined

include 'amrvacdef.f'

integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, iw, iB, ixGmin1,&
   ixGmin2,ixGmax1,ixGmax2
double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

double precision :: dx1,dx2,delydelx,gjjm1,delydelz
double precision :: Teb(ixGlo1:ixGhi1,ixGlo2:ixGhi2),pth(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),cg
double precision:: sigma, mid_pt, r_jet(ixGlo1:ixGhi1,ixGlo2:ixGhi2), jet_w,&
    jet_h, jet_cx, jet_cy !added
logical :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer :: ix1,ix2,idims,ixIntmin1,ixIntmin2,ixIntmax1,ixIntmax2
!----------------------------------------------------------------------------
oktest = index(teststr,'specialbound')>=1
if (oktest) write(unitterm,*) ' === specialbound  (in ) : ', 'ixO^L : ',&
   ixOmin1,ixOmin2,ixOmax1,ixOmax2

select case(iB)
case(3)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)=eqpar(BB1_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)=eqpar(BB3_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)=eqpar(BB2_)

   
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,tr1_)=w(ixOmin1:ixOmax1,ixOmax2+1,Dtr1_)&
        /w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
   enddo
  
   


   
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmax2+1,psi_)
   enddo
  
   

   
   ! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_) &
       +delydelx*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_))
     enddo
   enddo
  
   
   
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)=-w(ixOmin1:ixOmax1,ixOmax2&
      +dixB:ixOmax2+1:-1,m1_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)=-w(ixOmin1:ixOmax1,ixOmax2&
      +dixB:ixOmax2+1:-1,m2_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)=-w(ixOmin1:ixOmax1,ixOmax2&
      +dixB:ixOmax2+1:-1,m3_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
     w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
   enddo
   !This part I added, hopefully it works:
   !!This is where we add the jet
   !For FWHM
   sigma = 0.2d0
   jet_w = 0.05d0
   jet_h = 0.05d0
   jet_cx = (jet_w-jet_w)/2.0d0 !center pts x
   jet_cy = (jet_h-0.0d0)/2.0d0 !center pts y
   r_jet(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (x(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)-jet_cx)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)&
      -jet_cy)**2

   do ix2=ixOmin2,ixOmax2
   do ix1=ixOmin1,ixOmax1
      if (x(ix1,ix2,2).le. jet_h .and. x(ix1,ix2,1).le.jet_w&
         /2.0d0 .and. x(ix1,ix2,1).ge.-jet_w/2.0d0) then
         w(ix1,ix2,rho_) = ra(1)
 !w(ixO1,ixO2,p_) = pa(1)!10.0d0*pa(1)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
         w(ix1,ix2,v2_)  = (2.0d6/vunit)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
         w(ix1,ix2,tr1_) = 100.0d0
      endif
   end do
   end do

  
   

   patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
   call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,w,x,patchw)
case(4)
!! implementation of hydrostatic extrapolation at top boundary
   
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ix2-2,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ix2-1,b1_:b3_))
   enddo
  
   
   ! 2nd order CD for divB=0 to set normal B component better
   
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_) &
       -delydelx*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_))
     enddo
   enddo
  
   

   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,v1_)=w(ixOmin1:ixOmax1,ixOmin2-1,m1_)&
        /w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
     w(ixOmin1:ixOmax1,ix2,v2_)=w(ixOmin1:ixOmax1,ixOmin2-1,m2_)&
        /w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
     w(ixOmin1:ixOmax1,ix2,v3_)=w(ixOmin1:ixOmax1,ixOmin2-1,m3_)&
        /w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
   enddo
   !! obtain the thermal pressure in top layers and store in pth
   ixIntmin1=ixOmin1;ixIntmin2=ixOmin2;ixIntmax1=ixOmax1;ixIntmax2=ixOmax2;
   ixIntmin2=ixOmin2-dixB;ixIntmax2=ixOmin2-1;
   call getpthermal(w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixIntmin1,ixIntmin2,&
      ixIntmax1,ixIntmax2,pth)
   cg=dxlevel(2)*eqpar(grav2_)/two
   !! fill pth, rho ghost layers according to gravity stratification
   
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1,ixOmax1
       Teb(ix1,ix2-1)=pth(ix1,ix2-1)/w(ix1,ix2-1,rho_)
       Teb(ix1,ix2)=Teb(ix1,ix2-1)
       gjjm1=half*(SRadius**2/(SRadius+x(ix1,ix2,2))**2+SRadius**2/(SRadius&
          +x(ix1,ix2-1,2))**2)
       pth(ix1,ix2)=(pth(ix1,ix2-1)+cg*gjjm1*w(ix1,ix2-1,rho_))/(one&
          -cg*gjjm1/Teb(ix1,ix2))
       if(pth(ix1,ix2)<minp) pth(ix1,ix2)=pth(ix1,ix2-1)
       w(ix1,ix2,rho_)=pth(ix1,ix2)/Teb(ix1,ix2)
     end do
   end do
  
   
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,tr1_)=w(ixOmin1:ixOmax1,ixOmin2-1,Dtr1_)&
        /w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
   enddo
  
   


   
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmin2-1,psi_)
   enddo
  
   

   patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
   call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,w,x,patchw)
case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision :: bQgrid(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer :: iw
!-----------------------------------------------------------------------------

call addsource_gravSA(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,x)

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(out) :: ggrid(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!---------------------------------------------------------------------------
! calculate gravity
ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=eqpar(grav2_)*(SRadius/(SRadius&
   +x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2

end subroutine getggrav
!=============================================================================
subroutine addsource_gravSA(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
! gravity distribution along a magnetic loop (a circular arc)

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision :: ggrid(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer :: iw, idims
!---------------------------------------------------------------------------

call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,x)

! add sources from gravity
do iw= iwmin,iwmax
   select case (iw)
   case (m1_,m2_)
     ! dm_i/dt= +rho*g_i
      idims=iw-m0_
      if (abs(eqpar(grav0_+idims))>smalldouble) w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,m0_+idims)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
         +idims) +qdt*ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
         *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
   case (e_)
     ! de/dt= +g_i*m_i
      do idims=1,ndim
         if (abs(eqpar(grav0_+idims))>smalldouble) w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ee_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ee_) &
            +qdt*ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,m0_+idims)
      end do
   end select
end do

end subroutine addsource_gravSA
!=============================================================================
subroutine getdt_special(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,dtnew,dx1,dx2,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
    dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

call getdt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,&
   dtnew,dx1,dx2,x)

end subroutine getdt_special
!=============================================================================
subroutine getdt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,dtnew,dx1,dx2,x)

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   1:ndim), w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
double precision, intent(inout) :: dtnew

double precision:: dxinv(1:ndim), dtgrav
integer:: idims
!----------------------------------------------------------------------------

dxinv(1)=one/dx1;dxinv(2)=one/dx2;
dtgrav=bigdouble
do idims=1,ndim
   if(abs(eqpar(grav0_+idims))>zero)dtgrav=min(dtgrav,one/sqrt(abs(eqpar&
      (grav0_+idims))*dxinv(idims)))
enddo

dtnew=dtgrav

end subroutine getdt_grav
!=============================================================================
subroutine specialeta(w,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,&
   ixmax2, idirmin
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,7-2&
   *ndir:3), eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
   ixmin1,ixmin2,ixmax1,ixmax2,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
   ixmin2,ixmax1,ixmax2
double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
    x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iflag
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(out):: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop&
   (' iflag> nw, make change in parfile or in user file')

var(ixImin1:ixImax1,ixImin2:ixImax2) = zero

end subroutine specialvarforerrest
!=============================================================================
subroutine process_grid_usr(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
   ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in):: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw&
   +nwauxio),tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tmp_alt(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision                   :: normconv(0:nw+nwauxio)

double precision :: wloc(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
!-----------------------------------------------------------------------------

wloc(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)
if(saveprim)then
   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      p_)
 else
   call getpthermal(wloc,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,tmp)
endif

if(iprob==-1)then
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=wloc(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,rho_)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=tmp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)/wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
else
  ! output the temperature p/rho
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=tmp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)/wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)*Teunit
  ! output the plasma beta p*2/B**2
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=tmp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*two/(wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2&
     +wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2+wloc(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,b3_)**2)

!stuff I added
  !w(ixO^S,nw+3)=tmp(ixO^S)
  !w(ixO^S,nw+4)=wloc(ixO^S,rho_)
  !w(ixO^S,nw+5)= Teunit
  !w(ixO^S,nw+6) = (w(ixO^S,rho_)-wloc(ixO^S,rho_))/wloc(ixO^S,rho_)

! output the temperature p/rho
  !w(ixO^S,nw+6)=tmp_alt(ixO^S)/wloc(ixO^S,rho_)*Teunit
  ! output the plasma beta p*2/B**2
  !w(ixO^S,nw+7)=tmp_alt(ixO^S)*two/(^C&wloc(ixO^S,b^C_)**2+)
  !w(ixO^S,nw+8)=tmp_alt(ixO^S)
endif

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

if(iprob==-1)then
  primnames= TRIM(primnames)//' '//'rho T'
  wnames=TRIM(wnames)//' '//'rho T'
else
  primnames= TRIM(primnames)//' '//'T beta'
  wnames=TRIM(wnames)//' '//'T beta test'
endif

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
!-----------------------------------------------------------------------------

call mpistop("special B0 undefined")

end subroutine specialset_B0
!=============================================================================
subroutine bc_int(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g.
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)

end subroutine bc_int
!=============================================================================
subroutine printlog_special

! printlog: calculates volume averaged mean values

include 'amrvacdef.f'

logical :: fileopen
integer :: iigrid, igrid, level, nleafs_level(1:nlevelshi), iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2), volumeflat&
   (1:nlevelshi)
double precision :: tmpw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer :: numlevels, imon, idirmin
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+2+nlevelshi) :: dsum_send, dsum_recv
double precision :: wmeanmore
double precision :: wmaxte,wminte,wmaxvel,wmaxte_mype,wminte_mype,wmaxvel_mype
double precision :: invgminone
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------
volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero
nleafs_level(1:mxnest)=0

wmaxte=zero
wmaxte_mype=zero
wmaxvel=zero
wmaxvel_mype=zero
wminte=1.d20
wminte_mype=1.d20
invgminone=one/(eqpar(gamma_)-one)
wmeanmore=zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   nleafs_level(level)=nleafs_level(level)+1
   volumeflat(level)=volumeflat(level)+ (rnode(rpxmax1_,igrid)&
      -rnode(rpxmin1_,igrid))*(rnode(rpxmax2_,igrid)-rnode(rpxmin2_,igrid))
   if (slab) then
      dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=rnode(rpdx1_,igrid)&
         *rnode(rpdx2_,igrid)
   else
      dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=pgeo(igrid)%dvolume(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2)
      volume(level)=volume(level)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2))
   end if
   ! set dxlevel for use in gradient evaluation
   dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
   ! just use array for velocity
   tmpw(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dsqrt((pw(igrid)%w(ixMlo1:ixMhi1,&
      ixMlo2:ixMhi2,m1_)/pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))**2 &
      +(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,m2_)/pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))**2 +(pw(igrid)%w(ixMlo1:ixMhi1,&
      ixMlo2:ixMhi2,m3_)/pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))**2)
   wmaxvel_mype=max(wmaxvel_mype,maxval(tmpw(ixMlo1:ixMhi1,ixMlo2:ixMhi2)))

   wmean(Dtr1_)=wmean(Dtr1_)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,Dtr1_)/pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))


   wmean(psi_)=wmean(psi_)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,psi_))

   wmean(rho_)=wmean(rho_)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))
   ! kinetic energy in x
   wmean(m1_)=wmean(m1_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,m1_)**2)/pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))
   ! kinetic energy in y
   wmean(m2_)=wmean(m2_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,m2_)**2)/pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))
   ! kinetic energy in z
   wmean(m3_)=wmean(m3_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,m3_)**2)/pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_))
   ! total energy
   wmean(e_)=wmean(e_)+sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,e_))
   ! magnetic energy
   wmean(b1_)=wmean(b1_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      * (pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b1_)**2+pw(igrid)%w&
      (ixMlo1:ixMhi1,ixMlo2:ixMhi2,b2_)**2+pw(igrid)%w(ixMlo1:ixMhi1,&
      ixMlo2:ixMhi2,b3_)**2))
   ! magnetic energy in y
   wmean(b2_)=wmean(b2_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b2_)**2))
   ! magnetic energy in z
   wmean(b3_)=wmean(b3_)+half*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b3_)**2))
   call getpthermal(pw(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
      ixMlo1,ixMlo2,ixMhi1,ixMhi2,tmpw)
   wmeanmore=wmeanmore+invgminone*sum(dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      *tmpw(ixMlo1:ixMhi1,ixMlo2:ixMhi2))
   ! maximal temperature
   wmaxte_mype=max(wmaxte_mype,maxval(tmpw(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      /pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_)))
   ! minimal temperature
   wminte_mype=min(wminte_mype,minval(tmpw(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
      /pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,rho_)))
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

call MPI_REDUCE(wmaxte_mype,wmaxte,1,MPI_DOUBLE_PRECISION, MPI_MAX,0,icomm,&
   ierrmpi)
call MPI_REDUCE(wmaxvel_mype,wmaxvel,1,MPI_DOUBLE_PRECISION, MPI_MAX,0,icomm,&
   ierrmpi)

call MPI_REDUCE(wminte_mype,wminte,1,MPI_DOUBLE_PRECISION, MPI_MIN,0,icomm,&
   ierrmpi)

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
dsum_send(nw+2+numlevels:nw+2+numlevels)=wmeanmore
call MPI_REDUCE(dsum_send,dsum_recv,nw+2+numlevels,MPI_DOUBLE_PRECISION,&
    MPI_SUM,0,icomm,ierrmpi)
isum_send(1:numlevels)=nleafs_level(levmin:levmax)
call MPI_REDUCE(isum_send,isum_recv,numlevels,MPI_INTEGER, MPI_SUM,0,icomm,&
   ierrmpi)

if (mype==0) then

   wmean(1:nw)=dsum_recv(1:nw)
   wmeanmore=dsum_recv(nw+2+numlevels)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)
   nleafs_level(levmin:levmax)=isum_recv(1:numlevels)

   wmean=wmean/voltotal
   wmeanmore=wmeanmore/voltotal

   ! determine coverage in coordinate space
   volprob=(xprobmax1-xprobmin1)*(xprobmax2-xprobmin2)
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, MPI_INFO_NULL,log_fh,&
         ierrmpi)
      opened=.true.
      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout),&
          MPI_CHARACTER,status,ierrmpi)
      !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnameslog)-1
      write(wnameslog(i+3:i+14),"(a11)") "d1 d2 d3 d4"
      i=i+17
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "c",level
      end do
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "n",level
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a79)')"it   t  dt res ",wnameslog
         else
           write(line,'(a15,a79)')"it   t   dt    ",wnameslog
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a79)')"it res ",wnameslog
         else
           write(line,'(a7,a79)')"it     ",wnameslog
         endif
      end if

      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, status,&
         ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(e13.5))')it,t,dt,residual
      else
         write(line,'(i7,2(e13.5))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(e13.5))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)
   do iw=1,nw
      write(line,'(e13.5)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do
   do imon=1,1
      write(line,'(e13.5)')wmeanmore
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do
   write(line,'(e13.5)')wmaxte
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)

   write(line,'(e13.5)')wminte
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)
   write(line,'(e13.5)')wmaxvel
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
      ierrmpi)
   do level=1,mxnest
      write(line,'(e13.5)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(i6)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), MPI_CHARACTER,status,&
         ierrmpi)
   end do

end if

end subroutine printlog_special
!=============================================================================
subroutine userspecialconvert(qunitconvert)

include 'amrvacdef.f'
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type

integer  :: iigrid,igrid
logical  :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

if(mype==0)then
   print *,'converting to primitives, no specialvarout for integrals'
endif

do iigrid=1,igridstail; igrid=igrids(iigrid)
  !!call primitive(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x)
  call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     pw(igrid)%w,px(igrid)%x)
end do

call spatial_integral_w


patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=.false.
do iigrid=1,igridstail; igrid=igrids(iigrid)
  !!call conserve(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x,patchw)
  call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     pw(igrid)%w,px(igrid)%x,patchw)
end do

end subroutine userspecialconvert

!=============================================================================
subroutine mask_gridtrpos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,patchwi)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
logical, intent(inout)             :: patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision :: trtreshold
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
trtreshold=0.95d0
patchwi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   tr1_)>trtreshold)

return
end subroutine mask_gridtrpos
!=============================================================================
subroutine mask_gridtrzer(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,patchwi)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
logical, intent(inout)             :: patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision :: trtreshold,trtreshold2
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
trtreshold=0.05d0
trtreshold2=-0.05d0
patchwi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   tr1_)<trtreshold.and.w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_)>trtreshold2)

return
end subroutine mask_gridtrzer
!=============================================================================

subroutine spatial_integral_w

include 'amrvacdef.f'

double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2), timephy,xmom,ymom
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid

integer           :: nregions,ireg
integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,region
logical           :: patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2),alive
!-----------------------------------------------------------------------------
nregions=1


nregions=3

do ireg=1,nregions
 select case(ireg)
 case(1)
   region='fulldomain'

 case(2)
   region='trpos'
 case(3)
   region='trzer'

 end select
! number of integrals to perform
ni=12
allocate(integral_ipe(ni),integral_w(ni))
integral_ipe=0.d0
integral_w=0.d0

do iigrid=1,igridstail; igrid=igrids(iigrid);
  if(slab) then
    dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=rnode(rpdx1_,igrid)&
       *rnode(rpdx2_,igrid)
  else
    dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=pgeo(igrid)%dvolume(ixMlo1:ixMhi1,&
       ixMlo2:ixMhi2)
  end if
  dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
    myB0_face1 => pB0_face1(igrid)
    myB0_face2 => pB0_face2(igrid)
  end if
  typelimiter=typelimiter1(node(plevel_,igrid))
  typegradlimiter=typegradlimiter1(node(plevel_,igrid))
  patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=.false.
  select case(region)

  case('trpos')
     call mask_gridtrpos(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
        ixMhi2,pw(igrid)%w,px(igrid)%x,patchwi)
  case('trzer')
     call mask_gridtrzer(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
        ixMhi2,pw(igrid)%w,px(igrid)%x,patchwi)

  case('fulldomain')
     patchwi(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=.true.
  case default
     call mpistop("region not defined")
  end select
  integral_ipe(1)=integral_ipe(1)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,1,patchwi)
  integral_ipe(2)=integral_ipe(2)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,2,patchwi)
  integral_ipe(3)=integral_ipe(3)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,3,patchwi)
  integral_ipe(4)=integral_ipe(4)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,4,patchwi)
  integral_ipe(5)=integral_ipe(5)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,5,patchwi)
  integral_ipe(6)=integral_ipe(6)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,6,patchwi)
  integral_ipe(7)=integral_ipe(7)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,7,patchwi)
  integral_ipe(8)=integral_ipe(8)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,8,patchwi)
  integral_ipe(9)=integral_ipe(9)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
     ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,9,patchwi)
  integral_ipe(10)=integral_ipe(10)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,&
     ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,10,&
     patchwi)
  integral_ipe(11)=integral_ipe(11)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,&
     ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,11,&
     patchwi)
  integral_ipe(12)=integral_ipe(12)+ integral_grid(ixGlo1,ixGlo2,ixGhi1,&
     ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,pw(igrid)%w,px(igrid)%x,dvolume,12,&
     patchwi)
end do
call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
!!timephy=t*tunit
timephy=t
if(mype==0) then
  !!print *,'t=',timephy,'integral_w=',integral_w
  write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(region),".int"
  inquire(file=filename,exist=alive)
  if(alive) then
    open(unit=21,file=filename,form='formatted',status='old',access='append')
  else
    open(unit=21,file=filename,form='formatted',status='new')
  endif
  write(21,'(13(es12.4))') timephy,integral_w(1:12)
  close(21)
endif

deallocate(integral_ipe,integral_w)

enddo
end subroutine spatial_integral_w
!=============================================================================
function integral_grid(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,dvolume,intval,patchwi)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,intval
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim),dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
logical, intent(in) :: patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir) :: bvec,&
   current
double precision :: integral_grid,gammm1
integer :: ix1,ix2,idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------
gammm1=eqpar(gamma_)-one

integral_grid=0.d0
select case(intval)
 case(1)
  ! volume integration
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+dvolume(ix1,ix2)
  end do
  end do
 case(2)
  ! mass integration
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+w(ix1,ix2,rho_)&
        *dvolume(ix1,ix2)
  end do
  end do
 case(3)
  ! vertical velocities
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+dabs(w(ix1,ix2,v2_))&
        *dvolume(ix1,ix2)
  end do
  end do
 case(4)
  ! horizontal velocities
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+dsqrt(w(ix1,ix2,v1_)**2&
        +w(ix1,ix2,v3_)**2)*dvolume(ix1,ix2)
  end do
  end do
 case(5)
  ! T integration
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+(w(ix1,ix2,p_)&
        /w(ix1,ix2,rho_))*dvolume(ix1,ix2)
  end do
  end do
 case(6)
  ! magnetic energy integration
  if(B0field) then
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,1)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,2)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,3);
  else
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_);
  endif
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+ half*(bvec(ix1,ix2,1)&
        **2+bvec(ix1,ix2,2)**2+bvec(ix1,ix2,3)**2)*dvolume(ix1,ix2)
  end do
  end do
 case(7)
  ! kinetic energy integration
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+ half*w(ix1,ix2,rho_)&
        *(w(ix1,ix2,v0_+1)**2+w(ix1,ix2,v0_+2)**2+w(ix1,ix2,v0_&
        +3)**2)*dvolume(ix1,ix2)
  end do
  end do
 case(8)
  ! plasma beta integration
  if(B0field) then
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,1)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,2)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,3);
  else
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_);
  endif
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+ (two*w(ix1,ix2,p_)&
        /(bvec(ix1,ix2,1)**2+bvec(ix1,ix2,2)**2+bvec(ix1,ix2,3)&
        **2))*dvolume(ix1,ix2)
  end do
  end do
 case(9)
  ! J^2  integration (split case as well, current from B1 only))
  bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
     b1_)
  bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
     b2_)
  bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
     b3_);
  call curlvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,current,idirmin,1,ndir)
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+ (current(ix1,ix2,1)**2&
        +current(ix1,ix2,2)**2+current(ix1,ix2,3)**2)*dvolume(ix1,ix2)
  end do
  end do
 case(10)
  ! internal energy integration
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+(w(ix1,ix2,p_)&
        /gammm1)*dvolume(ix1,ix2)
  end do
  end do
 case(11)
  ! vertical flow with sign
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+w(ix1,ix2,v2_)&
        *dvolume(ix1,ix2)
  end do
  end do
 case(12)
  ! vertical B component
  if(B0field) then
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,1)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,2)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_)+myB0_cell%w(ixImin1:ixImax1,ixImin2:ixImax2,3);
  else
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b1_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b2_)
    bvec(ixImin1:ixImax1,ixImin2:ixImax2,3)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
       b3_);
  endif
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
     if(patchwi(ix1,ix2)) integral_grid=integral_grid+dabs(bvec(ix1,ix2,2))&
        *dvolume(ix1,ix2)
  end do
  end do
 case default
     call mpistop("intval not defined")
end select

return
end function integral_grid
!=============================================================================
! amrvacusr.t.promRTideal
!=============================================================================
