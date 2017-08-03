!=============================================================================
!=============================================================================
! INCLUDE:amrvacnul/specialini.t
!=============================================================================
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
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
!integer, intent(in)                :: ixI^LSUB1,ixO^LSUB1
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw&
   +nwauxio),tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision                   :: normconv(0:nw+nwauxio)


!logical patchw(ixG^T)
double precision :: wloc(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
!-----------------------------------------------------------------------------


wloc(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)
tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)
if(saveprim)then
   call getpthermal(wloc,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,tmp)
 else
   call getpthermal(wloc,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,tmp)
endif
! output the plasma beta p*2/B**2
 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
    *two/(wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+wloc(ixOmin1:ixOmax1,&
    ixOmin2:ixOmax2,b2_)**2+wloc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2)

!patchw(ixG^S)=.false.
!call conserve(ixG^L,ixO^L,w,x,patchw)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

!call mpistop("special wnames and primnames undefined")

! Example : as above in specialvar_output, assuming relativistic HD here...
primnames= TRIM(primnames)//' '//'beta'
wnames=TRIM(wnames)//' '//'beta'

end subroutine specialvarnames_output
!=============================================================================



!============================================================================
subroutine addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

integer :: iw, idims
!-----------------------------------------------------------------------------
! add sources from gravity
do iw= iwmin,iwmax
   select case (iw)
   case (m1_,m2_)
      ! dm_i/dt= +rho*g_i
      idims=iw-m0_
      if (abs(eqpar(grav0_+idims))>smalldouble) w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,m0_+idims)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
         +idims) +qdt*eqpar(grav0_+idims)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)

   case (e_)
      ! de/dt= +g_i*m_i
      do idims=1,ndim
         if (abs(eqpar(grav0_+idims))>smalldouble) &
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ee_)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ee_) &
              +qdt*eqpar(grav0_+idims)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 m0_+idims)
      end do

   end select
end do        

end subroutine addsource_grav
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
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/simple/specialsource.t
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
subroutine correctaux_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
integer, intent(inout)         :: patchierror(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

! correct solution from analytic case

end subroutine correctaux_usr
!==========================================================================================
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
logical, save :: firstinitglobal=.true.
!-----------------------------------------------------------------------------
!where eqpar are set
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=0.0d0 !this gives ideal MHD

! normilastion in terms of SI unit
k_B=1.3806d-23         ! J.K-1
m_p = 1.672621777e-27  ! kg
R_gas =k_B/m_p         ! J.k-1.kg-1
miu0=1.257d-6          ! H.m-1
Lunit=1.d6             ! m
UNIT_LENGTH=Lunit      ! m
Teunit=1.d5            ! K
nHunit=1.d15           ! m-3
mHunit=1.67262d-27     ! kg
runit= 1.4d0*mHunit*nHunit ! kg.m-3
UNIT_DENSITY=runit
punit=  2.3d0*nHunit*k_B*Teunit ! m-3.J.K-1.K = mb1.kg.sb2 = pa
Bunit= dsqrt(miu0*punit) !sqrt(H.m-2.kg.s-2) = sqrt(kg2.s-4.A-2)= kg.A.s-2 = T
vunit=  Bunit/dsqrt(miu0*runit) ! m/s
UNIT_VELOCITY=vunit
tunit=Lunit/vunit ! s
heatunit=punit/tunit ! mb1.kg.sb3 = pa.s-1
Ti = tunit! s

! units for convert
if(iprob==-1) then
  normvar(0) = one
else
  normvar(0) = UNIT_LENGTH
endif
normvar(rho_) = UNIT_DENSITY
normvar(v1_)   = UNIT_VELOCITY 
normvar(v2_)   = UNIT_VELOCITY 
normvar(v3_)   = UNIT_VELOCITY 
normvar(pp_)     = UNIT_VELOCITY**2 * UNIT_DENSITY
normvar(b1_)   = dsqrt(4.0d0*dpi*normvar(pp_)) 
normvar(b2_)   = dsqrt(4.0d0*dpi*normvar(pp_)) 
normvar(b3_)   = dsqrt(4.0d0*dpi*normvar(pp_)) 
normt = UNIT_LENGTH/UNIT_VELOCITY

eqpar(grav1_)=0.d0
eqpar(grav2_)= -275.4229*Lunit/vunit**2 !where [Lunit/vunit**2] = [s2/m],[s2/m]

!dr=(xprobmax2-xprobmin2)/dble(jmax) ! step size
dr=(xprobmax2-xprobmin2)/dble(jmax-6) ! step size

eqpar(BB1_)=0.d0 !Bx
eqpar(BB2_)=0.004/Bunit !By = 100G=0.01T, 40G=0.004T
eqpar(BB3_)=0.d0 !Bz

J_sp = 1.5d4/vunit !25 km s-1 (average speed of a spicule)

!this set up intail condtions
if(firstinitglobal) then
  call inithdstatic
  firstinitglobal=.false.
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
include 'amrvacdef.f'

real, dimension(jmax) :: rho, p, mu, Tem, z
integer :: i,ix,j,na,k
double precision:: res
!----------------------------------------------------------------------------
open (unit = 1, file ="data_256.0/data_p.dat", status='old')
open (unit = 2, file ="data_256.0/data_mu.dat", status='old')
open (unit = 3, file ="data_256.0/data_rho.dat", status='old')
open (unit = 4, file ="data_256.0/data_T.dat", status='old')
open (unit = 5, file ="data_256.0/data_Z.dat", status='old')

do i=1,jmax
 read(1,*) p(i) !dyn/cm2,!dyn/cm2
 read(2,*) mu(i) !dimension mean molecular weight
 read(3,*) rho(i) !g/cm3,!g/cm3
 read(4,*) Tem(i) !k
 read(5,*) Z(i) !Mm
end do

!Need to covert CGS to SI and then make dimensionless

do j=1,jmax
   Tea(j)=Tem(j)/Teunit
   rhoa(j)=rho(j)*1000.d0/runit
   pa(j)=0.1d0*p(j)/punit
   mua(j) = mu(j) !This is already dimensionless. This mu is the average mol wieght
   ya(j) = Z(j)*1000000/Lunit
enddo

!equation ideal gas law (igl): p = R_gas*rho*T
do k=1,jmax
   pigl(k) = ((R_gas/mu(k))*rho(k)*Tem(k))/punit
enddo

!k = jmax
iniene = pigl((k-1)-2)/(eqpar(gamma_)-1.0d0)-(eqpar(BB1_)*eqpar(BB1_)&
   +eqpar(BB2_)*eqpar(BB2_)+eqpar(BB3_)*eqpar(BB3_))/2.0d0
!print *, pigl(k-3)*punit

!ya(j-1)*Lunit this gives last element of matrix

J_d = rhoa(1) !jet density

! this creates txt file pruns with the output shown below
if (mype==0) then
 open(123,file='output_test',form='formatted')
 write(123,*) jmax, '  ya(ix)              ','pa(ix)                   ',&
    'pidgl              ','Temp(ix)                  ',&
    'rho(ix)                       ', 'mu(ix)           '
 do ix=1,jmax
    write(123,*) ya(ix), pa(ix)*punit, pigl(ix)*punit, Tea(ix)&
       *Teunit, rhoa(ix)*runit, mua(ix)
 enddo
 close(123)
endif
close(1)
close(2)
close(3)
close(4)
close(5)

end subroutine inithdstatic
!=============================================================================
subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
   ixmax1,ixmax2,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,ixmax1,&
   ixmax2
integer :: ix1,ix2,na,imode, jx
double precision:: res, lxsize, sigma, phase, mid_pt, eps, dy
double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
double precision:: rinlet(ixGlo1:ixGhi1,ixGlo2:ixGhi2), r_jet(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), p_bg(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision:: psi(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision :: jet_w, jet_h, jet_cx, jet_cy

logical, save:: first=.true., first_1=.true.
logical patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------



if (first) then
   if (mype==0) then
      !print *, 'ixmin2:', ixmin2, 'ixmax2:', ixmax2
      print *,'2.5D MHD jet simulation'
      print *, 'B0:', eqpar(BB2_)*Bunit, 'T'
      print *, 'Time unit:', Ti, ' s'
      print *, 'normilised Va speed: ', vunit, ' m s-1'
      print *, 'Jet Speed is: ', J_sp*vunit, ' m s-1'
      print *, 'Jet rho is: ', J_d*runit, ' kg m^2/s^2'
      print *, 'Rho0 = ', rhoa(1)*runit, 'g m^2/s^2'
      print *, 'Rhoz = ', rhoa(size(rhoa))*runit, 'g m^2/s^2'
   end if
   first=.false.
end if

w(ixmin1:ixmax1,ixmin2:ixmax2,v1_)=0.d0
w(ixmin1:ixmax1,ixmin2:ixmax2,v2_)=0.d0
w(ixmin1:ixmax1,ixmin2:ixmax2,v3_)=0.d0

eps = 0.05d0
phase = 3.0d0
sigma=0.02d0
mid_pt = (xprobmax2-xprobmin2)/2.0d0

dy = -abs(ya(3)-ya(2)) !distance between 2 pts is const thro out the array

!NOTE if you change the number of mxnest remember to change dr also.
!I should make a common var for this!!!!
do ix2=ixmin2,ixmax2
do ix1=ixmin1,ixmax1
   na=floor(((x(ix1,ix2,2)-(xprobmin2))/dr)+1.0d0+3.0d0) !(1)+3 for ghost cells
   w(ix1,ix2,rho_) = rhoa(na)
!! Useful checks to make sure correct values are taken from arrays
!   if (mype==0 .and. first_1) then
!      print *, 'check start value', na, w(ix^D,rho_)*runit, x(ix1,ix2,2)
!   endif
!   first_1=.false.
!   if (mype==7 .and. ix2 == ixmax2 .and. ix1 == ixmax1) then
!   print *, 'check last value', na , w(ix^D,rho_)*runit, x(ix1,ix2,2)
!   endif
end do
end do

patchw(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=.false.
call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,x,&
   patchw)
w(ixmin1:ixmax1,ixmax2+1,e_)=iniene !-2004.2883700596335
!print *, w(ix^S,e_)
call primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,&
   x)
!print *, w(ix^S,p_)

do ix1=ixmin1,ixmax1
do ix2=ixmax2,ixmin2,-1
! current issue is that I have not defined values in ghost cells so w(ixmin1:ixmax1,ixmax2+1,p_)=0
!Need to fix this as giving bart simpson at the end of my plot
   w(ix1,ix2,p_)=w(ix1,ix2+1,p_)+w(ix1,ix2,rho_)*dy*eqpar(grav2_)
enddo
enddo


!do ix_2=ixGlo2,ixGhi2
!do ix_1=ixGlo1+2,ixGhi1-2
!   w(ix_1,ix_2,rho_)=-(1.D0/eqpar(grav1_))*(1.D0/(12.D0*(x(ix_1+1,ix_2,1)-x(ix_1,ix_2,1))))*(w(ix_1+2,ix_2,p_)-8.D0*w(ix_1+1,ix_2,p_)+8.D0*w(ix_1-1,ix_2,p_)-w(ix_1-2,ix_2,p_))
!enddo
!enddo


!!This is where we add the jet
!For FWHM
sigma = 0.2d0
jet_w = 0.2d0
jet_h = 0.05d0
jet_cx = (jet_w-jet_w)/2.0d0 !center pts x
jet_cy = (jet_h-0.0d0)/2.0d0 !center pts y
r_jet(ixmin1:ixmax1,ixmin2:ixmax2) = (x(ixmin1:ixmax1,ixmin2:ixmax2,1)&
   -jet_cx)**2+(x(ixmin1:ixmax1,ixmin2:ixmax2,2)-jet_cy)**2

!do ix2=ixmin2,ixmax2
!do ix1=ixmin1,ixmax1
!   na=floor(((x(ix1,ix2,2)-(xprobmin2))/dr)+1)
!   if (x(ix^D,2).le. jet_h .and. x(ix^D,1).le.jet_w/2.0d0 .and. x(ix^D,1).ge.-jet_w/2.0d0) then
!      w(ix^D,rho_) = rhoa(na)
!      w(ix^D,p_) = pa(na)!10.0d0*pa(1)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
!      w(ix^D,v2_)  = (J_sp)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
!      w(ix^D,tr1_) = 100.0d0
!   endif
!end do
!end do



w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)  =eqpar(BB1_)
w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)  =eqpar(BB2_)
w(ixmin1:ixmax1,ixmin2:ixmax2,b3_)  =eqpar(BB3_)

!ixGmin1:ixGmax1,ixGmin2:ixGmax2 !whole domian, including ghost cells
!ixmin1,ixmax1,ixmin2,ixmax2 ! excludes ghosts cells


patchw(ixGmin1:ixGmax1,ixGmin2:ixGmax2)=.false.
call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,x,&
   patchw)

end subroutine initonegrid_usr
!=============================================================================
!INCLUDE:amrvacnul/specialbound.t
!=============================================================================
subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies
include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iw,iB
double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixO1,ixO2, ix1,ix2, ixIntmin1,&
   ixIntmin2,ixIntmax1,ixIntmax2

double precision :: dx1,dx2
double precision :: pth(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
logical :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!----------------------------------------------------------------------------
select case(iB)
case(3)
!-------------------------------
! This creates cont condition for special BC
!-------------------------------

   !yrange
   ixImin1=ixGmin1 != 1
   ixImax1=ixGmin1-1+dixB !=3
   !xrange
   ixImin2=ixGmin2 !=1
   ixImax2=ixGmax2 !=16

   ixIntmin1=ixOmin1;ixIntmin2=ixOmin2;ixIntmax1=ixOmax1;ixIntmax2=ixOmax2;
   ixIntmin2=ixOmin2;ixIntmax2=ixOmax2+1;
!   ixIntmin2=ixOmin2-1+dixB;ixIntmax2=ixOmax2+1;

!   print *, ixIntmin2, ixIntmax2

   call getpthermal(w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixIntmin1,ixIntmin2,&
      ixIntmax1,ixIntmax2,pth)
   do ix1=ixImin1,ixImax1
      w(ixImin2:ixImax2,ix1,rho_) = w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,p_) = pth(ixImin2:ixImax2,ixImax1+1)
      w(ixImin2:ixImax2,ix1,v1_) = w(ixImin2:ixImax2,ixImax1&
         +1,m1_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,v2_) = w(ixImin2:ixImax2,ixImax1&
         +1,m2_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,v3_) = w(ixImin2:ixImax2,ixImax1&
         +1,m3_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,b1_) = w(ixImin2:ixImax2,ixImax1+1,b1_)
      w(ixImin2:ixImax2,ix1,b2_) = w(ixImin2:ixImax2,ixImax1+1,b2_)
      w(ixImin2:ixImax2,ix1,b3_) = w(ixImin2:ixImax2,ixImax1+1,b3_)
   end do

   patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
   call conserve(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,w,x,patchw)
end select



end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x)

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

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!==========================================================================
!INCLUDE:amrvacnul/simple/specialsource.t
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
!-----------------------------------------------------------------------------

call addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

end subroutine specialsource
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

call getdt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,&
   dtnew,dx1,dx2,x)

end subroutine getdt_special
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
subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)=wB0(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0

!=============================================================================
! amrvacusr.t.nul
!=============================================================================
