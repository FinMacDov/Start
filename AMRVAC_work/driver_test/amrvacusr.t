!=============================================================================
!=============================================================================
! INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/simple/speciallog_mhd.t
INCLUDE:amrvacmodules/simple/gravity.t
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/simple/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
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
punit=  2.3d0*nHunit*k_B*Teunit ! m-3.J.K-1.K = m−1.kg.s−2 = pa
Bunit= dsqrt(miu0*punit) ! sqrt(H.m-2.kg.s-2) = sqrt(kg2.s-4.A-2)= kg.A.s-2 = T  
vunit=  Bunit/dsqrt(miu0*runit) ! m/s
UNIT_VELOCITY=vunit
tunit=Lunit/vunit ! s
heatunit=punit/tunit ! m−1.kg.s−3 = pa.s-1
Ti = tunit! s

! units for convert
if(iprob==-1) then
  normvar(0) = one 
else
  normvar(0) = UNIT_LENGTH
endif
normvar(rho_) = UNIT_DENSITY
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normvar(pp_)     = UNIT_VELOCITY**2 * UNIT_DENSITY
{^C&normvar(b^C_)   = dsqrt(4.0d0*dpi*normvar(pp_)) \}
normt = UNIT_LENGTH/UNIT_VELOCITY

eqpar(grav1_)=0.d0
eqpar(grav2_)= -275.4229*Lunit/vunit**2 !where [Lunit/vunit**2] = [s^2/m] 

dr=(xprobmax2-xprobmin2)/dble(jmax) ! step size
!dr=(xprobmax2-xprobmin2)/(dble(jmax)-4)

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
 read(1,*) p(i) !dyn/cm^2
 read(2,*) mu(i) !dimension mean molecular weight
 read(3,*) rho(i) !g/cm^3
 read(4,*) Tem(i) !k
 read(5,*) Z(i) !Mm
end do 

!Need to covert CGS to SI and then make dimensionless

do j=1,jmax
   Tea(j)=Tem(j)/Teunit 
   rhoa(j)=rho(j)*1000.d0/runit 
   pa(j)=0.1d0*p(j)/punit 
   mua(j) = mu(j) ! This is already dimensionless. This mu is the average mol wieght 
   ya(j) = Z(j)*1000000/Lunit
enddo

!equation ideal gas law (igl): p = R_gas*rho*T
do k=1,jmax
   pigl(k) = ((R_gas/mu(k))*rho(k)*Tem(k))/punit  
enddo

iniene = pigl(k-1)/(1-eqpar(gamma_))-(eqpar(BB1_)*eqpar(BB1_)+eqpar(BB2_)*eqpar(BB2_)+eqpar(BB3_)*eqpar(BB3_))/2
print *, iniene

!ya(j-1)*Lunit this gives last element of matrix

J_d = rhoa(1) !jet density

! this creates txt file pruns with the output shown below
if (mype==0) then
 open(123,file='output_test',form='formatted')
 write(123,*) jmax, '  ya(ix)              ','pa(ix)                   ','pidgl              ','Temp(ix)                  ','rho(ix)                       ', 'mu(ix)           '
 do ix=1,jmax
    write(123,*) ya(ix), pa(ix)*punit, pigl(ix)*punit, Tea(ix)*Teunit, rhoa(ix)*runit, mua(ix) 
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
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
integer :: ix^D,na,imode
double precision:: res, lxsize, sigma, phase, mid_pt, eps, dy
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision:: rinlet(ixG^T), r_jet(ixG^T), p_bg(ixG^T)
double precision:: psi(ixG^T)
DOUBLE PRECISION :: jet_w, jet_h, jet_cx, jet_cy

logical, save:: first=.true.
logical patchw(ixG^T)
!-----------------------------------------------------------------------------

{^IFONED call mpistop("This is a multi-D MHD problem") }

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

w(ix^S,v1_)=0.d0
w(ix^S,v2_)=0.d0
w(ix^S,v3_)=0.d0

eps = 0.05d0
phase = 3.0d0
sigma=0.02d0
mid_pt = (xprobmax2-xprobmin2)/2 

dy = -abs(ya(3)-ya(2))

do ix2=ixmin2,ixmax2
do ix1=ixmin1,ixmax1
   na=floor(((x(ix1,ix2,2)-(xprobmin2))/dr)+1)
   w(ix^D,rho_) = rhoa(na)
   w(ix^D,rho_) = pa(na)
end do
end do

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)
w(ixmin1:ixmax1,ixmax2,e_)=iniene !-2004.2883700596335
!print *, w(ix^S,e_)
call primitive(ixG^L,ix^L,w,x)
print *, w(ix^S,p_)

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

!For FWHM 
sigma = 0.2d0
jet_w = 0.2d0
jet_h = 0.05d0
jet_cx = (jet_w-jet_w)/2.0d0 !center pts x
jet_cy = (jet_h-0.0d0)/2.0d0 !center pts y
r_jet(ix^S) = (x(ix^S,1)-jet_cx)**2+(x(ix^S,2)-jet_cy)**2  

do ix2=ixmin2,ixmax2
do ix1=ixmin1,ixmax1
   na=floor(((x(ix1,ix2,2)-(xprobmin2))/dr)+1)
   if (x(ix^D,2).le. jet_h .and. x(ix^D,1).le.jet_w/2.0d0 .and. x(ix^D,1).ge.-jet_w/2.0d0) then
      w(ix^D,rho_) = rhoa(na)
      w(ix^D,p_) = pa(na)!10.0d0*pa(1)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
      w(ix^D,v2_)  = (J_sp)*dexp(-r_jet(ix1,ix2)/(sigma*sigma))
      w(ix^D,tr1_) = 100.0d0
   endif
end do
end do

w(ix^S,b1_)  =eqpar(BB1_)
w(ix^S,b2_)  =eqpar(BB2_)
w(ix^S,b3_)  =eqpar(BB3_)

!ixGmin1:ixGmax1,ixGmin2:ixGmax2 !whole domian, including ghost cells
!ixmin1,ixmax1,ixmin2,ixmax2 ! excludes ghosts cells


patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
!INCLUDE:amrvacnul/specialbound.t
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies
include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ixO^L,iw,iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
integer :: ixI^L, ixO^D, ix1,ix2, ixInt^L

double precision :: dx^D
double precision :: pth(ixG^T)
logical :: patchw(ixG^T)
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
 
   ixInt^L=ixO^L;
   ixIntmin2=ixOmin2;ixIntmax2=ixOmax2+1;
!   ixIntmin2=ixOmin2-1+dixB;ixIntmax2=ixOmax2+1;

!   print *, ixIntmin2, ixIntmax2

   call getpthermal(w,x,ixG^L,ixInt^L,pth)
   do ix1=ixImin1,ixImax1
      w(ixImin2:ixImax2,ix1,rho_) = w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,p_) = pth(ixImin2:ixImax2,ixImax1+1)
      w(ixImin2:ixImax2,ix1,v1_) = w(ixImin2:ixImax2,ixImax1+1,m1_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,v2_) = w(ixImin2:ixImax2,ixImax1+1,m2_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,v3_) = w(ixImin2:ixImax2,ixImax1+1,m3_)/w(ixImin2:ixImax2,ixImax1+1,rho_)
      w(ixImin2:ixImax2,ix1,b1_) = w(ixImin2:ixImax2,ixImax1+1,b1_)
      w(ixImin2:ixImax2,ix1,b2_) = w(ixImin2:ixImax2,ixImax1+1,b2_)
      w(ixImin2:ixImax2,ix1,b3_) = w(ixImin2:ixImax2,ixImax1+1,b3_)
   end do

   patchw(ixO^S)=.false.
   call conserve(ixG^L,ixO^L,w,x,patchw)
end select



end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

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

integer, intent(in) :: ixI^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!==========================================================================
!INCLUDE:amrvacnul/simple/specialsource.t
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

call addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------


end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest

!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0

!=============================================================================
! amrvacusr.t.nul
!=============================================================================
