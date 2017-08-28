!=============================================================================
subroutine errest
 
use mod_forest, only: refine, buffer
include 'amrvacdef.f'

integer :: igrid, iigrid, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision :: factor
!-----------------------------------------------------------------------------
if (igridstail==0) return

select case (errorestimate)
case (0) 
   ! all refinement solely based on user routine specialrefine_grid
case (1) 
   ! Richardson procedure: compare coarse-integrate versus integrate-coarse
   ! done with low order, dimensionally unsplit scheme typelow1
   ! For error estimate: compare 1st order, coarse 2*dx,2*dt step  with
   ! 1st order dt step, since otherwise wCT can not be filled by interpolation

   ! Note:when there are sources: only unsplit sources are taken into account

   ! Note: the low order solutions are obtained with dimsplit=F.
   !       When overall scheme uses dimsplit=T
   !       and courantpar>0.5, the low order can be unstable. In
   !       that case, simplify this scheme, set low order step of size dt
   !       and compare coarse with available solution at t_n. Enforce this
   !       through `skipfinestep' 


   if (skipfinestep) then
      factor=one
   else
      factor=two
   end if

   ! call the integrator on a grid twice as coarse.

   ixCoGmin1=1;ixCoGmin2=1;
   ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;

   call createCoarse(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2)
   call getbc(t,0.d0,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,pwCoarse,pwCoCo,&
      pgeoCoarse,pgeoCoCo,.true.,0,nwflux)
   call advectCoarse(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,factor)

   ! Now advance full grid and obtain relative error
   ! between coarse and full grid sln at time t_n+1
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call errest1_grid(igrid,pw(igrid)%w)
   end do
!$OMP END PARALLEL DO

case (2) 
   ! simply compare w_n-1 with w_n and trigger refinement on relative
   ! differences
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call compare1_grid(igrid,pwold(igrid)%w,pw(igrid)%w)
   end do
!$OMP END PARALLEL DO

case (3)
   ! Error estimation is based on Lohner's scheme
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_grid(igrid)
   end do
!$OMP END PARALLEL DO

case (4)
   ! Error estimation is based on Lohner's original scheme
!$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call lohner_orig_grid(igrid)
   end do
!$OMP END PARALLEL DO


case default
   call mpistop("Unknown error estimator")
end select

! enforce additional refinement on e.g. coordinate and/or time info here
if (nbufferx1/=0.or.nbufferx2/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,refine,&
   ngridshi*npe,MPI_LOGICAL,MPI_LOR, icomm,ierrmpi)
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call forcedrefine_grid(igrid,pw(igrid)%w)
end do
!$OMP END PARALLEL DO

if (nbufferx1/=0.or.nbufferx2/=0) buffer=.false.

end subroutine errest
!=============================================================================
subroutine lohner_grid(igrid)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, idims2, level
integer :: ixmin1,ixmin2,ixmax1,ixmax2, hxmin1,hxmin2,hxmax1,hxmax2, jxmin1,&
   jxmin2,jxmax1,jxmax2, h2xmin1,h2xmin2,h2xmax1,h2xmax2, j2xmin1,j2xmin2,&
   j2xmax1,j2xmax2, ix1,ix2
double precision :: epsilon, tolerance
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2) :: numerator,&
    denominator, error
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: tmp, tmp1, tmp2
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,&
      ixGlo2,ixGhi1,ixGhi2,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(pw(igrid)%w(jxmin1:jxmax1,&
             jxmin2:jxmax2,iflag))-dlog10(pw(igrid)%w(hxmin1:hxmax1,&
             hxmin2:hxmax2,iflag))
        else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2)=pw(igrid)%w(jxmin1:jxmax1,&
             jxmin2:jxmax2,iflag)-pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
             iflag)
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
        else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(jxmin1:jxmax1,jxmin2:jxmax2)&
             -tmp1(hxmin1:hxmax1,hxmin2:hxmax2)
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
         h2xmax1=ixMhi1-kr(1,idims2);h2xmax2=ixMhi2-kr(2,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
         j2xmax1=ixMhi1+kr(1,idims2);j2xmax2=ixMhi2+kr(2,idims2);
         numerator=numerator+(tmp(j2xmin1:j2xmax1,j2xmin2:j2xmax2)&
            -tmp(h2xmin1:h2xmax1,h2xmin2:h2xmax2))**2.0d0
      end do
   end do
   denominator=zero
   do idims=1,ndim
      if (iflag<=nw) then
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iflag)))
         else
          tmp=dabs(pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iflag))
         end if
      else
         if (logflag(iiflag)) then
          tmp=dabs(dlog10(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2)))
         else
          tmp=dabs(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2))
         end if
      end if
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
      tmp2(ixmin1:ixmax1,ixmin2:ixmax2)=tmp(jxmin1:jxmax1,jxmin2:jxmax2)&
         +tmp(hxmin1:hxmax1,hxmin2:hxmax2)
      hxmin1=ixMlo1-2*kr(1,idims);hxmin2=ixMlo2-2*kr(2,idims)
      hxmax1=ixMhi1-2*kr(1,idims);hxmax2=ixMhi2-2*kr(2,idims);
      jxmin1=ixMlo1+2*kr(1,idims);jxmin2=ixMlo2+2*kr(2,idims)
      jxmax1=ixMhi1+2*kr(1,idims);jxmax2=ixMhi2+2*kr(2,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,iflag))-dlog10(pw(igrid)%w&
             (ixMlo1:ixMhi1,ixMlo2:ixMhi2,iflag))) +dabs(dlog10(pw(igrid)%w&
             (ixMlo1:ixMhi1,ixMlo2:ixMhi2,iflag))-dlog10(pw(igrid)%w&
             (hxmin1:hxmax1,hxmin2:hxmax2,iflag)))
        else
           tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(pw(igrid)%w(jxmin1:jxmax1,&
              jxmin2:jxmax2,iflag)-pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
              iflag)) +dabs(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,iflag)&
              -pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,iflag))
        end if
      else
        if (logflag(iiflag)) then
          tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2))-dlog10(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2))) &
             +dabs(dlog10(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2))&
             -dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2)))
        else
           tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(tmp1(jxmin1:jxmax1,&
              jxmin2:jxmax2)-tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2)) &
              +dabs(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2)-tmp1(hxmin1:hxmax1,&
              hxmin2:hxmax2))
        end if
      end if
      do idims2=1,ndim
         h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
         h2xmax1=ixMhi1-kr(1,idims2);h2xmax2=ixMhi2-kr(2,idims2);
         j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
         j2xmax1=ixMhi1+kr(1,idims2);j2xmax2=ixMhi2+kr(2,idims2);
         denominator=denominator +(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)&
            +amr_wavefilter(level)*(tmp2(j2xmin1:j2xmax1,j2xmin2:j2xmax2)&
            +tmp2(h2xmin1:h2xmax1,h2xmin2:h2xmax2)))**2
      end do
   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.
tolerance=tol(level)
do ix2=ixMlo2,ixMhi2
do ix1=ixMlo1,ixMhi1

   if (error(ix1,ix2) >= tolerance) then
      refineflag(ix1,ix2) = .true.
   else if (error(ix1,ix2) <= tolratio(level)*tolerance) then
      coarsenflag(ix1,ix2) = .true.
   end if
end do
end do


if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level<mxnest) &
   refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level&
   >1) coarsen(igrid,mype)=.true.

end subroutine lohner_grid
!=============================================================================
subroutine lohner_orig_grid(igrid)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: iiflag, iflag, idims, level
integer :: ixmin1,ixmin2,ixmax1,ixmax2, hxmin1,hxmin2,hxmax1,hxmax2, jxmin1,&
   jxmin2,jxmax1,jxmax2, ix1,ix2
double precision :: epsilon
double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2) :: numerator,&
    denominator, error
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: dp, dm, dref, tmp1
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
epsilon=1.0d-6
level=node(plevel_,igrid)
ixmin1=ixMlo1;ixmin2=ixMlo2;ixmax1=ixMhi1;ixmax2=ixMhi2;

error=zero
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   numerator=zero
   denominator=zero
   if (iflag>nw)call specialvarforerrest(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,&
      ixGlo2,ixGhi1,ixGhi2,iflag,pw(igrid)%w,tmp1)
   do idims=1,ndim
      hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
      hxmax1=ixmax1-kr(1,idims);hxmax2=ixmax2-kr(2,idims);
      jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
      jxmax1=ixmax1+kr(1,idims);jxmax2=ixmax2+kr(2,idims);
      if (iflag<=nw) then
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(pw(igrid)%w(jxmin1:jxmax1,&
             jxmin2:jxmax2,iflag))-dlog10(pw(igrid)%w(ixmin1:ixmax1,&
             ixmin2:ixmax2,iflag))
          dm(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(pw(igrid)%w(ixmin1:ixmax1,&
             ixmin2:ixmax2,iflag))-dlog10(pw(igrid)%w(hxmin1:hxmax1,&
             hxmin2:hxmax2,iflag))
          dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(dlog10(pw(igrid)%w&
             (jxmin1:jxmax1,jxmin2:jxmax2,iflag)))+ 2.0d0 * &
             dabs(dlog10(pw(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,iflag))) &
             + dabs(dlog10(pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,iflag)))
        else
          dp(ixmin1:ixmax1,ixmin2:ixmax2)=pw(igrid)%w(jxmin1:jxmax1,&
             jxmin2:jxmax2,iflag)-pw(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
             iflag)
          dp(ixmin1:ixmax1,ixmin2:ixmax2)=pw(igrid)%w(ixmin1:ixmax1,&
             ixmin2:ixmax2,iflag)-pw(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
             iflag)
          dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=dabs(pw(igrid)%w(jxmin1:jxmax1,&
             jxmin2:jxmax2,iflag))+2.0d0*dabs(pw(igrid)%w(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,iflag)) +dabs(pw(igrid)%w(hxmin1:hxmax1,&
             hxmin2:hxmax2,iflag))
        end if
      else
        if (logflag(iiflag)) then
          dp(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2))-dlog10(tmp1(ixmin1:ixmax1,ixmin2:ixmax2))
          dm(ixmin1:ixmax1,ixmin2:ixmax2)=dlog10(tmp1(ixmin1:ixmax1,&
             ixmin2:ixmax2))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
          dref(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(dlog10(tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2)))+ 2.0d0 * dabs(dlog10(tmp1(ixmin1:ixmax1,&
             ixmin2:ixmax2))) + dabs(dlog10(tmp1(hxmin1:hxmax1,&
             hxmin2:hxmax2)))
        else
          dp(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(jxmin1:jxmax1,jxmin2:jxmax2)&
             -tmp1(ixmin1:ixmax1,ixmin2:ixmax2)
          dm(ixmin1:ixmax1,ixmin2:ixmax2)=tmp1(ixmin1:ixmax1,ixmin2:ixmax2)&
             -tmp1(hxmin1:hxmax1,hxmin2:hxmax2)
          dref(ixmin1:ixmax1,ixmin2:ixmax2)=dabs(tmp1(jxmin1:jxmax1,&
             jxmin2:jxmax2))+2.0d0*dabs(tmp1(ixmin1:ixmax1,ixmin2:ixmax2)) &
             +dabs(tmp1(hxmin1:hxmax1,hxmin2:hxmax2))
        end if
      end if

      numerator(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=numerator+(dp(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2)-dm(ixMlo1:ixMhi1,ixMlo2:ixMhi2))**2.0d0

      denominator(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=denominator &
         + (dabs(dp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)) + dabs(dm(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2)) + amr_wavefilter(level)*dref(ixMlo1:ixMhi1,&
         ixMlo2:ixMhi2))**2.0d0

   end do
   error=error+wflags(iiflag)*dsqrt(numerator/max(denominator,epsilon))
end do

refineflag=.false.
coarsenflag=.false.

do ix2=ixMlo2,ixMhi2
do ix1=ixMlo1,ixMhi1
   if (error(ix1,ix2) >= tol(level)) then
      refineflag(ix1,ix2) = .true.
   else if (error(ix1,ix2) <= tolratio(level)*tol(level)) then
      coarsenflag(ix1,ix2) = .true.
   end if
end do
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level<mxnest) &
   refine(igrid,mype)=.true.
if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level&
   >1) coarsen(igrid,mype)=.true.

end subroutine lohner_orig_grid
!=============================================================================
subroutine compare1_grid(igrid,wold,w)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
    w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)

integer :: ix1,ix2, iiflag, iflag, level
double precision :: epsilon
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
! identify the points to be flagged in two steps:
!  step I: compare w_n-1 with w_n solution, store flags in auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

epsilon=1.0d-6

refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.
coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.
level=node(plevel_,igrid)

do ix2=ixMlo2,ixMhi2 
do ix1=ixMlo1,ixMhi1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      averages(iflag) = w(ix1,ix2,iflag)-wold(ix1,ix2,iflag)
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (abs(wold(ix1,ix2,iflag))<smalldouble)then
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,ix2,&
            iflag))+epsilon)
      else
         error=error+wflags(iiflag)* abs(averages(iflag))/(abs(wold(ix1,ix2,&
            iflag)))
      end if
   end do
   if (error >= tol(level)) then
      refineflag(ix1,ix2) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(ix1,ix2) = .true.
   end if
end do
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level&
      >1) coarsen(igrid,mype)=.true.
end if

end subroutine compare1_grid
!=============================================================================
subroutine createCoarse(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2)
include 'amrvacdef.f'

integer, intent(in) :: ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2

integer :: iigrid, igrid
integer :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3
!-----------------------------------------------------------------------------
ixGmin1=ixCoGmin1;ixGmin2=ixCoGmin2;ixGmax1=ixCoGmax1;ixGmax2=ixCoGmax2;
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call createCoarse_grid(igrid,pwCoarse(igrid),pxCoarse(igrid),ixGmin1,&
      ixGmin2,ixGmax1,ixGmax2,pwold(igrid)%w,px(igrid)%x)
end do
!$OMP END PARALLEL DO

end subroutine createCoarse
!=============================================================================
subroutine createCoarse_grid(igrid,pwCo,pxCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2,wold,xold)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision :: wold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
    xold(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
type(walloc) pwCo
type(xalloc) pxCo

integer :: ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2
!-----------------------------------------------------------------------------
dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

! now coarsen by 2 in every direction - conservatively
! coarse grid dimension half its size
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmax1=ixCoGmax1-dixB
ixCoMmax2=ixCoGmax2-dixB;
call coarsen_grid(wold,xold,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
   ixMhi2,pwCo%w,pxCo%x,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixCoMmin1,&
   ixCoMmin2,ixCoMmax1,ixCoMmax2, pgeo(igrid),pgeoCoarse(igrid),&
   coarsenprimitive,.true.)

end subroutine createCoarse_grid
!=============================================================================
subroutine advectCoarse(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,factor)
include 'amrvacdef.f'

integer :: ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision, intent(in) :: factor

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call advectCoarse_grid(igrid,pwCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmax1,&
      ixCoGmax2,factor)
end do
!$OMP END PARALLEL DO

end subroutine advectCoarse
!=============================================================================
subroutine advectCoarse_grid(igrid,pwCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2,factor)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision, intent(in) :: factor
type(walloc) pwCo

double precision :: qdt, dx1,dx2
double precision, dimension(:,:,:), allocatable :: wCo1
double precision :: fC(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nwflux,&
   1:ndim)
integer :: level
!-----------------------------------------------------------------------------
dxlevel(1)=two*rnode(rpdx1_,igrid);dxlevel(2)=two*rnode(rpdx2_,igrid);
dx1=dxlevel(1);dx2=dxlevel(2);

! here we integrate on the coarse grid
allocate(wCo1(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw))
wCo1(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nwflux)=pwCo%w&
   (ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nwflux)

! 1st order scheme: do coarse time step of
!    size 2*dt starting from t_n-1 solution in pwCoarse
!    to arrive at t_n+1 (n-index from normal uncoarsened grid)
! result in pwCoarse: coarse solution at t_n+1

if (.not.slab) mygeo => pgeoCoarse(igrid)

qdt=factor*dt_grid(igrid)
level=node(plevel_,igrid)
call advect1_grid(typelow1(level),qdt,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
   1,ndim,t,wCo1,t, pwCo%w,wCo1,fC,dx1,dx2,pxCoarse(igrid)%x)

deallocate(wCo1)

end subroutine advectCoarse_grid
!=============================================================================
subroutine errest1_grid(igrid,w)

include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer :: level, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2
double precision :: dx1,dx2, qdt
double precision :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux,1:ndim),&
    wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

dx1=dx(1,level);dx2=dx(2,level);
dxlevel(1)=dx1;dxlevel(2)=dx2;

wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:nwflux)

if (.not.skipfinestep) then
   if (.not.slab) mygeo => pgeo(igrid)

   qdt=dt_grid(igrid)
   call advect1_grid(typelow1(level),qdt,ixGlo1,ixGlo2,ixGhi1,ixGhi2,1,ndim,t&
      +qdt,w, t+qdt,wFi,w,fC,dx1,dx2,px(igrid)%x)
end if

ixCoGmin1=1;ixCoGmin2=1;
ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;

call flagbadpoints(wFi,pwCoarse(igrid)%w,ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2,igrid,level)

end subroutine errest1_grid
!=============================================================================
subroutine flagbadpoints(w,wCo,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,igrid,&
   level)

! compare error between coarse and fine solution in wCo, w 
! We base the comparison on the physical field selected by the index flag_ 
!
! on entry:
!  w:       normal time integration 
!  wCo:     time integration on coarsened grid (2*dx)

use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in)         :: igrid, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
    level
double precision,intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
    wCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,nw)
double precision            :: specialvar(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
    specialvarCo(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2)

logical :: needgetaux
integer :: iCo1,iCo2, iFi1,iFi2, ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2,&
    iiflag, iflag
double precision :: average, error
double precision :: averages(nflag_)
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag, coarsenflag
!-----------------------------------------------------------------------------
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmax1=ixCoGmax1-dixB
ixCoMmax2=ixCoGmax2-dixB;

needgetaux=.false.
do iiflag=1,flags(nflag_); iflag=flags(iiflag);
  if (iflag>nwflux) needgetaux=.true.
end do
if (nwaux>0.and.needgetaux) then
   saveigrid=igrid
   if(.not.slab)mygeo=>pgeo(igrid)
   call getaux(.true.,w,px(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,&
      ixMhi1,ixMhi2,'flagbadpoints')
   call getaux(.true.,wCo,pxCoarse(igrid)%x,ixCoGmin1,ixCoGmin2,ixCoGmax1,&
      ixCoGmax2,ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2,'flagbadpointsCo')
end if

! identify the points to be flagged in two steps (needed!):
!  step I: compare coarse with fine solution, store flags in fine auxiliary
!  step II: transfer flags from auxiliary to refine and coarsen

refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.
coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = .false.


do iiflag=1,flags(nflag_); iflag=flags(iiflag);
   if (iflag>nw) then
      call specialvarforerrest(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
         ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,iflag,wCo,specialvarCo)
      call specialvarforerrest(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,&
         ixGhi1,ixGhi2,iflag,w,specialvar)
   end if
end do

iFi2 = ixMlo2
do iCo2 = ixCoMmin2,ixCoMmax2 
iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   average=zero
   error=zero
   do iiflag=1,flags(nflag_); iflag=flags(iiflag);
      if (slab) then
         if (iflag<=nw) averages(iflag)=sum(w(iFi1:iFi1+1,iFi2:iFi2&
            +1,iflag))/two**ndim
         if (iflag>nw)  averages(iflag)=sum(specialvar(iFi1:iFi1+1,iFi2:iFi2&
            +1))/two**ndim
      else
         if (iflag<=nw) averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1&
            +1,iFi2:iFi2+1) *w(iFi1:iFi1+1,iFi2:iFi2+1,iflag))&
            /pgeoCoarse(igrid)%dvolume(iCo1,iCo2)
         if (iflag>nw)  averages(iflag)=sum(pgeo(igrid)%dvolume(iFi1:iFi1&
            +1,iFi2:iFi2+1) *specialvar(iFi1:iFi1+1,iFi2:iFi2&
            +1))/pgeoCoarse(igrid)%dvolume(iCo1,iCo2)
      end if
      average=average+wflags(iiflag)*abs(averages(iflag))
      if (iflag<=nw) error=error+wflags(iiflag)*abs(averages(iflag)&
         -wCo(iCo1,iCo2,iflag))
      if (iflag> nw) error=error+wflags(iiflag)*abs(averages(iflag)&
         -specialvarCo(iCo1,iCo2))
   end do
   if (abs(average)>smalldouble) then
      error=error/average
   else
      write(unitterm,*)'Warning from flagbadpoints: zero average:',average
      write(unitterm,*)'   wCo(iCo^D,1:nw):',wCo(iCo1,iCo2,1:nw),' indices:',&
         iCo1,iCo2
      write(unitterm,*)'On grid:',igrid,' at level ',level
      write(unitterm,*)'   and grid indices : ',node(pig1_,igrid),node(pig2_,&
         igrid)
      write(unitterm,*)'cell indices : ',iCo1,iCo2
      call mpistop("")
   end if
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1,iFi2:iFi2+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1,iFi2:iFi2+1) = .true.
   end if
   iFi1 = iFi1+2
end do
   iFi2 = iFi2+2
end do

iFi2 = ixMlo2
do iCo2 = ixCoMmin2,ixCoMmax2 
iFi1 = ixMlo1
do iCo1 = ixCoMmin1,ixCoMmax1 
   if (error >= tol(level)) then
      refineflag(iFi1:iFi1+1,iFi2:iFi2+1) = .true.
   else if (error <= tolratio(level)*tol(level)) then
      coarsenflag(iFi1:iFi1+1,iFi2:iFi2+1) = .true.
   end if
   iFi1 = iFi1+2
end do
   iFi2 = iFi2+2
end do

if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2))) then
   if (level<mxnest) refine(igrid,mype)=.true.
end if
if (time_advance) then
   if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)).and.level&
      >1) coarsen(igrid,mype)=.true.
end if

end subroutine flagbadpoints
!=============================================================================
subroutine forcedrefine_grid(igrid,w)
use mod_forest, only: coarsen, refine, buffer
include 'amrvacdef.f'

integer, intent(in) :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer :: level
integer :: my_refine, my_coarsen
double precision :: qt
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

! initialize to 0
my_refine   = 0
my_coarsen  = 0

if (time_advance) then
   qt=t+dt
else
   if (errorestimate==1.or.errorestimate==2) then
      qt=t+dt
   else
      qt=t
   end if
end if
   
call specialrefine_grid(igrid,level,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,&
   ixMhi1,ixMhi2,qt,w,px(igrid)%x, my_refine,my_coarsen)

if (my_coarsen==1) then
   if (level>1) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
endif

if (my_coarsen==-1)then
   coarsen(igrid,mype)=.false.
end if

if (my_refine==1) then
   if (level<mxnest) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
   else
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
   end if
end if

if (my_refine==-1) then
  refine(igrid,mype)=.false.
end if

if (nbufferx1/=0.or.nbufferx2/=0) then
   if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
      refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=.true.
      call refinebuffer(igrid,refineflag)
   end if
end if

end subroutine forcedrefine_grid
!=============================================================================
subroutine forcedrefine_grid_io(igrid,w)
use mod_forest, only: coarsen, refine
include 'amrvacdef.f'

integer, intent(in)          :: igrid
double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer                   :: level, my_levmin, my_levmax
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: refineflag
!-----------------------------------------------------------------------------
level=node(plevel_,igrid)

if (level_io > 0) then
   my_levmin = level_io
   my_levmax = level_io
else
   my_levmin = max(1,level_io_min)
   my_levmax = min(mxnest,level_io_max)
end if


if (level>my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
elseif (level<my_levmin) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
end if

if (level==my_levmin .or. level==my_levmax) then
  refine(igrid,mype)=.false.
  coarsen(igrid,mype)=.false.
end if


if(refine(igrid,mype).and.level>=mxnest)refine(igrid,mype)=.false.
if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

end subroutine forcedrefine_grid_io
!=============================================================================
subroutine refinebuffer(igrid,refineflag)
use mod_forest, only: refine, buffer
include 'amrvacdef.f'

integer, intent(in) :: igrid
logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2), intent(in) :: refineflag

integer :: ishiftbuf1,ishiftbuf2, i1,i2, ixmin1,ixmin2,ixmax1,ixmax2,&
    ineighbor, ipe_neighbor, level
!-----------------------------------------------------------------------------
ishiftbuf1=ixMhi1-ixMlo1-nbufferx1+1;ishiftbuf2=ixMhi2-ixMlo2-nbufferx2+1;
do i2=-1,1
do i1=-1,1
   ixmin1=max(ixMlo1,ixMlo1+i1*ishiftbuf1)
   ixmin2=max(ixMlo2,ixMlo2+i2*ishiftbuf2);
   ixmax1=min(ixMhi1,ixMhi1+i1*ishiftbuf1)
   ixmax2=min(ixMhi2,ixMhi2+i2*ishiftbuf2);
   if (ixmax1<ixmin1.or.ixmax2<ixmin2) cycle
   if (any(refineflag(ixmin1:ixmax1,ixmin2:ixmax2))) then
      select case (neighbor_type(i1,i2,igrid))
      case (2)
         ineighbor=neighbor(1,i1,i2,igrid)
         ipe_neighbor=neighbor(2,i1,i2,igrid)
         if (.not.refine(ineighbor,ipe_neighbor)) then
            buffer(ineighbor,ipe_neighbor)=.true.
            refine(ineighbor,ipe_neighbor)=.true.
         end if
      case (3)
         level=node(plevel_,igrid)
         if (level<mxnest) then
            ineighbor=neighbor(1,i1,i2,igrid)
            ipe_neighbor=neighbor(2,i1,i2,igrid)
            if (.not.refine(ineighbor,ipe_neighbor)) then
               buffer(ineighbor,ipe_neighbor)=.true.
               refine(ineighbor,ipe_neighbor)=.true.
            end if
         end if
      end select
   end if
end do
end do

end subroutine refinebuffer
!=============================================================================
