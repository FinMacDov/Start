!#############################################################################
! module amrvacphys/- mhd

!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------
minrho= max(zero,smallrho)


minp  = max(zero,smallp)
smalle= minp/(eqpar(gamma_)-one)

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! set default values for entropy fixes

include 'amrvacdef.f'

integer :: il
!-----------------------------------------------------------------------------

eqpar(gamma_)=5.d0/3.d0





eqpar(eta_)=0.d0
eqpar(etahyper_)=0.d0
eqpar(etah_)=0.d0

do il=1,nw
   select case(il)
   case(fastRW_,fastLW_,slowRW_,slowLW_)
      entropycoef(il)=0.2d0
   case(alfvRW_,alfvLW_)
      entropycoef(il)=0.4d0
   case default
      entropycoef(il)= -one
   end select
end do

end subroutine initglobaldata
!=============================================================================
subroutine getaux(clipping,w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,subname)

! Calculate auxilary variables ixO^L from non-auxiliary entries in w
! clipping can be set to .true. to e.g. correct unphysical pressures,
! densities, v>c,  etc.

include 'amrvacdef.f'

logical                :: clipping
integer                :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision       :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
character(len=*)       :: subname
!-----------------------------------------------------------------------------

end subroutine getaux
!=============================================================================
subroutine checkw(checkprimitive,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,w,flag)

include 'amrvacdef.f'

logical :: checkprimitive
integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2)

double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------
flag(ixImin1:ixImax1,ixImin2:ixImax2)=.true.


if(checkprimitive)then
   tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)
else
   ! First calculate kinetic energy*2=m**2/rho
   tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,m3_)**2 )/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
   ! Add magnetic energy*2=b**2
   tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2&
      +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b3_)**2
   ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
   tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(eqpar(gamma_)-one)&
      *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)-half*tmp1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))
endif

flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   >=minp .and. w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)>=minrho)


end subroutine checkw
!=============================================================================
subroutine conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

integer, intent(in)    :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision       :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
logical                :: patchw(ixImin1:ixImax1,ixImin2:ixImax2)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

   ! Calculate total energy from pressure, kinetic and magnetic energy
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      p_)*invgam+&
      half*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)*(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,v1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)**2&
         +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)**2)+(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
         +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))

   ! Convert velocity to momentum
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_);

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_)

end where

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"conserve")

end subroutine conserve
!=============================================================================
subroutine conserven(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,patchw)

! Transform primitive variables into conservative ones

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
logical, intent(in)             :: patchw(ixImin1:ixImax1,ixImin2:ixImax2)
double precision       :: invgam
!-----------------------------------------------------------------------------

invgam=1.d0/(eqpar(gamma_)-one)
where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

   ! Calculate total energy from pressure, kinetic and magnetic energy
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      p_)*invgam+&
      half*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)*(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,v1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)**2&
         +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)**2)+(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
         +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))

   ! Convert velocity to momentum
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_);

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_)

end where

end subroutine conserven
!=============================================================================
subroutine primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x)

! Transform conservative variables into primitive ones

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
! .. local ..

integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)       :: patchierror

integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"primitive")

! Convert momentum to velocity
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_);

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=(eqpar(gamma_)-one)*(w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,e_)- &
       half*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,v2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)&
          **2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
       + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+ w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,b2_)**2+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))
if(strictsmall) then
  if(any(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)<minp)) then
    lowpindex=minloc(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_))
    lowpindex(1)=lowpindex(1)+ixOmin1-1;lowpindex(2)=lowpindex(2)+ixOmin2-1;
    write(*,*)'too small pressure = ',minval(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,p_)),' with limit=',&
    minp,' at x=',x(lowpindex(1),lowpindex(2),1:ndim),' array index&
       =',lowpindex,&
    ' where E_k=',half*(w(lowpindex(1),lowpindex(2),v1_)**2&
       +w(lowpindex(1),lowpindex(2),v2_)**2+w(lowpindex(1),lowpindex(2),v3_)&
       **2)*&
    w(lowpindex(1),lowpindex(2),rho_),&
    ' E_B=',half*(w(lowpindex(1),lowpindex(2),b1_)**2+w(lowpindex(1),&
       lowpindex(2),b2_)**2+w(lowpindex(1),lowpindex(2),b3_)**2),' E_total=',&
    w(lowpindex(1),lowpindex(2),p_)/(eqpar(gamma_)-one)+half*&
    (w(lowpindex(1),lowpindex(2),v1_)**2+w(lowpindex(1),lowpindex(2),v2_)**2&
       +w(lowpindex(1),lowpindex(2),v3_)**2)*w(lowpindex(1),lowpindex(2),&
       rho_)+&
    half*(w(lowpindex(1),lowpindex(2),b1_)**2+w(lowpindex(1),lowpindex(2),&
       b2_)**2+w(lowpindex(1),lowpindex(2),b3_)**2),' w(1:nwflux)=',&
    w(lowpindex(1),lowpindex(2),1:nwflux),' when t=',t,' it=',it
    call mpistop("=== primitive pressure problem===")
  end if
else
  if (strictgetaux) then
     where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)<minp)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=minp
     endwhere
  else
     where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)<minp)
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1
     elsewhere
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
     end where
     if (any(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)) &
   call correctaux(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,w,x,patchierror,'primitive')
 end if
end if


! We got rho, Dtr, now we can get the tracers:
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,Dtr1_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

end subroutine primitive
!=============================================================================
subroutine primitiven(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,patchw)

! Transform conservative variables into primitive ones

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
logical, intent(in),dimension(ixImin1:ixImax1,ixImin2:ixImax2)   :: patchw
!-----------------------------------------------------------------------------

where(.not.patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  ! Convert momentum to velocity
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     m1_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     m2_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     m3_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_);

  ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=(eqpar(gamma_)-one)&
     *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)- &
         half*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v1_)**2+w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,v2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v3_)&
            **2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
         + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+ w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b2_)**2+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
            **2))


! We got rho, Dtr, now we can get the tracers:
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tr1_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,Dtr1_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

end where

end subroutine primitiven
!=============================================================================
subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x)

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision,intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
!-----------------------------------------------------------------------------

w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rhos_)=(eqpar(gamma_)-one)&
   *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**(one-eqpar(gamma_)) &
               *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)-half&
                  *((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2&
                  +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)**2&
                  +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
                  /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                                   +(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                                      **2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                                      b2_)**2+w(ixOmin1:ixOmax1,&
                                      ixOmin2:ixOmax2,b3_)**2)))


end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
!-----------------------------------------------------------------------------

w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=(one/(eqpar(gamma_)-one))&
   *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**(eqpar(gamma_)-one)*&
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rhos_)+half*((w&
                (ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                m3_)**2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+&
             (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                b3_)**2))


end subroutine rhos_to_e
!=============================================================================
subroutine internalenergy(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,ie)

! get internal energy

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision                :: ie(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

ie(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
   **2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,m3_)**2 )/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
   + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2
! internal energy=e-0.5*(2ek+2eb)
ie(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)&
   -half*ie(ixOmin1:ixOmax1,ixOmin2:ixOmax2)



end subroutine internalenergy
!=============================================================================
subroutine getv(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,idims,v)

! Calculate v_idim=m_idim/rho within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, idims
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

v(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
   +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,idims,cmax,cmin,needcmin)

! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

include 'amrvacdef.f'

logical :: new_cmax,needcmin
integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, idims
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw), &
   cmax(ixImin1:ixImax1,ixImin2:ixImax2), cmin(ixImin1:ixImax1,&
   ixImin2:ixImax2)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

double precision :: csound2(ixImin1:ixImax1,ixImin2:ixImax2),&
    cfast2(ixImin1:ixImax1,ixImin2:ixImax2), AvMinCs2(ixImin1:ixImax1,&
   ixImin2:ixImax2), tmp(ixImin1:ixImax1,ixImin2:ixImax2)

!, vh(ixI^S), va2(ixI^S), Omegai(ixI^S)
!-----------------------------------------------------------------------------

!Direction independent part of getcmax:
call getcsound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,csound2)

if (B0field) then
   cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=( (w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b1_)+myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))**2&
      +(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+myB0%w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,2))**2+(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
      +myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3))**2 ) /w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)+csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=cfast2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-4.0d0*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *((myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)+w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b0_+idims))**2) /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
else
   cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=( w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
      +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2 )/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)+csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=cfast2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-4.0d0*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims))**2)/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)
endif

where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
   AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
end where

AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(AvMinCs2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2))


tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dsqrt(half*(cfast2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))



if(needcmin)then
  cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(+tmp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2) +(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
     +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)),zero)
  cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(-tmp(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2) +(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
     +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)),zero)
else
  cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
     +dabs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+idims)/w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,rho_))
end if

end subroutine getcmax
!=============================================================================
subroutine getpthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,p)

! Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw), p(ixImin1:ixImax1,&
   ixImin2:ixImax2)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)       :: patchierror
integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,'getpthermal')

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(eqpar(gamma_)-one)*(w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,e_)- &
       half*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
       +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))

! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)) then
    lowpindex=minloc(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    lowpindex(1)=lowpindex(1)+ixOmin1-1;lowpindex(2)=lowpindex(2)+ixOmin2-1;
    write(*,*)'too small pressure = ',minval(p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),' with limit=',minp,&
    ' at x=',x(lowpindex(1),lowpindex(2),1:ndim),' array index=',lowpindex,&
    ' where E_k=',half*(w(lowpindex(1),lowpindex(2),m1_)**2&
       +w(lowpindex(1),lowpindex(2),m2_)**2+w(lowpindex(1),lowpindex(2),m3_)&
       **2)/&
    w(lowpindex(1),lowpindex(2),rho_),&
    ' E_B=',half*(w(lowpindex(1),lowpindex(2),b1_)**2+w(lowpindex(1),&
       lowpindex(2),b2_)**2+w(lowpindex(1),lowpindex(2),b3_)&
       **2),' E_total=',w(lowpindex(1),lowpindex(2),e_),&
    ' w(1:nwflux)=',w(lowpindex(1),lowpindex(2),1:nwflux),&
    ' when t=',t,' it=',it
    call mpistop("=== strictsmall in getpthermal ===")
  end if
else
  if (strictgetaux) then
     where(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)
        p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=minp
     endwhere
  else
     where(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1
     elsewhere
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
     end where
     if (any(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0))then
       call correctaux(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x,patchierror,'getpthermal')
       where(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
         p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(eqpar(gamma_)-one)&
            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)- &
        half*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
           /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
           +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))
       end where
     end if
 end if
end if


end subroutine getpthermal
!=============================================================================
subroutine getcsound2prim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,csound2)

! Calculate the square of the thermal sound speed csound2 within ixO^L
! from the primitive variables in w.
! csound2=gamma*p/rho

include 'amrvacdef.f'

integer, intent(in)             :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixImin1,&
   ixImin2,ixImax1,ixImax2
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=eqpar(gamma_)*w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)


end subroutine getcsound2prim
!=============================================================================
subroutine getcsound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,csound2)

! Calculate the square of the thermal sound speed csound2 within ixO^L.
! csound2=gamma*p/rho

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

call getpthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,csound2)
csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=eqpar(gamma_)*csound2&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)


end subroutine getcsound2
!=============================================================================
subroutine getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,p)

! Calculate total pressure within ixO^L including magnetic pressure
! p=(g-1)*e-0.5*(g-1)*m**2/rho+(1-0.5*g)*b**2

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2)
!.. local ..

double precision :: gamma

integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)       :: patchierror
integer, dimension(ndim)       :: lowpindex
!-----------------------------------------------------------------------------
if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,'getptotal')


gamma=eqpar(gamma_)
p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(one-half*gamma)*( w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
   +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2 )+(gamma-one)*&
  (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)-half*(w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)**2&
     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)/w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,rho_))

if(strictsmall) then
  if(any(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)) then
    lowpindex=minloc(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    lowpindex(1)=lowpindex(1)+ixOmin1-1;lowpindex(2)=lowpindex(2)+ixOmin2-1;
    write(*,*)'too small pressure = ',minval(p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),' at x=',&
    x(lowpindex(1),lowpindex(2),1:ndim),lowpindex,' with limit&
       =',minp,' where E_k=',&
    half*(w(lowpindex(1),lowpindex(2),m1_)**2+w(lowpindex(1),lowpindex(2),&
       m2_)**2+w(lowpindex(1),lowpindex(2),m3_)**2)/w(lowpindex(1),&
       lowpindex(2),rho_),' E_B=',&
    half*(w(lowpindex(1),lowpindex(2),b1_)**2+w(lowpindex(1),lowpindex(2),&
       b2_)**2+w(lowpindex(1),lowpindex(2),b3_)**2),'E_total&
       =',w(lowpindex(1),lowpindex(2),e_),&
    ' w(1:nwflux)=',w(lowpindex(1),lowpindex(2),1:nwflux),' when t=',t,' it&
       =',it
    call mpistop("=== strictsmall in getptotal ===")
  end if
else
  if (strictgetaux) then
     where(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)
        p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=minp
     endwhere
  else
     where(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<minp)
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1
     elsewhere
       patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0
     end where
     if (any(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0))then
       call correctaux(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,w,x,patchierror,'getptotal')
       where(patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0)
          p(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(one-half*gamma)&
             *( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
             **2 )+(gamma-one)*&
            (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)-half*(w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
               **2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
               /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
       end where
     end if
 end if
end if




end subroutine getptotal
!=============================================================================
subroutine getfluxforhllc(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iw, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)
!.. local ..
logical :: transport
integer          :: idir
double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
   ixImin2:ixImax2)
!-----------------------------------------------------------------------------
transport=.true.

if (B0field) then
   if (iw==m0_+idims .or. iw==e_) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      =myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b1_)+myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)&
      *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+myB0%w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
end if

select case (iw)
   ! f_i[rho]=v_i*rho
   case (rho_)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero

  case (tr1_)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero

   ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
   case (m1_)
      if (idims==1) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp2)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=tmp2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
            =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)+tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_) &
                                -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                                   +idims)*myB0%w(ixOmin1:ixOmax1,&
                                   ixOmin2:ixOmax2,1)
      end if
   case (m2_)
      if (idims==2) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp2)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=tmp2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
            =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)+tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) &
                                -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                                   +idims)*myB0%w(ixOmin1:ixOmax1,&
                                   ixOmin2:ixOmax2,2)
      end if
   case (m3_)
      if (idims==3) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp2)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=tmp2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
            =f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)+tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) &
                                -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                                   +idims)*myB0%w(ixOmin1:ixOmax1,&
                                   ixOmin2:ixOmax2,3)
      end if

   ! f_i[e]=v_i*e+(m_i*ptotal-b_i*(b_k*m_k))/rho
   case (e_)
      call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,tmp2)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,m0_+idims)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)- &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)*( w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,b1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,m2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_) ))/w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,rho_)

      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+ tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
                          *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                             +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                          -( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
                             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) )&
                             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                          *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)
      end if

   ! f_i[b_k]=v_i*b_k-m_k/rho*b_i
   case (b1_)
      if (idims==1) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero 

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if
   case (b2_)
      if (idims==2) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero 

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if
   case (b3_)
      if (idims==3) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero 

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if

   case default
      call mpistop("Error in getflux: unknown flow variable!")
end select

end subroutine getfluxforhllc
!=============================================================================
subroutine getflux(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iw, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision,intent(out)    :: f(ixImin1:ixImax1,ixImin2:ixImax2)
!.. local ..
logical :: transport
double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)
integer          :: idirmin, idir
!-----------------------------------------------------------------------------
transport=.true.

if (B0field) then
   if (iw==m0_+idims .or. iw==e_) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      =myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b1_)+myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)&
      *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+myB0%w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
end if

select case (iw)
   ! f_i[rho]=v_i*rho
   case (rho_)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero

  case (tr1_)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero

   ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
   case (m1_)
      if (idims==1) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,f)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_) &
                          -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                             +idims)*myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      end if
   case (m2_)
      if (idims==2) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,f)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b2_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) &
                          -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                             +idims)*myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
      end if
   case (m3_)
      if (idims==3) then
         call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,f)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
            +idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
         if (B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b3_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)
      end if
      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
            *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) &
                          -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
                             +idims)*myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)
      end if

   ! f_i[e]=v_i*e+(m_i*ptotal-b_i*(b_k*m_k))/rho
   case (e_)
      call getptotal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,f)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         m0_+idims)*f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)- &
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)*( w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,b1_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,m2_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_) ))/w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,rho_)

      if (B0field) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+ tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
                          *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                             +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                          -( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
                             +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
                             *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) )&
                             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                          *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)
      end if

   ! f_i[b_k]=v_i*b_k-m_k/rho*b_i
   case (b1_)
      if (idims==1) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if
   case (b2_)
      if (idims==2) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if
   case (b3_)
      if (idims==3) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero

         transport=.false.
      else
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= -w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,b0_+idims)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
            /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=f(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) &
                     +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                        +idims)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                        *myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3) &
                     -myB0%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)&
                        *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
                        /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if

      end if

   case default
      call mpistop("Error in getflux: unknown flow variable!")
end select

end subroutine getflux
!=============================================================================
subroutine addsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
logical, intent(in)             :: qsourcesplit
!.. local ..
double precision:: dx1,dx2
!-----------------------------------------------------------------------------

dx1=dxlevel(1);dx2=dxlevel(2);
if(qsourcesplit .eqv. ssplitresis) then
! Sources for resistivity in eqs. for e, B1, B2 and B3
if(dabs(eqpar(eta_))>smalldouble)then
   if (.not.slab) call mpistop("no resistivity in non-slab geometry")
   if(compactres)then
      call addsource_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
   else
      call addsource_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
   endif
endif

if (eqpar(etahyper_)>0.d0)then
   call addsource_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
end if
endif


if(qsourcesplit .eqv. ssplitdivb) then
! Sources related to div B
select case (typedivbfix)

case ('powel')
   call addsource_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
case ('janhunen')
   call addsource_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
case ('linde')
   call addsource_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
case ('lindejanhunen')
   call addsource_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
   call addsource_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
case ('lindepowel')
   call addsource_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
   call addsource_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)
end select
endif

end subroutine addsource
!=============================================================================
subroutine addsource_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add resistive source to w within ixO 
! Uses 3 point stencil (1 neighbour) in each direction, non-conservative
! If the fourthorder precompiler flag is set, uses fourth order central difference for the laplacian. Then the stencil is 5 (2 neighbours).  
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in)    :: dx1,dx2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!.. local ..
integer :: ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,idirmin,iw,idims,&
   jxOmin1,jxOmin2,jxOmax1,jxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,ix


double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
   ixImin2:ixImax2)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7&
   -2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2)
double precision :: gradeta(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
!-----------------------------------------------------------------------------


! Calculating resistive sources involve one extra layer
ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

if(ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop&
   ("Error in addsource_res1: Non-conforming input limits")

! Calculate current density and idirmin
call getcurrent(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,idirmin,current)

if(eqpar(eta_)>zero)then
   eta(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(eta_)
   gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)=zero
else
   call specialeta(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idirmin,x,current,eta)
   ! assumes that eta is not function of current?
   do idims=1,ndim
      call gradient(eta,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,idims,tmp)
      gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
   enddo
endif

do idir=1,ndir

   ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp

   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
   tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
      b0_+idir)
   do idims=1,ndim
      jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
      jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+&
           (tmp2(jxOmin1:jxOmax1,jxOmin2:jxOmax2)-2.0d0*tmp2(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)+tmp2(hxOmin1:hxOmax1,hxOmin2:hxOmax2))&
              /dxlevel(idims)**2
   enddo



   ! Multiply by eta
   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
   if(eqpar(eta_)<zero)then
      do jdir=1,ndim; do kdir=idirmin,3
         if(lvc(idir,jdir,kdir)/=0)then
            if(lvc(idir,jdir,kdir)==1)then
               tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2)-gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
            else
               tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2)+gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
            endif
         endif
      enddo; enddo
   endif

   ! Add sources related to eta*laplB-grad(eta) x J to B and e
   do iw=1,nw
      if(iw==b0_+idir)then
         ! dB_idir/dt+=tmp
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

      else if(iw==e_)then
         ! de/dt+=B.tmp
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
            *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir)

      endif
   end do  ! iw
enddo ! idir



! de/dt+=eta*J**2
tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
do idir=idirmin,3
   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      +current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)**2
enddo
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)&
   +qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"addsource_res1")

end subroutine addsource_res1
!=============================================================================
subroutine addsource_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add resistive source to w within ixO 
! Uses 5 point stencil (2 neighbours) in each direction, conservative

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in)    :: dx1,dx2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!.. local ..
integer :: ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,idirmin,iw,idims,&
   idirmin1

double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
   ixImin2:ixImax2)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7&
   -2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2),curlj(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:3)
double precision :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,1:3),&
   tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;

if(ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop&
   ("Error in addsource_res2: Non-conforming input limits")

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
! Determine exact value of idirmin while doing the loop.
call getcurrent(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idirmin,current)

if(eqpar(eta_)>zero)then
   eta(ixmin1:ixmax1,ixmin2:ixmax2)=eqpar(eta_)
else
   call specialeta(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idirmin,x,current,eta)
endif

! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
tmpvec(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
do jdir=idirmin,3
   tmpvec(ixmin1:ixmax1,ixmin2:ixmax2,jdir)=current(ixmin1:ixmax1,&
      ixmin2:ixmax2,jdir)*eta(ixmin1:ixmax1,ixmin2:ixmax2)*qdt
enddo
call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,curlj,idirmin1,1,3)
 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
    b1_)-curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_-b0_)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     b2_)-curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_-b0_)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     b3_)-curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_-b0_)


! de/dt= +div(B x Jeta)
tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixmin1:ixmax1,ixmin2:ixmax2)=wCT(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
         +jdir)*current(ixmin1:ixmax1,ixmin2:ixmax2,kdir)*eta(ixmin1:ixmax1,&
         ixmin2:ixmax2)*qdt
      if(lvc(idir,jdir,kdir)==1)then
         tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,idir)=tmpvec2(ixmin1:ixmax1,&
            ixmin2:ixmax2,idir)+tmp(ixmin1:ixmax1,ixmin2:ixmax2)
      else
         tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,idir)=tmpvec2(ixmin1:ixmax1,&
            ixmin2:ixmax2,idir)-tmp(ixmin1:ixmax1,ixmin2:ixmax2)
      endif
   endif
enddo; enddo; enddo
!select case(typediv)
!case("central")
   call divvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,tmp)
!case("limited")
 !call divvectorS(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,tmp)
!end select
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)&
   +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"addsource_res2")

end subroutine addsource_res2
!=============================================================================
subroutine addsource_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add Hyper-resistive source to w within ixO 
! Uses 9 point stencil (4 neighbours) in each direction.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in)    :: dx1,dx2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!.. local ..
double precision                :: current(ixImin1:ixImax1,ixImin2:ixImax2,7&
   -2*ndir:3)
double precision                :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:3),tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,1:3),tmp(ixImin1:ixImax1,&
   ixImin2:ixImax2),ehyper(ixImin1:ixImax1,ixImin2:ixImax2,1:3)
integer                         :: ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,&
   idirmin,idirmin1
!-----------------------------------------------------------------------------
ixmin1=ixOmin1-3;ixmin2=ixOmin2-3;ixmax1=ixOmax1+3;ixmax2=ixOmax2+3;
if(ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop&
   ("Error in addsource_hyperres: Non-conforming input limits")

call getcurrent(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idirmin,current)
tmpvec(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
do jdir=idirmin,3
   tmpvec(ixmin1:ixmax1,ixmin2:ixmax2,jdir)=current(ixmin1:ixmax1,&
      ixmin2:ixmax2,jdir)
enddo

ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,tmpvec2,idirmin1,1,3)

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
tmpvec(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
call curlvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,tmpvec,idirmin1,1,3)
ehyper(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir) = - tmpvec(ixmin1:ixmax1,&
   ixmin2:ixmax2,1:ndir)*eqpar(etahyper_)

ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
call curlvector(ehyper,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,tmpvec2,idirmin1,1,3)

 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
    b1_)-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*qdt
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     b2_)-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)*qdt
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     b3_)-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)*qdt



! de/dt= +div(B x Ehyper)
ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   tmpvec2(ixmin1:ixmax1,ixmin2:ixmax2,idir) = tmpvec(ixmin1:ixmax1,&
      ixmin2:ixmax2,idir)&
        + lvc(idir,jdir,kdir)*wCT(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
           +jdir)*ehyper(ixmin1:ixmax1,ixmin2:ixmax2,kdir)
enddo; enddo; enddo
tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
call divvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,tmp)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)&
   +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*qdt

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"addsource_hyperres")

end subroutine addsource_hyperres
!=============================================================================
subroutine getcurrent(w,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idirmin,current)

! Calculate idirmin and the idirmin:3 components of the common current array
! make sure that dxlevel(^D) is set correctly.
include 'amrvacdef.f'

integer, parameter:: idirmin0=7-2*ndir
integer :: ixmin1,ixmin2,ixmax1,ixmax2, idirmin, ixImin1,ixImin2,ixImax1,&
   ixImax2
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7&
   -2*ndir:3),bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
!-----------------------------------------------------------------------------

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
!bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
call curlvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,current,idirmin,idirmin0,ndir)

end subroutine getcurrent
!=============================================================================
subroutine getdt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,dtnew,dx1,dx2,x)

! If resistivity is not zero, check diffusion time limit for dt

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(out)   :: dtnew
double precision, intent(in)    :: dx1,dx2
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
!.. local ..
integer :: idirmin,idims
double precision :: dxarr(ndim)
double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7&
   -2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2) 
!double precision :: {#IFDEF HALL , dthall }
!-----------------------------------------------------------------------------
dtnew=bigdouble

dxarr(1)=dx1;dxarr(2)=dx2;
dxlevel(1)=dx1;dxlevel(2)=dx2;
if(eqpar(eta_)>zero)then
   dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/eqpar(eta_)
else if(eqpar(eta_)<zero)then
   call getcurrent(w,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idirmin,current)
   call specialeta(w,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idirmin,x,current,eta)
   dtnew=bigdouble
   do idims=1,ndim
      dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixmin1:ixmax1,&
         ixmin2:ixmax2)/dxarr(idims)**2)))
   enddo
endif

if(eqpar(etahyper_)>zero)then
   dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/eqpar(etahyper_),dtnew)
end if 

!{#IFDEF HALL
! This is now covered in cmax, so no need.
!if(eqpar(etah_)>zero)then
!   call getdthall(w,x,ixI^L,ix^L,dx^D,dthall)
!   dtnew=min(dtnew,dthall)
!end if
!}
end subroutine getdt
!=============================================================================
subroutine ppmflatcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,w,&
   d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,&
   ixRmax1,ixRmax2
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
   d2w(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)
double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2),&
   dp(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

if(useprimitive)then
 drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =eqpar(gamma_)*dabs(d2w&
    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))&
              /min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,rho_),w(ixRmin1:ixRmax1,&
                 ixRmin2:ixRmax2,rho_))
 dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dabs(d2w(ixOmin1:ixOmax1,&
    ixOmin2:ixOmax2,p_))/min(w(ixLmin1:ixLmax1,ixLmin2:ixLmax2,p_),&
    w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,p_))
end if


end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,ixLmax1,&
   ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,ixRRmin2,ixRRmax1,&
   ixRRmax2,idims,w,drho,dp,dv)

include 'amrvacdef.f'

! based on Mignone and Miller and Collela 2002
! PPM flattening at shocks: we use total pressure and not thermal pressure 

integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,&
   ixLmin2,ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,ixRRmin1,ixRRmin2,&
   ixRRmax1,ixRRmax2
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)

double precision, intent(inout) :: drho(ixImin1:ixImax1,ixImin2:ixImax2),&
   dp(ixImin1:ixImax1,ixImin2:ixImax2),dv(ixImin1:ixImax1,ixImin2:ixImax2)
double precision :: ptot(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------


if(useprimitive)then
   ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
   ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      p_)+half*( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2 )
   where (dabs(ptot(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2)-ptot&
      (ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2))>smalldouble)
      drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dabs((ptot(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2)-ptot(ixLmin1:ixLmax1,ixLmin2:ixLmax2))&
                        /(ptot(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2)&
                           -ptot(ixLLmin1:ixLLmax1,ixLLmin2:ixLLmax2)))
   elsewhere
      drho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
   end where

   !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
   !  use "dp" to save squared sound speed, assume primitive in w
   dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(eqpar(gamma_)*w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

   dp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = dabs(ptot(ixRmin1:ixRmax1,&
      ixRmin2:ixRmax2)-ptot(ixLmin1:ixLmax1,ixLmin2:ixLmax2))&
                /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)*dp(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2))
   ! recycle ptot to store v
   ptot(ixImin1:ixImax1,ixImin2:ixImax2)= w(ixImin1:ixImax1,ixImin2:ixImax2,&
      v0_+idims)
   call gradient(ptot,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,idims,dv)
end if



end subroutine ppmflatsh
!=============================================================================
subroutine addgeometry(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!.. local ..
integer          :: iw,idir, h1xmin1,h1xmin2,h1xmax1,h1xmax2, h2xmin1,h2xmin2,&
   h2xmax1,h2xmax2
double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)
logical          :: angmomfix=.false.
!-----------------------------------------------------------------------------

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('cylindrical')
   do iw=1,nwflux
      select case (iw)
      ! s[mr]=(ptotal-Bphi**2+mphi**2/rho)/radius
      case (mr_)
         call getptotal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
            /x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero


      end select

      ! Divide by radius and add to w
      if (iw==mr_) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
             /x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      endif

   end do
case ('spherical')
   h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2);h1xmax1=ixOmax1-kr(1,1)
   h1xmax2=ixOmax2-kr(1,2);  h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
   h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);
   do iw=1,nwflux
      select case (iw)
      ! s[m1]=((mtheta**2+mphi**2)/rho+2*ptotal-(Btheta**2+Bphi**2))/r
      case (m1_)
         call getptotal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp)
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               1)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)+myB0_cell%w&
               (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)*wCT(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,b2_)+myB0_cell%w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
         end if
         ! For nonuniform Cartesian grid this provides hydrostatic equil.
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) &
            *(mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
            -mygeo%surfaceC1(h1xmin1:h1xmax1,h1xmin2:h1xmax2)) &
            /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
               +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)**2&
                  /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                  -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2 &
               +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2&
                  /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                  -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2
         if (B0field.and.ndir>1) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)-2.0d0*myB0_cell%w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,2) *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
               -2.0d0*myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3) &
               *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
         end if

      ! s[m2]=-(mr*mtheta/rho-Br*Btheta)/r
      !       + cot(theta)*(mphi**2/rho+(p+0.5*B**2)-Bphi**2)/r
      case (m2_)


         call getptotal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,tmp)
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               1)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)+myB0_cell%w&
               (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)*wCT(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,b2_)+myB0_cell%w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)
         end if
         ! This will make hydrostatic p=const an exact solution
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
                     *(mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                        -mygeo%surfaceC2(h2xmin1:h2xmax1,h2xmin2:h2xmax2)) &
                     /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)


         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,m1_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
            /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                     -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                        *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               1)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) &
                                 +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                                    *myB0_cell%w(ixOmin1:ixOmax1,&
                                    ixOmin2:ixOmax2,2)
         end if



         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2&
            /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                        -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
                           **2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)) &
                                           /dsin(x(ixOmin1:ixOmax1,&
                                              ixOmin2:ixOmax2,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)-2.0d0*myB0_cell%w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
                       *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))&
                          /dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
         end if

      ! s[m3]=-(mphi*mr/rho-Bphi*Br)/r
      !       -cot(theta)*(mtheta*mphi/rho-Btheta*Bphi)/r
      case (m3_)
         if (.not.angmomfix) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,m3_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
               /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                     -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
                        *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_))  &
                   -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                      *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
                      /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) &
                     -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
                        *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)) &
                   *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))&
                      /dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
            if (B0field) then
               tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2)+myB0_cell%w(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,1)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  b3_) &
                          +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
                             *myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                             3)  &
                          +(myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)&
                             *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) &
                            +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
                               *myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                               3)) &
                   *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))&
                      /dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
            end if
         end if



      ! s[b2]=(mr*Btheta-mtheta*Br)/rho/r
      !       + cot(theta)*psi/r
      case (b2_)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,m1_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) &
                    -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                       *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_))&
                       /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
               *myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) &
                       -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
                          *myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))&
                          /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if



      ! s[b3]=(mr*Bphi-mphi*Br)/rho/r
      !       -cot(theta)*(mphi*Btheta-mtheta*Bphi)/rho/r
      case (b3_)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,m1_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_) &
                 -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)*wCT&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_))/wCT&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)  &
                -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)*wCT&
                   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_) &
                 -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)*wCT&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_))*dcos(x&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)) &
                               /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                                  *dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
               *myB0_cell%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3) &
               -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)*myB0_cell%w&
                  (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))/wCT(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,rho_) &
               -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)*myB0_cell%w&
                  (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) &
                -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)*myB0_cell%w&
                   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,3))*dcos(x&
                   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)) &
                               /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
                                  *dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
         end if

      end select
      ! Divide by radius and add to w
      if (iw==m1_.or.iw==m2_.or.iw==b2_ .or.iw==b3_ .or.(iw==m3_&
         .and..not.angmomfix)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
         =w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)+qdt*tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
   end do
end select

end subroutine addgeometry
!=============================================================================
! end module amrvacphys/- mhd
!#############################################################################
