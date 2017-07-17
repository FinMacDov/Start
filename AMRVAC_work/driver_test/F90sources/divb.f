
!=============================================================================
subroutine addsource_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add divB related sources to w within ixO
! corresponding to Powel
include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, wCT(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in) :: dx1,dx2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
integer :: iw
double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,divb)

! e = e - qdt (v . b) * div b
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)-&
     qdt*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)*wCT(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,b1_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
        *wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+wCT(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,m3_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_))&
        /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)


! b = b - qdt v * div b
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

! m = m - qdt b * div b
if (B0field) then

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m1_)-qdt*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)+myB0_cell%w&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)) &
        *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m2_)-qdt*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+myB0_cell%w&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)) &
        *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m3_)-qdt*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)+myB0_cell%w&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)) &
        *divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
else

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m1_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)*divb(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m2_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)*divb(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      m3_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)*divb(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)
end if

! since this option changes energy: smallvalues call
if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,"addsource_powel")

end subroutine addsource_powel
!=============================================================================
subroutine addsource_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add divB related sources to w within ixO
! corresponding to Janhunen, just the term in the induction equation.
include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, wCT(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, intent(in) :: dx1,dx2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
!.. local ..
double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,divb)

! b = b - qdt v * div b
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)&
   -qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)/wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end subroutine addsource_janhunen
!=============================================================================
subroutine addsource_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,dx1,dx2)

! Add Linde's divB related sources to wnew within ixO
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
double precision, intent(in)    :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(in)    :: dx1,dx2
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
integer :: iw, idims, ixmin1,ixmin2,ixmax1,ixmax2, ixpmin1,ixpmin2,ixpmax1,&
   ixpmax2, i1,i2, iside
double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2),graddivb&
   (ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

! Calculate div B
ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
call getdivb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,ixmax2,&
   divb)
! for AMR stability, retreat one cell layer from the boarders of level jump
ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmax1=ixOmax1;ixpmax2=ixOmax2;
do idims=1,ndim
  select case(idims)
   case(1)
      do iside=1,2
        i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
        if(leveljump(i1,i2)) then
          if(iside==1) then
            ixpmin1=ixOmin1-i1
          else
            ixpmax1=ixOmax1-i1
          end if 
        end if
      end do
   
   case(2)
      do iside=1,2
        i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
        if(leveljump(i1,i2)) then
          if(iside==1) then
            ixpmin2=ixOmin2-i2
          else
            ixpmax2=ixOmax2-i2
          end if 
        end if
      end do
   
  end select
end do

! Add Linde's diffusive terms
do idims=1,ndim
   ! Calculate grad_idim(divb)
   select case(typegrad)
   case("central")
     call gradient(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
        ixpmax1,ixpmax2,idims,graddivb)
   case("limited")
     call gradientS(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
        ixpmax1,ixpmax2,idims,graddivb)
   end select
   !ixmin^D=ixpmin^D;
   !ixmax^D=merge(ixmin^D,ixmax^D,kr(idims,^D)==1);
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=merge(ixmax^D,ixmin^D,kr(idims,^D)==1);
   !ixmax^D=ixpmax^D;
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=ixpmin^D+kr(idims,^D);
   !ixmax^D=ixpmax^D-kr(idims,^D);
   !call gradientS(divb,ixI^L,ix^L,idims,graddivb)

   ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
   if (slab) then
      graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
         ixpmin2:ixpmax2)*divbdiff/(1.0d0/dxlevel(1)**2+1.0d0/dxlevel(2)**2)
   else
      graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
         ixpmin2:ixpmax2)*divbdiff /(1.0d0/mygeo%dx(ixpmin1:ixpmax1,&
         ixpmin2:ixpmax2,1)**2+1.0d0/mygeo%dx(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
         2)**2)
   end if
   do iw= iwmin,iwmax
      if (iw==b0_+idims) then
         ! B_idim += eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,iw)=w(ixpmin1:ixpmax1,&
            ixpmin2:ixpmax2,iw)+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)

      else if (iw==e_ .and. typedivbdiff=='all') then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,iw)=w(ixpmin1:ixpmax1,&
            ixpmin2:ixpmax2,iw)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,b0_&
            +idims)*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)

      end if
   end do
end do

if(fixsmall) call smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,&
   ixpmin2,ixpmax1,ixpmax2,"addsource_linde")

end subroutine addsource_linde
!=============================================================================
subroutine getdivb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,divb)

! Calculate div B within ixO

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
    ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2)

double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndir)


!-----------------------------------------------------------------------------

bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
   ixImin2:ixImax2,b0_+1:b0_+ndir)


select case(typediv)
case("central")
  call divvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,divb)
case("limited")
  call divvectorS(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divb)
end select

end subroutine getdivb
!=============================================================================
