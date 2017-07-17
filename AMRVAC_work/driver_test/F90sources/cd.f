!#############################################################################
! module cd
! Centered difference scheme
!=============================================================================
subroutine centdiff(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)

! Advance the iws flow variables from t to t+qdt within ixO^L by centered 
! differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1,dx2
double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

double precision :: v(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ndim), f(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
double precision :: dxinv(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2
logical :: transport
!-----------------------------------------------------------------------------

! An extra layer is needed in each direction for which fluxes are added.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
   ixmax1=ixmax1+kr(idims,1);ixmax2=ixmax2+kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) then
   call mpistop("Error in CentDiff: Non-conforming input limits")
end if

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;

! Add fluxes to w
do idims= idimmin,idimmax
   ixmin1=ixOmin1-kr(idims,1);ixmin2=ixOmin2-kr(idims,2)
   ixmax1=ixOmax1+kr(idims,1);ixmax2=ixOmax2+kr(idims,2); ixCmin1=ixmin1
   ixCmin2=ixmin2; ixCmax1=ixOmax1;ixCmax2=ixOmax2;
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2); 
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

   call getv(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idims,f)
   v(ixmin1:ixmax1,ixmin2:ixmax2,idims)=f(ixmin1:ixmax1,ixmin2:ixmax2)

   do iw=1,nwflux
      ! Get non-transported flux
      call getflux(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
         ixmax2,iw,idims,f,transport)
      ! Add transport flux
      if (transport) f(ixmin1:ixmax1,ixmin2:ixmax2)=f(ixmin1:ixmax1,&
         ixmin2:ixmax2)+v(ixmin1:ixmax1,ixmin2:ixmax2,idims)&
         *wCT(ixmin1:ixmax1,ixmin2:ixmax2,iw)
      ! Center flux to interface
      fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=half*(f(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+f(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
      if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=dxinv(idims)&
            *fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,idims)&
            -fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=-qdt*mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw,1)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
                 (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                    hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=-qdt*mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw,2)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
                 (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                    hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2)
         end select
      end if
   end do    !next iw
end do       !next idims

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
   .false.)

end subroutine centdiff
!=============================================================================
subroutine centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,w,fC,dx1,dx2,x)

! Advance the flow variables from t to t+qdt within ixO^L by
! fourth order centered differencing in space 
! for the dw/dt+dF_i(w)/dx_i=S type equation.
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1,dx2
double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)             ::  xi
double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim)

double precision :: v(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ndim), f(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)      :: vLC, vRC,&
    phi, cmaxLC, cmaxRC

double precision :: dxinv(1:ndim), dxdim
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
    hxCmin1,hxCmin2,hxCmax1,hxCmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2,&
    kkxCmin1,kkxCmin2,kkxCmax1,kkxCmax2, kkxRmin1,kkxRmin2,kkxRmax1,kkxRmax2
logical :: transport, new_cmax, patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------
! two extra layers are needed in each direction for which fluxes are added.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do

if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) then
   call mpistop("Error in CentDiff4: Non-conforming input limits")
end if

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;

if (useprimitive) then
   ! primitive limiting:
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif

! Add fluxes to w
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      end select
   end if

   ixmin1=ixOmin1-2*kr(idims,1);ixmin2=ixOmin2-2*kr(idims,2)
   ixmax1=ixOmax1+2*kr(idims,1);ixmax2=ixOmax2+2*kr(idims,2); 
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

   ixCmin1=hxOmin1;ixCmin2=hxOmin2; ixCmax1=ixOmax1;ixCmax2=ixOmax2;
   hxCmin1=ixCmin1-kr(idims,1);hxCmin2=ixCmin2-kr(idims,2)
   hxCmax1=ixCmax1-kr(idims,1);hxCmax2=ixCmax2-kr(idims,2); 
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2); 
   kxCmin1=ixCmin1+2*kr(idims,1);kxCmin2=ixCmin2+2*kr(idims,2)
   kxCmax1=ixCmax1+2*kr(idims,1);kxCmax2=ixCmax2+2*kr(idims,2); 

   kkxCmin1=ixImin1;kkxCmin2=ixImin2; kkxCmax1=ixImax1-kr(idims,1)
   kkxCmax2=ixImax2-kr(idims,2);
   kkxRmin1=kkxCmin1+kr(idims,1);kkxRmin2=kkxCmin2+kr(idims,2)
   kkxRmax1=kkxCmax1+kr(idims,1);kkxRmax2=kkxCmax2+kr(idims,2);
   wRC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,1:nwflux)=wCT(kkxRmin1:kkxRmax1,&
      kkxRmin2:kkxRmax2,1:nwflux)
   wLC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,1:nwflux)=wCT(kkxCmin1:kkxCmax1,&
      kkxCmin2:kkxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,1:ndim) = x(kkxCmin1:kkxCmax1,&
      kkxCmin2:kkxCmax2,1:ndim)
   xi(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,idims) = half* ( &
      x(kkxRmin1:kkxRmax1,kkxRmin2:kkxRmax2,idims)+x(kkxCmin1:kkxCmax1,&
      kkxCmin2:kkxCmax2,idims) )


   dxdim=-qdt/dxinv(idims)
   call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wCT,wCT,wLC,wRC,x,.false.,&
      dxdim)

   ! get auxiliaries for L and R states
   if (nwaux>0.and.(.not.(useprimitive))) then
         call getaux(.true.,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'cd4_wLC')
         call getaux(.true.,wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'cd4_wRC')
   end if

   ! Calculate velocities from upwinded values
   new_cmax=.true.
   call getcmax(new_cmax,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
      ixCmax1,ixCmax2,idims,cmaxLC,vLC,.false.)
   call getcmax(new_cmax,wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
      ixCmax1,ixCmax2,idims,cmaxRC,vLC,.false.)
   ! now take the maximum of left and right states
   vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cmaxRC(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

   ! Calculate velocities for centered values
   call getv(wCT,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idims,f)
   v(ixmin1:ixmax1,ixmin2:ixmax2,idims)=f(ixmin1:ixmax1,ixmin2:ixmax2)

if (useprimitive) then
   ! primitive limiting:
   ! this call ensures wCT is primitive with updated auxiliaries
   patchw(ixImin1:ixImax1,ixImin2:ixImax2) = .false.
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x,patchw)
endif

   do iw=1,nwflux
      ! Get non-transported flux
      call getflux(wCT,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
         ixmax1,ixmax2,iw,idims,f,transport)
      ! Add transport flux
      if (transport) f(ixmin1:ixmax1,ixmin2:ixmax2)=f(ixmin1:ixmax1,&
         ixmin2:ixmax2)+v(ixmin1:ixmax1,ixmin2:ixmax2,idims)&
         *wCT(ixmin1:ixmax1,ixmin2:ixmax2,iw)
      ! Center flux to interface
      ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
      fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=(-f(kxCmin1:kxCmax1,&
         kxCmin2:kxCmax2)+7.0d0*(f(jxCmin1:jxCmax1,jxCmin2:jxCmax2)&
         +f(ixCmin1:ixCmax1,ixCmin2:ixCmax2))-f(hxCmin1:hxCmax1,&
         hxCmin2:hxCmax2))/12.0d0
      ! add rempel dissipative flux, only second order version for now
      fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,iw,idims)-half*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) &
         *(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,iw))

      if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=dxinv(idims)&
            *fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)
         ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,idims)&
            -fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=-qdt*mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw,1)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
                 (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                    hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=-qdt*mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw,2)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
                 (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                    hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2)
         end select
      end if
   end do    !next iw
end do       !next idims

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,w,x,&
   .false.)

end subroutine centdiff4
!=============================================================================
! end module cd
!#############################################################################
