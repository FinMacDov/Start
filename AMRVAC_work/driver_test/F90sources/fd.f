!=============================================================================
subroutine fd(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

include 'amrvacdef.f'

character(len=*), intent(in)                                     :: method
double precision, intent(in)                                     :: qdt, qtC,&
    qt, dx1,dx2
integer, intent(in)                                              :: ixImin1,&
   ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
    intent(in)            :: x

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    intent(inout)           :: wCT, wnew, wold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,1:ndim),&
    intent(out)  :: fC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)                      :: fCT
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)                          :: fm, fp, fmR, fpL
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                  &
                :: v
double precision                                                 :: &
   dxinv(1:ndim), dxdims
logical                                                          :: transport
integer                                                          :: idims, iw,&
    ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,&
   hxOmin2,hxOmax1,hxOmax2, ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2
!-----------------------------------------------------------------------------


dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
do idims= idimmin,idimmax

   select case (idims)
      case (1) 
      dxdims = dx1
      case (2) 
      dxdims = dx2
   end select
   if (B0field) then
      myB0 => myB0_cell
   end if

   ! Get fluxes for the whole grid (mesh+dixB)
    ixCmin1 = ixOmin1 - dixB * kr(idims,1)
     ixCmin2 = ixOmin2 - dixB * kr(idims,2)
    ixCmax1 = ixOmax1 + dixB * kr(idims,1)
     ixCmax2 = ixOmax2 + dixB * kr(idims,2)

   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixmax1=ixOmax1;ixmax2=ixOmax2; ixmin1=hxOmin1;ixmin2=hxOmin2;


   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmax1=ixCmax1;ixCRmax2=ixCmax2;


   ! Calculate velocities for transport fluxes
   call getv(wCT,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCRmin1,ixCRmin2,ixCRmax1,&
      ixCRmax2,idims,v)
   
   do iw=1,nwflux
      call getfluxforhllc(wCT,x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCRmin1,ixCRmin2,&
         ixCRmax1,ixCRmax2,iw,idims,fCT,transport)
      if (transport) fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) &
         = fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) + v(ixCRmin1:ixCRmax1,&
         ixCRmin2:ixCRmax2) * wCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw)
      ! Lax-Friedrich splitting:
      fp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) = half * (fCT&
         (ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) + tvdlfeps * cmax_global &
         * wCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw))
      fm(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) = half * (fCT&
         (ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw) - tvdlfeps * cmax_global &
         * wCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,iw))
   end do ! iw loop
  
   ! now do the reconstruction of fp and fm:
   call reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idims,fp,fpL,dxdims)
   call reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,idims,fm,fmR,dxdims)

   do iw=1,nwflux
      if (slab) then
         fC(ixmin1:ixmax1,ixmin2:ixmax2,iw,idims) = dxinv(idims) &
            * (fpL(ixmin1:ixmax1,ixmin2:ixmax2,iw) + fmR(ixmin1:ixmax1,&
            ixmin2:ixmax2,iw))
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixmin1:ixmax1,ixmin2:ixmax2,iw,1)=-qdt*mygeo%surfaceC1&
               (ixmin1:ixmax1,ixmin2:ixmax2) * (fpL(ixmin1:ixmax1,&
               ixmin2:ixmax2,iw) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,iw))
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixmin1:ixmax1,ixmin2:ixmax2,iw,2)=-qdt*mygeo%surfaceC2&
               (ixmin1:ixmax1,ixmin2:ixmax2) * (fpL(ixmin1:ixmax1,&
               ixmin2:ixmax2,iw) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,iw))
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if
   end do ! iw loop

end do !idims loop

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)


end subroutine fd
!=============================================================================
subroutine reconstructL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,w,wLC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
   iLmin2,iLmax1,iLmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    dxdims

double precision, intent(out)   :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 

double precision                :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
    dwC(ixImin1:ixImax1,ixImin2:ixImax2)
integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,&
   ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, kxCmin1,kxCmin2,&
   kxCmax1,kxCmax2, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
      iLmax2,idims,w,wLC)
case default 

   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2);

   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux) = w(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
   jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);

   ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
   ixCmin2=iLmin2-kr(idims,2);
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

   do iw=1,nwflux
      dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenL'
      if(savetypelimiter=='cada')  typelimiter='cadaL'
      if(savetypelimiter=='cada3') typelimiter='cada3L'
      call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
         +half*ldw(iLmin1:iLmax1,iLmin2:iLmax2)
   end do
end select

end subroutine reconstructL
!=============================================================================
subroutine reconstructR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,w,wRC,dxdims)
include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
   iLmin2,iLmax1,iLmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    dxdims

double precision, intent(out)   :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 

double precision                :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
    dwC(ixImin1:ixImax1,ixImin2:ixImax2)
integer                         :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,&
   ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, kxCmin1,kxCmin2,&
   kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2, iw
character*79                    :: savetypelimiter
!-----------------------------------------------------------------------------
select case (typelimiter)
case ('mp5')
   call MP5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
      iLmax2,idims,w,wRC)
case default 

   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=w(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)

   jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
   jxRmax1=iLmax1+kr(idims,1);jxRmax2=iLmax2+kr(idims,2);
   ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=iLmin1-kr(idims,1)
   ixCmin2=iLmin2-kr(idims,2);
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

   do iw=1,nwflux
      dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
      
      savetypelimiter=typelimiter
      if(savetypelimiter=='koren') typelimiter='korenR'
      if(savetypelimiter=='cada')  typelimiter='cadaR'
      if(savetypelimiter=='cada3') typelimiter='cada3R'
      call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,ldw,dxdims)
      typelimiter=savetypelimiter

      wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
         -half*ldw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)
   end do
end select

end subroutine reconstructR
!=============================================================================
