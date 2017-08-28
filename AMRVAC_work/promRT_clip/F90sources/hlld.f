subroutine hlld(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

! method=='hlld'  --> 2nd order HLLD scheme.
! method=='hlld1' --> 1st order HLLD scheme.
! method=='hlldd' --> 2nd order HLLD+tvdlf scheme.
! method=='hlldd1'--> 1st order HLLD+tvdlf scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1,dx2
integer, intent(in)                                  :: ixImin1,ixImin2,&
   ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
    intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)             ::  xi
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)               :: wCT, wnew, wold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
   1:ndim)  :: fC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)            :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                 &
   :: vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                 &
   :: cmaxC,cminC

double precision, dimension(1:ndim)                :: dxinv
double precision                                   :: dxdim

integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                          &
   :: patchf
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
    kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2
logical :: transport, new_cmax, CmaxMeanState, logiB
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw

!=== specific to HLLD and HLLDD ===!
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)     :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)     :: whll, Fhll, fD
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   -1:1)         :: lambdaD 
!-----------------------------------------------------------------------------

call mpistop('hlld is still just a dummy, not implemented...')

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hlld1' .and. method/='hlldd1')call mpistop&
   ("Error in hlld: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in hlld : Nonconforming input limits")

if ((method/='hlld'.and.method/='hlldd1').and.useprimitive) then  
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      end select
   end if

   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;

   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 (eq.4.38a,b)
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

   ! enlarged for ppm purposes
   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2); ![1,15]
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2); ![2,16]

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:ndim) = x(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) = half* ( x(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) )

   ! for hlld and hlldd (second order schemes): apply limiting
   if (method=='hlld'.or.method=='hlldd') then
      dxdim=-qdt/dxinv(idims)
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wold,wCT,&
            wLC,wRC,x,.true.,dxdim)
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wCT,wCT,wLC,&
            wRC,x,.false.,dxdim)
      case ('original')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
            ixCmax1,ixCmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wnew,wCT,&
            wLC,wRC,x,.true.,dxdim)
      case default
         call mpistop("Error in hlld: no such base for limiter")
      end select
   end if

   ! For the high order hlld scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)= half&
            *(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)&
            +wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hlld_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)=two*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)-wRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,1:nwflux)
   else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive))) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hlld_wLC')
            call getaux(.true.,wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hlld_wRC')
         end if
         new_cmax=.true.
         ! to save memory, use cmaxC and lambdaD for cmacRC and cmaxLC respectively
         ! to save memory, use vLC   and vRC      for cminRC and cminLC respectively
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,vLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cminC,vRC,.true.)
         ! now take the maximum of left and right states
         cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cminC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2),cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
         cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=min(vRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2),vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
   end if

   patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
   where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
   elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
   endwhere

   if ((nwaux>0) .and. CmaxMeanState .and.(.not.(useprimitive))) then
      call getaux(.true.,wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,'hlld_wLC_B')
      call getaux(.true.,wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,'hlld_wRC_B')
   end if

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB)) call &
      getv(wLC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
      idims,vLC)
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB)) call &
      getv(wRC,xi,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
      idims,vRC)


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) call getfluxforhllc(wLC,xi,ixGlo1,&
        ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idims,fLC,&
        transport)
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) call getfluxforhllc(wRC,xi,ixGlo1,&
        ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,iw,idims,fRC,&
        transport)
     if (transport) then
       if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB&
          .and.(iw==b0_+idims .or.iw==psi_))) fLC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)&
          +vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*wLC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)
       if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB&
          .and.(iw==b0_+idims .or. iw==psi_))) fRC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)&
          +vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*wRC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)
     end if
   end do
   ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4 
   if(method=='hlldd' .or. method=='hlldd1') call diffuse_hlldd(ixGlo1,ixGlo2,&
      ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wLC,wRC,fLC,fRC,&
      patchf)

   !---- calculate speed lambda at CD ----!
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)) call getlD(wLC,wRC,fLC,&
      fRC,cminC,cmaxC,idims,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,&
      ixCmax1,ixCmax2, whll,Fhll,lambdaD,patchf)

   ! now patchf may be -1 or 1 due to getlD 
   if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
      !======== flux at intermediate state ========!
      call getwD(wLC,wRC,whll,x,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaD,cminC,&
         cmaxC,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
         idims,fD)
   endif ! Calculate the CD flux


   do iw=1,nwflux
     if (logiB.and.(iw==b0_+idims .or.iw==psi_)) then
       fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw) = half*((fLC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)) &
          -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
          dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)))*(wRC(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
     else
       where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-3)
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))<=2)
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fD(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
        ! fallback option, reducing to HLL flux
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=Fhll(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==4)
        ! fallback option, reducing to TVDLF flux
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw) = half*((fLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)) &
           -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
           dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)))*(wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
       endwhere
     endif

     if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=dxinv(idims)&
            *fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
     else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=-qdt*mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=-qdt*mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method/='hlld1'.and.method/='hlldd1').and.useprimitive) then  
   patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x,patchw)
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
!if (sourceunsplit) call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
!                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine hlld
