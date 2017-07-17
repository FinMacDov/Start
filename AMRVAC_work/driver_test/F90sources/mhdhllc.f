!##############################################################################
! module amrvacphys.mhdhllc.t 
!=============================================================================
subroutine diffuse_hllcd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idims,wLC,wRC,fLC,fRC,patchf)

! when method is hllcd or hllcd1 then: 

! this subroutine is to enforce regions where we AVOID HLLC
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    intent(in)      :: wRC,wLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
   intent(in) :: fLC, fRC

integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(inout)        :: patchf

integer                                           :: ixOO1,ixOO2,TxOOmin1,&
   TxOOmin2,TxOOmax1,TxOOmax2
!-----------------------------------

! In a user-controlled region around any point with flux sign change between
! left and right, ensure fallback to TVDLF
do ixOO1= ixOmin1,ixOmax1
do ixOO2= ixOmin2,ixOmax2
  
  TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idims,1), ixOmin1);
  TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idims,1), ixOmax1);
  
  
  TxOOmin2= max(ixOO2 - nxdiffusehllc*kr(idims,2), ixOmin2);
  TxOOmax2= min(ixOO2 + nxdiffusehllc*kr(idims,2), ixOmax2);
  
  if(abs(patchf(ixOO1,ixOO2)) == 1 .or. abs(patchf(ixOO1,ixOO2)) == 4)Then
     if(any(fRC(ixOO1,ixOO2,1:nwflux)*fLC(ixOO1,ixOO2,1:nwflux)<&
        -smalldouble))Then
       where(Abs(patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2))==1)
         patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2) = 4
       endwhere
     endif
  endif
enddo
enddo

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2, whll,Fhll,lambdaCD,patchf)

! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    intent(in)      :: wLC,wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(in):: fLC,fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: cmax,cmin

integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(inout)        :: patchf

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(out) :: Fhll,whll
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(out)            :: lambdaCD


logical         , dimension(ixImin1:ixImax1,ixImin2:ixImax2)     :: &
   Cond_patchf
double precision                       :: Epsilon
integer                                :: iw
!--------------------------------------------
!--------------------------------------------
! on entry, patch is preset to contain values from -2,1,2,4
!      -2: take left flux, no computation here
!      +2: take right flux, no computation here
!      +4: take TVDLF flux, no computation here
!       1: compute the characteristic speed for the CD

Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(abs(patchf(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2))==1)

do iw=1,nwflux
  where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  !============= compute HLL flux ==============!
  Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)= (cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
     -cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*fRC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,iw) + cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
     *cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*(wRC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)))&
     /(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2))
  !======== compute intermediate HLL state =======!
  whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = (cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
     -cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wLC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,iw)+fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)&
     -fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))/(cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  endwhere
enddo

! deduce the characteristic speed at the CD
where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=whll(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,m0_+idims)/whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
end where


where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
  ! double check whether obtained speed is in between min and max speeds given
  ! and identify in which part of the Riemann fan the time-axis is
  where(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero.and.lambdaCD&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2)>zero.and.lambdaCD(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)<cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -1
  elsewhere(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>zero.and.lambdaCD&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero.and.lambdaCD(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)>cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  1
  elsewhere(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2).or.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
     <= cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = zero
    ! we will fall back to HLL flux case in this degeneracy
    patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  3
  endwhere
endwhere


where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)== 3)
  Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.false.
end where


! handle the specific case where the time axis is exactly on the CD 
if(any(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==zero.and.Cond_patchf&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2)))then
  ! determine which sector (forward or backward) of the Riemann fan is smallest
  ! and select left or right flux accordingly
  where(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==zero.and.Cond_patchf&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    where(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=cmax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  1
    elsewhere
      patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -1
    endwhere
  endwhere
endif

!-----------------------------------------------------------------------!
! eigenvalue lambda for contact is near zero: decrease noise by this trick
if(flathllc)then
  Epsilon=1.0d-6
  where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2).and. dabs(lambdaCD&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2))/max(cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2),Epsilon)< Epsilon  .and. dabs(lambdaCD(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2))/max(dabs(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
     Epsilon)< Epsilon)
    lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) =  zero
  end where
end if 
!-----------------------------------------------------------------------!

return
end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
   ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! reference Li S., JCP, 203, 2005, 344-357
! reference T. Miyoski, Kusano JCP, 2008, 2005.

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    intent(in)      :: wRC,wLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(in):: whll, Fhll
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: vRC, vLC,lambdaCD
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: cmax,cmin
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(in):: fRC,fLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
   intent(out):: f

integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: patchf

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)        :: wCD,wSub
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)    :: fSub
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)             :: &
   vSub,cspeed,pCD,VdotBCD
integer                                        :: iw
!--------------------------------------------

!-------------- auxiliary Speed and array-------------!
where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == 1)
  cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cmax(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = vRC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
elsewhere(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == -1)
  cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cmin(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = vLC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
endwhere

do iw=1,nw
  if(iw /= b0_+idims)then
    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == 1)
      wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  wRC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw)
    elsewhere(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == -1)
      wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  wLC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw)
    endwhere
  endif
enddo

do iw=1,nwflux
  if(iw /= b0_+idims)then
    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == 1)
      fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  fRC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw)
    elsewhere(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == -1)
      fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) =  fLC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw)
    endwhere
  endif
enddo

where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
  wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = wSub(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,rho_)*(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
     -vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/(cspeed(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)-lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

 wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,Dtr1_)   = wSub(ixOmin1:ixOmax1,&
    ixOmin2:ixOmax2,Dtr1_)&
                   *(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                      -vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
                      /(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                      -lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) 

endwhere

!==== Magnetic field ====!
do iw =b1_,b3_
  where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
    ! case from eq 31
    wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = whll(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,iw)
  endwhere
enddo 

!------- Momentum ------!
do iw =1, 3
  if(iw /= idims)then
     where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
       ! eq. 21 22
      wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+iw)=(cspeed(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
         +iw)-fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+iw)  &
         -wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)*wCD(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,b0_+iw)) /(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
         -lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
     end where
  else
     where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
        ! eq. 20
        wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+iw) =  wCD(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,rho_) * lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
     endwhere
  endif
enddo

where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)

   VdotBCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (whll(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,m1_)*whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
      +whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)*whll(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b2_)+whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)&
      *whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_))/whll(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)
   ! Eq 17
   pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)  = wsub(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,rho_)*(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      -vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
                  *(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                     -vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
                  +fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                     +idims)- wsub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
                     +idims)*vSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                  + wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims)**2.0d0
   ! Eq 31
   wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = (cspeed(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2) * wSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) &
                 - fSub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) &
                    +lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                    *pCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                 -VdotBCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wCD&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims))&
                 /(cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-lambdaCD&
                    (ixOmin1:ixOmax1,ixOmin2:ixOmax2))

end where

do iw=1,nwflux
 if(iw /= b0_+idims)then
   where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2))==1)
       ! f_i=fsub+lambda (wCD-wSub)
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=fsub(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw)+cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
          *(wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)-wsub(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,iw))
   endwhere
 else
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero
 end if
end do

return
end subroutine getwCD
!=============================================================================
! end module amrvacphys.mhdhllc.t 
!##############################################################################
