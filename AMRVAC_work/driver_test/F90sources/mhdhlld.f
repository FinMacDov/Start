!##############################################################################
! module amrvacphys.mhdhlld.t 
!=============================================================================
! Aim : HLLD
! Description :
!=============================================================================
subroutine diffuse_hlldd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
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
call mpistop('difuse_hlldd is just a dummy for hlld, not yet implemented')

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

end subroutine diffuse_hlldd
!=============================================================================
subroutine getlD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2, whll,Fhll,lambdaD,patchf)

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
    intent(inout) :: Fhll,whll
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   -1:1), intent(inout)     :: lambdaD


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

call mpistop('getlD is just a dummy for hlld, not yet implemented')
return
end subroutine getlD
!=============================================================================
subroutine getwD(wLC,wRC,whll,x,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaD,cmin,cmax,&
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
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
    intent(in)    :: x
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(in):: whll, Fhll
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: vRC, vLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   -1:1), intent(in)      ::lambdaD
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: cmax,cmin
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
    intent(in)  :: fRC,fLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux),&
   intent(inout)  :: f
integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2),&
    intent(in)           :: patchf

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)        :: wCD,wD,wLD,wRD,wSub
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)    :: fSub
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)             :: &
   vSub,cspeed,pCD,ptotLC,ptotRC,VdotBCD,SignB
integer                                        :: iw
!--------------------------------------------
call mpistop('getwD is just a dummy for hlld, not yet implemented')
return
end subroutine getwD
!=============================================================================
! end module amrvacphys.mhdhlld.t 
!##############################################################################
