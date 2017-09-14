!##############################################################################
! module mhd/correctaux.t version July 2009 
! merged with mhdiso, September 2012 by Oliver Porth
!=============================================================================
subroutine correctaux(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,x,patchierror,subname)

include 'amrvacdef.f'

integer, intent(in)         :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
integer, intent(in)         :: patchierror(ixImin1:ixImax1,ixImin2:ixImax2)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in)      :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
integer        :: iw, kxOmin1,kxOmin2,kxOmax1,kxOmax2, ix1,ix2, i
logical        :: patchw(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

do ix2= ixOmin2,ixOmax2
do ix1= ixOmin1,ixOmax1
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getpthermal/total or primitive call
    if (patchierror(ix1,ix2)/=0) then
        ! verify in cube with border width nflatgetaux the presence
        ! of cells where all went ok
        do i=1,nflatgetaux
           kxOmin1= max(ix1-i,ixOmin1);
            kxOmax1= min(ix1+i,ixOmax1);
           kxOmin2= max(ix2-i,ixOmin2);
            kxOmax2= min(ix2+i,ixOmax2);
           ! in case cells are fine within smaller cube than 
           ! the userset nflatgetaux: use that smaller cube
           if (any(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2)==0)) exit
        end do
        if (any(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2)==0))then
           ! within surrounding cube, cells without problem were found
           if (subname/='primitive') then
              patchw(kxOmin1:kxOmax1,kxOmin2:kxOmax2)=(patchierror&
                 (kxOmin1:kxOmax1,kxOmin2:kxOmax2)/=0)
              call primitiven(ixImin1,ixImin2,ixImax1,ixImax2,kxOmin1,kxOmin2,&
                 kxOmax1,kxOmax2,w,patchw)
           endif
           ! faulty cells are corrected by averaging here
           ! only average those which were ok and replace faulty cells
           do iw = 1,nw
              if (iw/=b1_.and.iw/=b2_.and.iw/=b3_) then
                   w(ix1,ix2,iw)=sum(w(kxOmin1:kxOmax1,kxOmin2:kxOmax2,iw),&
                      patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2)==0)&
                      /count(patchierror(kxOmin1:kxOmax1,kxOmin2:kxOmax2)==0)
              end if
           end do
           if (subname/='primitive') then
              ! in addition to those switched to primitive variables
              ! above, also switch the corrected variables
              patchw(ix1,ix2)=.false.
              call conserven(ixImin1,ixImin2,ixImax1,ixImax2,kxOmin1,kxOmin2,&
                 kxOmax1,kxOmax2,w,patchw)
           endif
        else
           ! no cells without error were found in cube of size nflatgetaux
           ! --> point of no recovery
           print*,'Getaux error:',patchierror(ix1,ix2),'ix^D=',ix1,ix2
           print*,'New ','rho=',w(ix1,ix2,rho_),'m=', w(ix1,ix2,m1_),w(ix1,ix2,m2_),w(ix1,ix2,m3_)
           print*,'B=',w(ix1,ix2,b1_),w(ix1,ix2,b2_),w(ix1,ix2,b3_)

           print*,'pressure ', (eqpar(gamma_)-one)*(w(ix1,ix2,e_)- &
                           half*((w(ix1,ix2,m1_)**2+w(ix1,ix2,m2_)**2&
                              +w(ix1,ix2,m3_)**2)/w(ix1,ix2,rho_)&
                                 +w(ix1,ix2,b1_)**2+w(ix1,ix2,b2_)**2&
                                    +w(ix1,ix2,b3_)**2))
           print*,'e=',w(ix1,ix2,e_)

           print*,'position ', x(ix1,ix2, 1:ndim),' time ',t,it
           print*,'Called from: ',subname
           if (patchierror(ix1,ix2)<0) then
               call mpistop("-correctaux from smallvalues-----")
           else
               call mpistop("-correctaux from primitive-getpthermal-total-")
           end if
        end if
    end if
enddo
enddo

end subroutine correctaux
!=============================================================================
subroutine smallvalues(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,subname)

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
character(len=*), intent(in)    ::subname
!.. local ..
integer                         :: posvec(ndim)
integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)       :: patchierror
double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
    Te(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------


pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(eqpar(gamma_)-one)*(w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,e_)- &
       half*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
       + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+ w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,b2_)**2+ w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2))
if(smallT>0.d0) then
  Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
     /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
else
  Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
end if
if(strictsmall) then
  if(smallT>0.d0 .and. any(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2) <=smallT)) then
    print *,'SMALLVALUES of temperature under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    posvec(1)=posvec(1)+ixOmin1-1;posvec(2)=posvec(2)+ixOmin2-1;
    write(*,*)'minimum temperature= ', minval(Te(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),' with limit=',smallT,&
    ' at x=', x(posvec(1),posvec(2),1:ndim),' array index=',posvec,&
       ' where rho=',&
    w(posvec(1),posvec(2),rho_),', velocity v=',&
    w(posvec(1),posvec(2),m1_)/w(posvec(1),posvec(2),rho_),w(posvec(1),&
       posvec(2),m2_)/w(posvec(1),posvec(2),rho_),w(posvec(1),posvec(2),m3_)&
       /w(posvec(1),posvec(2),rho_),&
    ', and magnetic field B=',w(posvec(1),posvec(2),b1_),w(posvec(1),&
       posvec(2),b2_),w(posvec(1),posvec(2),b3_),&
    ' w(1:nwflux)=',w(posvec(1),posvec(2),1:nwflux)
    call mpistop("Smallvalues of temperature with strictsmall=T failed")
  endif
  if(any(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) <=minp)) then
    print *,'SMALLVALUES of pressure under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    posvec(1)=posvec(1)+ixOmin1-1;posvec(2)=posvec(2)+ixOmin2-1;
    write(*,*)'minimum pressure = ', minval(pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),' with limit=',minp,&
    ' at x=', x(posvec(1),posvec(2),1:ndim),' array index=',posvec,&
       ' where rho=',&
    w(posvec(1),posvec(2),rho_),', velocity v=',&
    w(posvec(1),posvec(2),m1_)/w(posvec(1),posvec(2),rho_),w(posvec(1),&
       posvec(2),m2_)/w(posvec(1),posvec(2),rho_),w(posvec(1),posvec(2),m3_)&
       /w(posvec(1),posvec(2),rho_),&
    ', and magnetic field B=',w(posvec(1),posvec(2),b1_),w(posvec(1),&
       posvec(2),b2_),w(posvec(1),posvec(2),b3_),&
    ' w(1:nwflux)=',w(posvec(1),posvec(2),1:nwflux)
    call mpistop("Smallvalues of pressure with strictsmall=T failed")
  endif
  if(any(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) <=smalle)) then
    print *,'SMALLVALUES of energy under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
    posvec(1)=posvec(1)+ixOmin1-1;posvec(2)=posvec(2)+ixOmin2-1;
    write(*,*)'minimum e =', minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)),&
       ' with limit=',smalle,&
    ' at x=', x(posvec(1),posvec(2),1:ndim),' array index=',posvec,&
       ' where E_k=',&
    half*(w(posvec(1),posvec(2),m1_)**2+w(posvec(1),posvec(2),m2_)**2&
       +w(posvec(1),posvec(2),m3_)**2)/w(posvec(1),posvec(2),rho_),&
    ' E_total=',w(posvec(1),posvec(2),e_),&
    ' w(1:nwflux)=',w(posvec(1),posvec(2),1:nwflux)
    call mpistop("Smallvalues of energy with strictsmall=T failed")
  endif
  if(any(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) <=minrho)) then
    print *,'SMALLVALUES of density under strictsmall problem From:  ', &
            subname,' iteration=', it,' time=',t
    posvec(1:ndim)=minloc(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
    posvec(1)=posvec(1)+ixOmin1-1;posvec(2)=posvec(2)+ixOmin2-1;
    write(*,*)'minimum rho =', minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)),' with limit=',minrho,&
    ' at x=', x(posvec(1),posvec(2),1:ndim),' array index=',posvec,&
       ' where E_k=',&
    half*(w(posvec(1),posvec(2),m1_)**2+w(posvec(1),posvec(2),m2_)**2&
       +w(posvec(1),posvec(2),m3_)**2)/w(posvec(1),posvec(2),rho_),&
    ' E_total=',w(posvec(1),posvec(2),e_),&
    ' w(1:nwflux)=',w(posvec(1),posvec(2),1:nwflux)
    call mpistop("Smallvalues of density with strictsmall=T failed")
  endif
else
  if(strictgetaux)then
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) < minrho)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=minrho
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_) =zero
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_) =zero
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_) =zero;
    end where
    where(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < minp)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=minp/(eqpar(gamma_)-one)+&
       ((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
          +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2)*half
    end where
    where(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < smallT)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=smallT*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)/(eqpar(gamma_)-one)+&
       ((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)**2+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,m2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m3_)**2)&
          /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
          +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2)*half
    end where
  else
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) < minrho &
       .or. w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) < smalle&
          .or. pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < minp &
             .or. Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < smallT)
      patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-1
    elsewhere
      patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0
    end where
    call correctaux(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,patchierror,subname)
  end if
end if



end subroutine smallvalues
!=============================================================================
! end module mhd/correctaux.t
!##############################################################################
