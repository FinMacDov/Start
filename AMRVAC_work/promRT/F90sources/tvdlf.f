!============================================================================
subroutine PPMlimitervar(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idims,q,qCT,qLC,qRC)

! references:
! Mignone et al 2005, ApJS 160, 199,
! Miller and Colella 2002, JCP 183, 26
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,&
   ixmin2,ixmax1,ixmax2, idims
double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2),&
   qCT(ixImin1:ixImax1,ixImin2:ixImax2)

double precision, intent(inout) :: qRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
   qLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)  :: dqC,d2qC,ldq
double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2)  :: qMin,qMax,tmp

integer   :: lxCmin1,lxCmin2,lxCmax1,lxCmax2,lxRmin1,lxRmin2,lxRmax1,lxRmax2
integer   :: ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,ixLmax1,&
   ixLmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
   ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2
integer   :: hxLmin1,hxLmin2,hxLmax1,hxLmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
   hxRmin1,hxRmin2,hxRmax1,hxRmax2
integer   :: kxLLmin1,kxLLmin2,kxLLmax1,kxLLmax2,kxLmin1,kxLmin2,kxLmax1,&
   kxLmax2,kxCmin1,kxCmin2,kxCmax1,kxCmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
   kxRRmin1,kxRRmin2,kxRRmax1,kxRRmax2
!--------------------------------------------------------------------------
ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
ixOmax1=ixmax1+kr(idims,1);ixOmax2=ixmax2+kr(idims,2); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
ixLmax1=ixOmax1-kr(idims,1);ixLmax2=ixOmax2-kr(idims,2); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmin2=ixLmin2-kr(idims,2)
ixLLmax1=ixLmax1-kr(idims,1);ixLLmax2=ixLmax2-kr(idims,2); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
ixRmax1=ixOmax1+kr(idims,1);ixRmax2=ixOmax2+kr(idims,2); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmin2=ixRmin2+kr(idims,2)
ixRRmax1=ixRmax1+kr(idims,1);ixRRmax2=ixRmax2+kr(idims,2); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmax1=ixmax1;hxCmax2=ixmax2; !hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
hxLmax1=hxCmax1-kr(idims,1);hxLmax2=hxCmax2-kr(idims,2); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
hxRmax1=hxCmax1+kr(idims,1);hxRmax2=hxCmax2+kr(idims,2); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1;kxCmin2=ixLLmin2; kxCmax1=ixRmax1;kxCmax2=ixRmax2; !kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmin2=kxCmin2-kr(idims,2)
kxLmax1=kxCmax1-kr(idims,1);kxLmax2=kxCmax2-kr(idims,2); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmin2=ixLLmin2-kr(idims,2);lxCmax1=ixRRmax1
lxCmax2=ixRRmax2; !ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
lxRmax1=lxCmax1+kr(idims,1);lxRmax2=lxCmax2+kr(idims,2); !lxR=[iMmin1-3,ixMmax1+4]


dqC(lxCmin1:lxCmax1,lxCmin2:lxCmax2)=q(lxRmin1:lxRmax1,lxRmin2:lxRmax2)&
   -q(lxCmin1:lxCmax1,lxCmin2:lxCmax2)
! Eq. 64,  Miller and Colella 2002, JCP 183, 26
d2qC(kxCmin1:kxCmax1,kxCmin2:kxCmax2)=half*(q(kxRmin1:kxRmax1,&
   kxRmin2:kxRmax2)-q(kxLmin1:kxLmax1,kxLmin2:kxLmax2))
where(dqC(kxCmin1:kxCmax1,kxCmin2:kxCmax2)*dqC(kxLmin1:kxLmax1,&
   kxLmin2:kxLmax2)>zero)
   ! Store the sign of d2qC in qMin
   qMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2)= sign(one,d2qC(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2))
   ! Eq. 65,  Miller and Colella 2002, JCP 183, 26
   ldq(kxCmin1:kxCmax1,kxCmin2:kxCmax2)= qMin(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2)*min(dabs(d2qC(kxCmin1:kxCmax1,kxCmin2:kxCmax2)),2.0d0&
      *dabs(dqC(kxLmin1:kxLmax1,kxLmin2:kxLmax2)),2.0d0*dabs(dqC&
      (kxCmin1:kxCmax1,kxCmin2:kxCmax2)))
elsewhere
   ldq(kxCmin1:kxCmax1,kxCmin2:kxCmax2)=zero
end where

! Eq. 66, Miller and Colella 2002, JCP 183, 26
qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   +half*dqC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+(ldq(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)-ldq(ixRmin1:ixRmax1,ixRmin2:ixRmax2))/6.0d0
if (flatppm)then
   qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qRC(ixLLmin1:ixLLmax1,&
      ixLLmin2:ixLLmax2)+(half*dqC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)&
      +(ldq(ixLmin1:ixLmax1,ixLmin2:ixLmax2)-ldq(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))/6.0d0)
else
   qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2) &
      -(half*dqC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)+(ldq(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2)-ldq(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/6.0d0)
endif

! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
call extremaq(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,kxCmax1,kxCmax2,&
   qCT,1,qMax,qMin)

! Eq. B8, page 217, Mignone et al 2005, ApJS
qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=max(qMin(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2),min(qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
   qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)))
qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMin(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2),min(qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
   qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

! Eq. B9, page 217, Mignone et al 2005, ApJS
where((qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)-qCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2))*(qCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-qLC&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2))<=zero)
   qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end where

qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   -qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2))*(qCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)-(qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+qRC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2))/2.0d0)
qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   -qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2))**2.0d0/6.0d0
tmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2)=qRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2)

! Eq. B10, page 218, Mignone et al 2005, ApJS
where(qMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2)>qMin(hxRmin1:hxRmax1,&
   hxRmin2:hxRmax2))
   qRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2)= 3.0d0*qCT(hxRmin1:hxRmax1,&
      hxRmin2:hxRmax2)-2.0d0*qLC(hxRmin1:hxRmax1,hxRmin2:hxRmax2)
end where
! Eq. B11, page 218, Mignone et al 2005, ApJS
where(qMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2)<-qMin(hxCmin1:hxCmax1,&
   hxCmin2:hxCmax2))
   qLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2)= 3.0d0*qCT(hxCmin1:hxCmax1,&
      hxCmin2:hxCmax2)-2.0d0*tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2)
end where

end subroutine PPMlimitervar
!============================================================================
subroutine PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idims,w,wCT,wLC,wRC)

! references:
! Mignone et al 2005, ApJS 160, 199, 
! Miller and Colella 2002, JCP 183, 26 
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,&
   ixmin2,ixmax1,ixmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
   wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision, intent(inout) :: wRC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw),&
   wLC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw) 

double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)  :: dwC,d2wC,&
   ldw
double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nwflux)  :: wMin,&
   wMax,tmp
double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2) :: aa, ab, ac, dv
double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim) ::  exi

integer   :: lxCmin1,lxCmin2,lxCmax1,lxCmax2,lxRmin1,lxRmin2,lxRmax1,lxRmax2
integer   :: ixLLmin1,ixLLmin2,ixLLmax1,ixLLmax2,ixLmin1,ixLmin2,ixLmax1,&
   ixLmax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
   ixRRmin1,ixRRmin2,ixRRmax1,ixRRmax2
integer   :: hxLmin1,hxLmin2,hxLmax1,hxLmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
   hxRmin1,hxRmin2,hxRmax1,hxRmax2
integer   :: kxLLmin1,kxLLmin2,kxLLmax1,kxLLmax2,kxLmin1,kxLmin2,kxLmax1,&
   kxLmax2,kxCmin1,kxCmin2,kxCmax1,kxCmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
   kxRRmin1,kxRRmin2,kxRRmax1,kxRRmax2
integer   :: iw, idimss

double precision, parameter :: betamin=0.75d0, betamax=0.85d0,Zmin&
   =0.25d0, Zmax=0.75d0,eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0
!--------------------------------------------------------------------------
ixOmin1=ixmin1-kr(idims,1);ixOmin2=ixmin2-kr(idims,2)
ixOmax1=ixmax1+kr(idims,1);ixOmax2=ixmax2+kr(idims,2); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmin2=ixOmin2-kr(idims,2)
ixLmax1=ixOmax1-kr(idims,1);ixLmax2=ixOmax2-kr(idims,2); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmin2=ixLmin2-kr(idims,2)
ixLLmax1=ixLmax1-kr(idims,1);ixLLmax2=ixLmax2-kr(idims,2); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
ixRmax1=ixOmax1+kr(idims,1);ixRmax2=ixOmax2+kr(idims,2); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmin2=ixRmin2+kr(idims,2)
ixRRmax1=ixRmax1+kr(idims,1);ixRRmax2=ixRmax2+kr(idims,2); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmin2=ixOmin2;hxCmax1=ixmax1;hxCmax2=ixmax2; !hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmin2=hxCmin2-kr(idims,2)
hxLmax1=hxCmax1-kr(idims,1);hxLmax2=hxCmax2-kr(idims,2); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmin2=hxCmin2+kr(idims,2)
hxRmax1=hxCmax1+kr(idims,1);hxRmax2=hxCmax2+kr(idims,2); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1;kxCmin2=ixLLmin2; kxCmax1=ixRmax1;kxCmax2=ixRmax2; !kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmin2=kxCmin2-kr(idims,2)
kxLmax1=kxCmax1-kr(idims,1);kxLmax2=kxCmax2-kr(idims,2); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmin2=ixLLmin2-kr(idims,2);lxCmax1=ixRRmax1
lxCmax2=ixRRmax2; !ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmin2=lxCmin2+kr(idims,2)
lxRmax1=lxCmax1+kr(idims,1);lxRmax2=lxCmax2+kr(idims,2); !lxR=[iMmin1-3,ixMmax1+4]

dwC(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)=w(lxRmin1:lxRmax1,&
   lxRmin2:lxRmax2,1:nwflux)-w(lxCmin1:lxCmax1,lxCmin2:lxCmax2,1:nwflux)
! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
d2wC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=half*(w(kxRmin1:kxRmax1,&
   kxRmin2:kxRmax2,1:nwflux)-w(kxLmin1:kxLmax1,kxLmin2:kxLmax2,1:nwflux))
where(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)*dwC(kxLmin1:kxLmax1,&
   kxLmin2:kxLmax2,1:nwflux)>zero)
   ! Store the sign of dwC in wMin
   wMin(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)= sign(one,&
      d2wC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux))
   ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
   ldw(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)= wMin(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)*min(dabs(d2wC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
      1:nwflux)),2.0d0*dabs(dwC(kxLmin1:kxLmax1,kxLmin2:kxLmax2,1:nwflux)),&
      2.0d0*dabs(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)))
elsewhere
   ldw(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=zero
endwhere

! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux)+half*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   1:nwflux)+(ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)&
   -ldw(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nwflux))/6.0d0
if(flatppm)then 
   wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wRC(ixLLmin1:ixLLmax1,&
      ixLLmin2:ixLLmax2,1:nwflux)+half*dwC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
      1:nwflux)+(ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)&
      -ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))/6.0d0
else
   wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wRC(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,1:nwflux)-(half*dwC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
      1:nwflux)+(ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)&
      -ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))/6.0d0)
endif
 
! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
call extremaw(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,kxCmax1,kxCmax2,&
   wCT,1,wMax,wMin)

! Eq. B8, page 217, Mignone et al 2005, ApJS
wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=max(wMin(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux),min(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   1:nwflux),wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux))) 
wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=max(wMin(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux),min(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   1:nwflux),wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)))
  

! Eq. B9, page 217, Mignone et al 2005, ApJS
where((wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)-wCT(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux))*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)&
   -wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))<=zero)
   wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux)
   wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux)
end where

wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux))&
   *(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)-(wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux)+wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux))&
   /2.0d0)
wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=(wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,1:nwflux)-wRC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux))&
   **2.0d0/6.0d0
tmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nwflux)=wRC(ixLmin1:ixLmax1,&
   ixLmin2:ixLmax2,1:nwflux)
! Eq. B10, page 218, Mignone et al 2005, ApJS
where(wMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,1:nwflux)>wMin(hxRmin1:hxRmax1,&
   hxRmin2:hxRmax2,1:nwflux))
   wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,1:nwflux)= 3.0d0*wCT(hxRmin1:hxRmax1,&
      hxRmin2:hxRmax2,1:nwflux)-2.0d0*wLC(hxRmin1:hxRmax1,hxRmin2:hxRmax2,&
      1:nwflux)
endwhere
! Eq. B11, page 218, Mignone et al 2005, ApJS
where(wMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,1:nwflux)<-wMin(hxCmin1:hxCmax1,&
   hxCmin2:hxCmax2,1:nwflux))
   wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,1:nwflux)= 3.0d0*wCT(hxCmin1:hxCmax1,&
      hxCmin2:hxCmax2,1:nwflux)-2.0d0*tmp(hxLmin1:hxLmax1,hxLmin2:hxLmax2,&
      1:nwflux)
endwhere

! flattening at the contact discontinuities
if(flatcd)then
 call ppmflatcd(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,kxCmax1,&
    kxCmax2,kxLmin1,kxLmin2,kxLmax1,kxLmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
    wCT,d2wC,aa,ab)
 if(any(kappa*aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2)>=ab(kxCmin1:kxCmax1,&
    kxCmin2:kxCmax2)))then
  do iw=1,nwflux
   where(kappa*aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2)>=ab(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2).and. dabs(dwC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw))&
      >smalldouble)
      wMax(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw) = wCT(kxRmin1:kxRmax1,&
         kxRmin2:kxRmax2,iw)-2.0d0*wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw)&
         +wCT(kxLmin1:kxLmax1,kxLmin2:kxLmax2,iw)
   end where

   where(wMax(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)*wMax(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,iw)<zero .and.dabs(wCT(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
      iw)-wCT(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw))-eps*min(dabs(wCT&
      (ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)),dabs(wCT(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,iw)))>zero .and. kappa*aa(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)>=ab(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      .and. dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))>smalldouble)

     ac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixLLmin1:ixLLmax1,&
        ixLLmin2:ixLLmax2,iw)-wCT(ixRRmin1:ixRRmax1,ixRRmin2:ixRRmax2,iw)&
        +4.0d0*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))/(12.0d0&
        *dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw))
     wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=max(zero,min(eta1&
        *(ac(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-eta2),one))
   elsewhere
     wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=zero
   end where

   where(wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)>zero)
     wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = wLC(hxCmin1:hxCmax1,&
        hxCmin2:hxCmax2,iw)*(one-wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw))&
        +(wCT(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)+half*ldw(hxCmin1:hxCmax1,&
        hxCmin2:hxCmax2,iw))*wMin(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)
   end where
   where(wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)>zero)
     wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = wRC(hxCmin1:hxCmax1,&
        hxCmin2:hxCmax2,iw)*(one-wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw))&
        +(wCT(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)-half*ldw(hxRmin1:hxRmax1,&
        hxRmin2:hxRmax2,iw))*wMin(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)
   end where
  end do
 endif
endif

! flattening at the shocks
if(flatsh)then
  ! following MILLER and COLELLA 2002 JCP 183, 26
  kxCmin1=ixmin1-2;kxCmin2=ixmin2-2; kxCmax1=ixmax1+2;kxCmax2=ixmax2+2; !kxC=[ixMmin1-2,ixMmax1+2]
  do idimss=1,ndim
   kxLmin1=kxCmin1-kr(idimss,1);kxLmin2=kxCmin2-kr(idimss,2)
   kxLmax1=kxCmax1-kr(idimss,1);kxLmax2=kxCmax2-kr(idimss,2); !kxL=[ixMmin1-3,ixMmax1+1]
   kxRmin1=kxCmin1+kr(idimss,1);kxRmin2=kxCmin2+kr(idimss,2)
   kxRmax1=kxCmax1+kr(idimss,1);kxRmax2=kxCmax2+kr(idimss,2); !kxR=[ixMmin1-1,ixMmax1+3]
   kxLLmin1=kxLmin1-kr(idimss,1);kxLLmin2=kxLmin2-kr(idimss,2)
   kxLLmax1=kxLmax1-kr(idimss,1);kxLLmax2=kxLmax2-kr(idimss,2); !kxLL=[ixMmin-4,ixMmax]
   kxRRmin1=kxRmin1+kr(idimss,1);kxRRmin2=kxRmin2+kr(idimss,2)
   kxRRmax1=kxRmax1+kr(idimss,1);kxRRmax2=kxRmax2+kr(idimss,2); !kxRR=[ixMmin,ixMmax+4]

   call ppmflatsh(ixImin1,ixImin2,ixImax1,ixImax2,kxCmin1,kxCmin2,kxCmax1,&
      kxCmax2,kxLLmin1,kxLLmin2,kxLLmax1,kxLLmax2,kxLmin1,kxLmin2,kxLmax1,&
      kxLmax2,kxRmin1,kxRmin2,kxRmax1,kxRmax2,kxRRmin1,kxRRmin2,kxRRmax1,&
      kxRRmax2,idimss,wCT,aa,ab,dv)

   ! eq. B17, page 218, Mignone et al 2005, ApJS (had(Xi1))
   ac(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = max(zero,min(one,(betamax&
      -aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2))/(betamax-betamin)))
   ! eq. B18, page 218, Mignone et al 2005, ApJS (had(Xi1))
   ! recycling aa(ixL^S)
   where (dabs(dv(kxCmin1:kxCmax1,kxCmin2:kxCmax2))<smalldouble)
      aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = max(ac(kxCmin1:kxCmax1,&
         kxCmin2:kxCmax2), min(one,(Zmax-ab(kxCmin1:kxCmax1,kxCmin2:kxCmax2))&
         /(Zmax-Zmin)))
   elsewhere
      aa(kxCmin1:kxCmax1,kxCmin2:kxCmax2) = one
   endwhere

    
     call extremaa(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
        ixOmax2,aa,1,exi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,idimss))
  enddo
   ab(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(exi(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1),exi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
  ! recycling wMax
  do iw=1,nwflux
     where(dabs(ab(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-one)>smalldouble)
       wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw) = (one-ab(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)
     endwhere

     where(dabs(ab(hxCmin1:hxCmax1,hxCmin2:hxCmax2)-one)>smalldouble)
       wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = ab(hxCmin1:hxCmax1,&
          hxCmin2:hxCmax2)*wLC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)&
          +wMax(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)
     endwhere

     where(dabs(ab(hxRmin1:hxRmax1,hxRmin2:hxRmax2)-one)>smalldouble)
       wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw) = ab(hxRmin1:hxRmax1,&
          hxRmin2:hxRmax2)*wRC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iw)&
          +wMax(hxRmin1:hxRmax1,hxRmin2:hxRmax2,iw)
     endwhere
  enddo
endif

end subroutine PPMlimiter
!============================================================================
subroutine MP5limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,w,wCT,wLC,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
   iLmin2,iLmax1,iLmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
   wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
   wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2, iLmmmin1,&
   iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2, iLppmin1,&
   iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,iLpppmax2
integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
   idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
   idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
   iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
   ieppmax1,ieppmax2
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)  :: f, fmp,&
    fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)  :: wRCtmp,&
    wLCtmp
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: tmp, tmp2,&
    tmp3, a, b, c
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2)       :: flagL, flagR
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:
! range to process:
!iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

!{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! also add one ghost zone!
!   {iL^L=iL^L^LADD1;}
!}

! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);

f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 1.0d0/60.0d0 * (2.0d0&
   * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,1:nwflux) - 13.0d0&
   * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux) + 47.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iw)&
      -w(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
      -w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iw))
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw) &
      + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw) &
      + b(iLmin1:iLmax1,iLmin2:iLmax2)
end do ! iw loop

! get dm4:
idmax1=iLmax1;idmax2=iLmax2; idmin1=iLmin1-kr(idims,1)
idmin2=iLmin2-kr(idims,2);
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);

iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2); iemin1=idmin1
iemin2=idmin2;
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);

d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iepmin1:iepmax1,iepmin2:iepmax2,&
   1:nwflux)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,1:nwflux)+w(iemmin1:iemmax1,&
   iemmin2:iemmax2,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,iw)&
      -d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idpmin1:idpmax1,idpmin2:idpmax2,&
      iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp)
   a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp2)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
      tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))/2.0d0&
   -dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)) &
   + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = max(min(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),min(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = min(max(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),max(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
end do



! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-fmp(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux)) .le. eps)
   wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
elsewhere
   wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
end where

! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);

f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 1.0d0/60.0d0 * (2.0d0&
   * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,1:nwflux) - 13.0d0&
   * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux) + 47.0d0&
   * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
   1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
      -w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw)-w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iw))
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw) + b(iLmin1:iLmax1,iLmin2:iLmax2)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2); idmin1=iLmin1
idmin2=iLmin2;
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);

iemax1=idmax1;iemax2=idmax2; iemin1=idmin1-kr(idims,1)
iemin2=idmin2-kr(idims,2);
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
ieppmax1=iepmax1+kr(idims,1);ieppmax2=iepmax2+kr(idims,2);

d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iemin1:iemax1,iemin2:iemax2,&
   1:nwflux)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,1:nwflux)&
   +w(ieppmin1:ieppmax1,ieppmin2:ieppmax2,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,iw)&
      -d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmmin1:idmmax1,idmmin2:idmmax2,&
      iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp)
   a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp2)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
      tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))/2.0d0&
   -dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
   1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = max(min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = min(max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
end do

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)&
   -fmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))  .le. eps)
   wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = f(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
elsewhere
   wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
end where

! Since limiter not TVD, negative pressures or densities could result.  
! Fall back to flat interpolation (minmod would also work). 
call checkw(useprimitive,ixGlo1,ixGlo2,ixGhi1,ixGhi2,iLmin1,iLmin2,iLmax1,&
   iLmax2,wLCtmp,flagL)
call checkw(useprimitive,ixGlo1,ixGlo2,ixGhi1,ixGhi2,iLmin1,iLmin2,iLmax1,&
   iLmax2,wRCtmp,flagR)

do iw=1,nwflux
   where (flagL(iLmin1:iLmax1,iLmin2:iLmax2).and.flagR(iLmin1:iLmax1,&
      iLmin2:iLmax2))
      wLC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
         iw)
      wRC(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
         iw)
   end where
end do


end subroutine MP5limiter
!============================================================================
subroutine MP5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,w,wLC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
   iLmin2,iLmax1,iLmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2, iLmmmin1,&
   iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2, iLppmin1,&
   iLppmin2,iLppmax1,iLppmax2
integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
   idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
   idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
   iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
   ieppmax1,ieppmax2
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)  :: f, fmp,&
    fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: tmp, tmp2,&
    tmp3, a, b, c
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:


iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);

f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 1.0d0/60.0d0 * (2.0d0&
   * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,1:nwflux) - 13.0d0&
   * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux) + 47.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
   1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iw)&
      -w(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
      -w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iw))
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw) &
      + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw) &
      + b(iLmin1:iLmax1,iLmin2:iLmax2)
end do ! iw loop

! get dm4:
idmax1=iLmax1;idmax2=iLmax2; idmin1=iLmin1-kr(idims,1)
idmin2=iLmin2-kr(idims,2);
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);

iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2); iemin1=idmin1
iemin2=idmin2;
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);

d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iepmin1:iepmax1,iepmin2:iepmax2,&
   1:nwflux)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,1:nwflux)+w(iemmin1:iemmax1,&
   iemmin2:iemmax2,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,iw)&
      -d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idpmin1:idpmax1,idpmin2:idpmax2,&
      iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp)
   a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp2)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
      tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))/2.0d0&
   -dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)) &
   + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = max(min(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),min(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = min(max(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),max(w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
end do


! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-fmp(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux)) .le. eps)
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,&
      1:nwflux)
elsewhere
   wLC(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
end where


end subroutine MP5limiterL
!============================================================================
subroutine MP5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
   iLmax2,idims,w,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
   iLmin2,iLmax1,iLmax2, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
! .. local ..
integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2, iLpmin1,&
   iLpmin2,iLpmax1,iLpmax2, iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,&
   iLpppmin2,iLpppmax1,iLpppmax2
integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
   idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
   idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
   iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
   ieppmax1,ieppmax2
integer                         :: iw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)  :: f, fmp,&
    fmin, fmax, ful, dm4, d, fmd, flc, flim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: tmp, tmp2,&
    tmp3, a, b, c
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------
! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);

f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 1.0d0/60.0d0 * (2.0d0&
   * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,1:nwflux) - 13.0d0&
   * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux) + 47.0d0&
   * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,&
   iLmin2:iLmax2,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
   1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = w(iLmin1:iLmax1,iLmin2:iLmax2,iw)&
      -w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = alpha*(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw)-w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,iw))
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,tmp)
   fmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2)
   ful(iLmin1:iLmax1,iLmin2:iLmax2,iw) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
      iw) + b(iLmin1:iLmax1,iLmin2:iLmax2)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2); idmin1=iLmin1
idmin2=iLmin2;
idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);

iemax1=idmax1;iemax2=idmax2; iemin1=idmin1-kr(idims,1)
iemin2=idmin2-kr(idims,2);
iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
ieppmax1=iepmax1+kr(idims,1);ieppmax2=iepmax2+kr(idims,2);

d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iemin1:iemax1,iemin2:iemax2,&
   1:nwflux)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,1:nwflux)&
   +w(ieppmin1:ieppmax1,ieppmin2:ieppmax2,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,iw)&
      -d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmmin1:idmmax1,idmmin2:idmmax2,&
      iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp)
   a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
   b(idmin1:idmax1,idmin2:idmax2) = d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,a,&
      b,tmp2)
   call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,idmax2,&
      tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,idmin2:idmax2)
end do

! get fmd:
fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (w(iLmin1:iLmax1,iLmin2:iLmax2,&
   1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))/2.0d0&
   -dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
   1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)

fmin(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = max(min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),min(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

fmax(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = min(max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   fmd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)),max(w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux),ful(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
   flc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1,iLmin2:iLmax2) = fmin(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   b(iLmin1:iLmax1,iLmin2:iLmax2) = f(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   c(iLmin1:iLmax1,iLmin2:iLmax2) = fmax(iLmin1:iLmax1,iLmin2:iLmax2,iw)
   call median(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,iLmax2,a,&
      b,c,tmp)
   flim(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,iLmin2:iLmax2)
end do

! check case
where ((f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)-w(iLpmin1:iLpmax1,&
   iLpmin2:iLpmax2,1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)&
   -fmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))  .le. eps)
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,&
      1:nwflux)
elsewhere
   wRC(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flim(iLmin1:iLmax1,&
      iLmin2:iLmax2,1:nwflux)
end where

end subroutine Mp5limiterR
!============================================================================
subroutine minmod(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,a,b,minm)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
    b(ixImin1:ixImax1,ixImin2:ixImax2)
double precision, intent(out):: minm(ixImin1:ixImax1,ixImin2:ixImax2)

minm(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (sign(one,a(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2))+sign(one,b(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))&
   /2.0d0 * min(abs(a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),abs(b(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)))

end subroutine minmod
!============================================================================
subroutine median(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,a,b,c,med)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
    b(ixImin1:ixImax1,ixImin2:ixImax2), c(ixImin1:ixImax1,ixImin2:ixImax2)
double precision, intent(out):: med(ixImin1:ixImax1,ixImin2:ixImax2)
! .. local ..
double precision             :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
   tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = b(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
   - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = c(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
   - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

med(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = a(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
   + (sign(one,tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2))+sign(one,&
   tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))/2.0d0 * min(abs(tmp1&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2)),abs(tmp2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)))

end subroutine median
!============================================================================
!=============================================================================
subroutine hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax,qtC,wCT,qt,wnew,dx1,dx2,x)

! The non-conservative Hancock predictor for TVDLF

! on entry:
! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBdixB

! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2

! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2


! FCT not implemented here

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1,dx2, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: vLC, vRC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
logical :: transport
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw
!-----------------------------------------------------------------------------

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======Hancock predictor: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmin2=ixImin2;ixtestmax1=ixImax1
   ixtestmax2=ixImax2;
!   ixtestmin1=400;ixtestmax1=400+2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
endif

! Expand limits in each idims direction in which fluxes are added
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
   ixmax1=ixmax1+kr(idims,1);ixmax2=ixmax2+kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in Hancock: Nonconforming input limits")

if (useprimitive) then  
   if(oktest.and.mype==0) then
      do iw=1,nw 
         print *,'iw,wCT before useprimitive:',iw,wCT(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
      enddo
   endif
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 

if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wCT: before dimensional loop',iw,wCT(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
    enddo
endif

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
dxdim(1)=dx1;dxdim(2)=dx2;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      case (2)
         myB0 => myB0_face2
      end select
   end if

   ! Calculate w_j+g_j/2 and w_j-g_j/2
   ! First copy all variables, then upwind wLC and wRC.
   ! wLC is to the left of ixO, wRC is to the right of wCT.
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

   wRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux)
   wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wCT(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1:nwflux)

   if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wRC: before upwindLR',iw,wRC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
    enddo
    do iw=1,nw 
       print *,'iw,wLC: before upwindLR',iw,wLC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
    enddo
   endif

   call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,idims,wCT,wCT,wLC,wRC,x,.false.,&
      dxdim(idims))

   if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wRC: after upwindLR',iw,wRC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
    enddo
    do iw=1,nw 
       print *,'iw,wLC: after upwindLR',iw,wLC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
    enddo
   endif

   if (nwaux>0.and.(.not.(useprimitive))) then
   !!if (nwaux>0) then
      call getaux(.true.,wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,'hancock_wLC')
      call getaux(.true.,wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,&
         hxOmin2,hxOmax1,hxOmax2,'hancock_wRC')
   end if

   ! Calculate vLC and vRC velocities
   call getv(wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,hxOmin2,hxOmax1,&
      hxOmax2,idims,vRC)
   call getv(wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,idims,vLC)

   ! Advect w(iw)
   do iw=1,nwflux
      ! Calculate the fLC and fRC fluxes
      call getflux(wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,hxOmin2,&
         hxOmax1,hxOmax2,iw,idims,fRC,transport)
      call getflux(wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,iw,idims,fLC,transport)
      if (transport) then
         fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)=fRC(hxOmin1:hxOmax1,&
            hxOmin2:hxOmax2)+vRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
            *wRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw)
         fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=fLC(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+vLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
            *wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)
      end if

      if(oktest.and.mype==0) then
           print *,'iw: updating wnew',iw,wnew(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
           print *,'iw: updating wnew with fluxes fLC',fLC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2)
           print *,'iw: updating wnew with fluxes fRC',fRC(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2)
           if(.not.slab) then
             print *,'volumes:',mygeo%dvolume(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2)
             print *,'surfaceC1:',mygeo%surfaceC1(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2)
           endif
      endif

      ! Advect w(iw)
      if (slab) then
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+dxinv(idims)* (fLC(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)-fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))
      else
         select case (idims)
         case (1)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)-qdt/mygeo%dvolume(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) &
                  *(mygeo%surfaceC1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                     *fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
                   -mygeo%surfaceC1(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
                      *fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))
         case (2)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)-qdt/mygeo%dvolume(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2) &
                  *(mygeo%surfaceC2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
                     *fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2) &
                   -mygeo%surfaceC2(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
                      *fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))
         end select
      end if
   end do
end do ! next idims

if (useprimitive) then  
    patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
    call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
       ixImax2,wCT,x,patchw)
    if(oktest.and.mype==0) then
      do iw=1,nw 
         print *,'iw,wCT after conserve:',iw,wCT(ixtestmin1:ixtestmax1,ixtestmin2:ixtestmax2,iw)
      enddo
    endif
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixImin1,ixImin2,ixImax1,ixImax2,'hancock_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine hancock
!============================================================================
subroutine upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
   ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,idims,w,wCT,wLC,wRC,x,needprim,&
   dxdim)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixLmin1,ixLmin2,&
   ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, idims
logical, intent(in) :: needprim
double precision, intent(in) :: dxdim
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: w, wCT
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: x

integer :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
    jxCmin1,jxCmin2,jxCmax1,jxCmax2, iw, ixtestmin1,ixtestmin2,ixtestmax1,&
   ixtestmax2
double precision :: wLtmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    wRtmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision :: ldw(ixImin1:ixImax1,ixImin2:ixImax2), dwC(ixImin1:ixImax1,&
   ixImin2:ixImax2)
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: flagL, flagR
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw, patchwLC,&
    patchwRC

character*79 :: savetypelimiter
!-----------------------------------------------------------------------------

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======in upwindLR'
   ixtestmin1=ixImin1;ixtestmin2=ixImin2;ixtestmax1=ixImax1
   ixtestmax2=ixImax2;
!   ixtestmin1=400;ixtestmax1=400+2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
endif

! Transform w,wL,wR to primitive variables
if (needprim.and.useprimitive) then
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,w,x)
end if

if(typelimiter/='ppm' .and. typelimiter /= 'mp5')then
 jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
 jxRmax1=ixRmax1+kr(idims,1);jxRmax2=ixRmax2+kr(idims,2);
 ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=ixLmin1-kr(idims,1)
 ixCmin2=ixLmin2-kr(idims,2);
 jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
 jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

 do iw=1,nwflux
   if (loglimit(iw)) then
      w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=dlog10(w(ixCmin1:jxCmax1,&
         ixCmin2:jxCmax2,iw))
      wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=dlog10(wLC(ixLmin1:ixLmax1,&
         ixLmin2:ixLmax2,iw))
      wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=dlog10(wRC(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2,iw))
   end if

   dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,iw)&
      -w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)

   savetypelimiter=typelimiter
   if(savetypelimiter=='koren') typelimiter='korenL'
   if(savetypelimiter=='cada')  typelimiter='cadaL'
   if(savetypelimiter=='cada3') typelimiter='cada3L'
   call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
      ixCmax1,ixCmax2,iw,idims,ldw,dxdim)

   wLtmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=wLC(ixLmin1:ixLmax1,&
      ixLmin2:ixLmax2,iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2)
   if(savetypelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
        ixCmax1,ixCmax2,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
        ixCmax1,ixCmax2,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
        ixCmax1,ixCmax2,iw,idims,ldw,dxdim)
   endif
   wRtmp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=wRC(ixRmin1:ixRmax1,&
      ixRmin2:ixRmax2,iw)-half*ldw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)
   typelimiter=savetypelimiter

   if (loglimit(iw)) then
      w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=10.0d0**w(ixCmin1:jxCmax1,&
         ixCmin2:jxCmax2,iw)
      wLtmp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=10.0d0**wLtmp(ixLmin1:ixLmax1,&
         ixLmin2:ixLmax2,iw)
      wRtmp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=10.0d0**wRtmp(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2,iw)
   end if
 end do

 call checkw(useprimitive,ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
    ixLmax1,ixLmax2,wLtmp,flagL)
 call checkw(useprimitive,ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
    ixRmax1,ixRmax2,wRtmp,flagR)

 do iw=1,nwflux
   where (flagL(ixLmin1:ixLmax1,ixLmin2:ixLmax2).and.flagR(ixRmin1:ixRmax1,&
      ixRmin2:ixRmax2))
      wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=wLtmp(ixLmin1:ixLmax1,&
         ixLmin2:ixLmax2,iw)
      wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=wRtmp(ixRmin1:ixRmax1,&
         ixRmin2:ixRmax2,iw)
   end where

   if (loglimit(iw)) then
      where (.not.(flagL(ixLmin1:ixLmax1,ixLmin2:ixLmax2).and.flagR&
         (ixRmin1:ixRmax1,ixRmin2:ixRmax2)))
         wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=10.0d0**wLC(ixLmin1:ixLmax1,&
            ixLmin2:ixLmax2,iw)
         wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=10.0d0**wRC(ixRmin1:ixRmax1,&
            ixRmin2:ixRmax2,iw)
      end where
   end if
 enddo
else if (typelimiter .eq. 'ppm') then
 call PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,&
    idims,w,wCT,wLC,wRC)
else
 call MP5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
    ixLmax2,idims,w,wCT,wLC,wRC)
endif

! Transform w,wL,wR back to conservative variables
if (useprimitive) then
   if(needprim)then 
      patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
      call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
         ixImax2,w,x,patchw)
   endif
   patchwLC(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
   patchwRC(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
      ixLmax2,wLC,x,patchwLC)
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,ixRmax1,&
      ixRmax2,wRC,x,patchwRC)
end if

end subroutine upwindLR
!============================================================================
subroutine dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
   ixCmax1,ixCmax2,iw,idims,ldw,dxdim)

! Limit the centered dwC differences within ixC for iw in direction idim.
! The limiter is chosen according to typelimiter.

! Note that this subroutine is called from upwindLR (hence from methods 
! like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
! but also from the gradientS and divvectorS subroutines in geometry.t
! Accordingly, the typelimiter here corresponds to one of typelimiter1
! or one of typegradlimiter1.

! note: there is no iw dependence here...

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,ixCmin2,&
   ixCmax1,ixCmax2, iw, idims
double precision, intent(in) :: dxdim
double precision, intent(in) :: dwC(ixImin1:ixImax1,ixImin2:ixImax2)
double precision, intent(out) :: ldw(ixImin1:ixImax1,ixImin2:ixImax2)

double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)
integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2
double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12

! cada limiter parameter values
double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma=1.6d0
! full third order cada limiter
double precision :: rdelinv
double precision :: ldwA(ixImin1:ixImax1,ixImin2:ixImax2),ldwB&
   (ixImin1:ixImax1,ixImin2:ixImax2),tmpeta(ixImin1:ixImax1,ixImin2:ixImax2)
double precision, parameter :: cadepsilon=1.d-14, invcadepsilon&
   =1.d14,cadradius=0.1d0
!-----------------------------------------------------------------------------

! Contract indices in idim for output.
ixOmin1=ixCmin1+kr(idims,1);ixOmin2=ixCmin2+kr(idims,2); ixOmax1=ixCmax1
ixOmax2=ixCmax2;
hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

! Store the sign of dwC in tmp
tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sign(one,dwC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2))
rdelinv=one/(cadradius*dxdim)**2

select case (typelimiter)
case ('minmod')
   ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min(dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2)))
case ('woodward')
   ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=two*tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)* max(zero,min(dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2),tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*quarter*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
      +dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))))
case ('mcbeta')
   ! Woodward and Collela limiter, with factor beta
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min(mcbeta*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),mcbeta&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*half&
      *(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)+dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))))
case ('superbee')
   ! Roes superbee limiter (eq.3.51i)
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min(two*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2)),min(dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),two&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2)))
case ('vanleer')
  ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
  ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=two*max(dwC(hxOmin1:hxOmax1,&
     hxOmin2:hxOmax2)*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero) &
     /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+dwC(hxOmin1:hxOmax1,&
     hxOmin2:hxOmax2)+qsmall)
case ('albada')
  ! Albada limiter (eq.3.51g) with delta2=1D.-12
  ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
     *(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2+qsmall)+dwC(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)**2&
     +qsmall))/(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2+dwC(hxOmin1:hxOmax1,&
     hxOmin2:hxOmax2)**2+qsmall2)
case ('korenR')
   ! Barry Koren Right variant
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min(two*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),two&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2),(two*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)))*third))
case ('korenL')
   ! Barry Koren Left variant
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min(two*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),two&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2),(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+two*dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)))*third))
case ('cadaR')
   ! Cada Right variant
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min((two*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)))*third, max(-cadalfa*dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)),                     min(cadbeta*dabs(dwC&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2)),                  (two&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)*tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))&
      *third, cadgamma*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)))))
case ('cadaL')
   ! Cada Left variant
   ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min((two*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
      +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2))*third, max(-cadalfa*tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2),&
                           min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2),                  (two&
      *dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))+tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*third, cadgamma&
      *dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))))))
case ('cada3R')
   tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)**2)*rdelinv
   ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(two*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*third
   ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min((two*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)&
      *tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)))*third, max(-cadalfa*dabs(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)),                     min(cadbeta*dabs(dwC&
      (ixOmin1:ixOmax1,ixOmin2:ixOmax2)),                  (two&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)*tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))&
      *third, cadgamma*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)))))
   where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<=one-cadepsilon)
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=ldwA(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
   elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=one+cadepsilon)
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=ldwB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
   elsewhere
     tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmpeta(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-one)*invcadepsilon
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*( (one-tmp(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2) +(one&
        +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*ldwB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))
   endwhere
case ('cada3L')
   tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2)**2)*rdelinv
   ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(two*dwC(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*third
   ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      * max(zero,min((two*dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
      +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,&
      hxOmin2:hxOmax2))*third, max(-cadalfa*tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2),&
                           min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
      *dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2),                  (two&
      *dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))+tmp(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*third, cadgamma&
      *dabs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2))))))
   where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<=one-cadepsilon)
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=ldwA(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
   elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)>=one+cadepsilon)
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=ldwB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)
   elsewhere
     tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmpeta(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)-one)*invcadepsilon
     ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*( (one-tmp(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2) +(one&
        +tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*ldwB(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2))
   endwhere
case default
   write(*,*)'Unknown limiter:',typelimiter
   call mpistop("Error in dwLimiter: No such TVD limiter")
end select

end subroutine dwlimiter2
!=============================================================================
subroutine tvdlf(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

!> method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
!> method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.

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
   1:ndim)        :: fC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: fLC, fRC,&
    vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC,&
    cmaxRC, cmaxLC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
   ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
    ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
logical :: transport, new_cmax, CmaxMeanState
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======tvdlf: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmin2=ixImin2;ixtestmax1=ixImax1
   ixtestmax2=ixImax2;
   ixtestmin1=1;ixtestmax1=2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
endif

!!call mpistop("tijdelijke stop")

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='tvdlf1')call mpistop&
   ("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in TVDLF: Nonconforming input limits")


if ((method=='tvdlf').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 


dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
dxdim(1)=dx1;dxdim(2)=dx2;
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
!   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;


   ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;


   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmax1=ixCmax1;ixCRmax2=ixCmax2;

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:ndim) = x(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) = half* ( x(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) )


   ! for tvdlf (second order scheme): apply limiting
   if (method=='tvdlf') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wold,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wCT,&
            wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wnew,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
      ! determine mean state and store in wLC
      wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)= half*(wLC&
         (ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)+wRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1:nwflux))
      ! get auxilaries for mean state
      if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'tvdlf_cmaxmeanstate')
      end if
      new_cmax=.true.
      call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
         ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cmaxLC,.false.)
      
      ! We regain wLC for further use
      wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)=two*wLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1:nwflux)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1:nwflux)
      if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'tvdlf_wLC_A')
      endif
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'tvdlf_wRC_A')
      end if
   else
      ! get auxilaries for L and R states
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         !!if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'tvdlf_wLC')
         call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,'tvdlf_wRC')
      end if
      new_cmax=.true.
      call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
         ixCmin2,ixCmax1,ixCmax2,idims,cmaxLC,cmaxC,.false.)
      call getcmax(new_cmax,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
         ixCmin2,ixCmax1,ixCmax2,idims,cmaxRC,cmaxC,.false.)
      ! now take the maximum of left and right states
      cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cmaxRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
   end if
   
   ! Calculate velocities for transport fluxes
   call getv(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vLC)
   call getv(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vRC)



! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,idims)



   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call getflux(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,fLC,transport)
      call getflux(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,fRC,transport)
      if (transport) then
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*(fLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      ! Add TVDLF dissipation to the flux
      if ((.not.BnormLF) .and. (iw==b0_+idims .or.iw==psi_) .and. b0_>0) then
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.d0
      else
         ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=-tvdlfeps*cmaxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)*half*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)&
            -wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
      end if
      ! fLC contains physical+dissipative fluxes
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
         end select
      end if

   end do ! Next iw

end do ! Next idims



!Now update the state:
do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)=dxinv(idims)&
            *fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,1)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,2)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method=='tvdlf').and.useprimitive) then  
    patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
    call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
       ixImax2,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixImin1,ixImin2,ixImax1,ixImax2,'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine tvdlf
!=============================================================================
subroutine hll(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

!> method=='hll'  --> 2nd order HLL scheme.
!> method=='hll1' --> 1st order HLL scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1,dx2
integer, intent(in)                                  :: ixImin1,ixImin2,&
   ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
    intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:ndim)             :: xi
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nw)               :: wCT, wnew, wold
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
   1:ndim)  :: fC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: fLC, fRC,&
    vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC,&
    cmaxRC, cmaxLC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cminC,&
    cminRC, cminLC
double precision, dimension(1:ndim)     :: dxinv, dxdim
integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)               :: patchf
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
   ixCRmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2,&
    kxRmin1,kxRmin2,kxRmax1,kxRmax2
logical :: transport, new_cmax, CmaxMeanState, logiB
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hll1')call mpistop("Error in hll: Unsplit dim. and original is limited")


! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in hll : Nonconforming input limits")

if (method=='hll'.and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
dxdim(1)=dx1;dxdim(2)=dx2;
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


   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2); 
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2)
   

   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmax1=ixCmax1;ixCRmax2=ixCmax2;

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:ndim) = x(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) = half* ( x(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) )


   ! for hll (second order scheme): apply limiting
   if (method=='hll') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wold,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wCT,&
            wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wnew,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hll: no such base for limiter")
      end select
   end if

   ! For the high order hll scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)= half&
            *(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)&
            +wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hll_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)=two*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)-wRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,1:nwflux)
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hll_wLC_B')
         endif
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
            call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hll_wRC_B')
         end if
    else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
         !!if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hll_wLC')
            call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hll_wRC')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxLC,cminLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxRC,cminRC,.true.)
         ! now take the maximum of left and right states 
         ! S.F. Davis, SIAM J. Sci. Statist. Comput. 1988, 9, 445
         cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(cmaxRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
         cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=min(cminRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2),cminLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
   end if 

   patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
   where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
   elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
   endwhere

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB)) call &
      getv(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vLC)
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB)) call &
      getv(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vRC)


! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,idims)



   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) then 
        call getflux(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
           ixCmax1,ixCmax2,iw,idims,fLC,transport)
        if (transport) fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
           =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
     end if
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) then 
        call getflux(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
           ixCmax1,ixCmax2,iw,idims,fRC,transport)
        if (transport) fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
           =fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)+vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
     end if


     if (logiB.and.(iw==b0_+idims .or.iw==psi_)) then
        if (BnormLF) then
           ! flat B norm using tvdlf
           fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)= half*((fLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)) &
              -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
              dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)))&
              *(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2,iw)))
         else
           fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=zero
         endif
     else
       where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)
         ! Add hll dissipation to the flux
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(cmaxC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            -cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2) +tvdlfeps*cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*(wRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))&
            /(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-cminC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2))
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)== 2)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
       endwhere
     endif

     if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)
     else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
         end select
     end if

   end do ! Next iw
end do ! Next idims



do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)=dxinv(idims)&
            *fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,1)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,2)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if (method=='hll'.and.useprimitive) then  
   patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixImin1,ixImin2,ixImax1,ixImax2,'hll_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine hll
!=============================================================================
subroutine hllc(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

! method=='hllc'  --> 2nd order HLLC scheme.
! method=='hllc1' --> 1st order HLLC scheme.
! method=='hllcd' --> 2nd order HLLC+tvdlf scheme.
! method=='hllcd1'--> 1st order HLLC+tvdlf scheme.

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

double precision, dimension(1:ndim)                :: dxinv, dxdim

integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                          &
   :: patchf
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
   ixCRmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2,&
    kxRmin1,kxRmin2,kxRmax1,kxRmax2
logical :: transport, new_cmax, CmaxMeanState, logiB, firstordermethod
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw

!=== specific to HLLC and HLLCD ===!
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)     :: fLC, fRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux)     :: whll, Fhll, fCD
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)              :: &
   lambdaCD 
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)
logiB=.false.

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hllc1' .and. method/='hllcd1')call mpistop&
   ("Error in hllc: Unsplit dim. and original is limited")



! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in hllc : Nonconforming input limits")

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 
firstordermethod=(method=='hllc1'.or.method=='hllcd1')

dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
dxdim(1)=dx1;dxdim(2)=dx2;
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


   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

   ! enlarged for ppm purposes
   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2); 
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2)
   

   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmax1=ixCmax1;ixCRmax2=ixCmax2;

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:ndim) = x(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) = half* ( x(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) )


   ! for hllc and hllcd (second order schemes): apply limiting
   if (method=='hllc'.or.method=='hllcd') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wold,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wCT,&
            wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wnew,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hllc: no such base for limiter")
      end select
   end if

   ! For the high order hllc scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)= half&
            *(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)&
            +wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hllc_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)=two*wLC&
            (ixCmin1:ixCmax1,ixCmin2:ixCmax2,1:nwflux)-wRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,1:nwflux)
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hllc_wLC_B')
         end if
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
            call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hllc_wRC_B')
         end if
   else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
         !!if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hllc_wLC')
            call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
               ixCmin2,ixCmax1,ixCmax2,'hllc_wRC')
         end if
         new_cmax=.true.
         ! to save memory, use cmaxC and lambdaCD for cmacRC and cmaxLC respectively
         ! to save memory, use vLC   and vRC      for cminRC and cminLC respectively
         call getcmax(new_cmax,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,vLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,lambdaCD,vRC,.true.)
         ! now take the maximum of left and right states
         cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=max(lambdaCD(ixCmin1:ixCmax1,&
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

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB)) call &
      getv(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vLC)
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB)) call &
      getv(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vRC)


! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,idims)


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/= 2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) then
        call getfluxforhllc(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
           ixCmin2,ixCmax1,ixCmax2,iw,idims,fLC,transport)
        if (transport)  fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)&
           =fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)+vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
     end if
     if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=-2).or.(logiB&
        .and.(iw==b0_+idims .or.iw==psi_))) then
        call getfluxforhllc(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
           ixCmin2,ixCmax1,ixCmax2,iw,idims,fRC,transport)
        if (transport)   fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)&
           =fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)+vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
     end if
   end do

   ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4 
   if(method=='hllcd' .or. method=='hllcd1') call diffuse_hllcd(ixImin1,&
      ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wLC,wRC,&
      fLC,fRC,patchf)

   !---- calculate speed lambda at CD ----!
   if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)) call getlCD(wLC,wRC,&
      fLC,fRC,cminC,cmaxC,idims,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
      ixCmin2,ixCmax1,ixCmax2, whll,Fhll,lambdaCD,patchf)

   ! now patchf may be -1 or 1 due to getlCD 
   if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
      !======== flux at intermediate state ========!
      call getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cminC,&
         cmaxC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
         ixCmax2,idims,fCD)
   endif ! Calculate the CD flux

   do iw=1,nwflux
     if (logiB.and.(iw==b0_+idims .or.iw==psi_)) then
        if (BnormLF) then
           ! flat B norm using tvdlf
           fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw) = half*((fLC&
              (ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)+fRC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2,iw)) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2),dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)))&
              *(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,&
              ixCmin2:ixCmax2,iw)))
        else
           fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=zero
        end if
     else
       where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))==1)
        fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fCD(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
       elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==2)
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
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,iw)
     else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,iw)
         end select
      end if

   end do ! Next iw

end do ! Next idims





!Now update the state
do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)=dxinv(idims)&
            *fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,idims)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,1)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,2)=-qdt*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw) &
              + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then  
   patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
   call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixImin1,ixImin2,ixImax1,ixImax2,'hllc_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine hllc
!=============================================================================
subroutine tvdmusclf(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,idimmin,idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,dx2,x)

! method=='tvdmu'  --> 2nd order (may be 3rd order in 1D) TVD-MUSCL scheme.
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.)
! FCT not implemented here.
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
   1:ndim)        :: fC

double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC, wRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: fLC, fRC,&
    vLC, vRC
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC,&
    cmaxRC, cmaxLC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
   ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,kxRmax2,&
    ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
logical :: transport, new_cmax, CmaxMeanState
logical, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: patchw
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======tvdmusl: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmin2=ixImin2;ixtestmax1=ixImax1
   ixtestmax2=ixImax2;
   ixtestmin1=1;ixtestmax1=2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmin2,ixtestmax1,ixtestmax2
endif

!!call mpistop("tijdelijke stop")

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='tvdmu1')call mpistop&
   ("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
   ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
end do
if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call mpistop("Error in TVDMUSCLF: Nonconforming input limits")


if ((method=='tvdmu').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
      ixImax2,wCT,x)
endif 


dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
dxdim(1)=dx1;dxdim(2)=dx2;
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

   kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
   kxCmax2=ixImax2-kr(idims,2);
   kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
   kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

   ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmax1=ixCmax1;ixCRmax2=ixCmax2;

   wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,1:nwflux)
   wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nwflux)=wCT(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:ndim) = x(kxCmin1:kxCmax1,&
      kxCmin2:kxCmax2,1:ndim)
   xi(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) = half* ( x(kxRmin1:kxRmax1,&
      kxRmin2:kxRmax2,idims)+x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idims) )


   ! for tvdlf and tvdmu (second order schemes): apply limiting
   if (method=='tvdmu') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wold,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wCT,&
            wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
            ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wnew,&
            wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! handle all other methods than tvdlf, namely tvdmu and tvdmu1 here
   if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdmu1')) then
      !!if (nwaux>0) then
      call getaux(.true.,wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
         ixCmin2,ixCmax1,ixCmax2,'tvdlf_wLC_B')
      call getaux(.true.,wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
         ixCmin2,ixCmax1,ixCmax2,'tvdlf_wRC_B')
   end if
   

   ! Calculate velocities for transport fluxes
   call getv(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vLC)
   call getv(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
      ixCmax2,idims,vRC)



! Solve the Riemann problem for the linear 2x2 system for normal
! B-field and GLM_Psi according to Dedner 2002:
call glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
   ixCmax2,idims)


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call getflux(wLC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,fLC,transport)
      call getflux(wRC,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
         ixCmax1,ixCmax2,iw,idims,fRC,transport)
      if (transport) then
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
         fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=fRC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2)+vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
            *wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=half*(fLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      if (slab) then
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=dxinv(idims)&
            *fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,iw)+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
            idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,1)=-qdt*mygeo%surfaceC1&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,1)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,1))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         case (2)
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,2)=-qdt*mygeo%surfaceC2&
               (ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+ &
              (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,2)-fC(hxOmin1:hxOmax1,&
                 hxOmin2:hxOmax2,iw,2))/mygeo%dvolume(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)
         end select
      end if

   end do ! Next iw

   ! For the MUSCL scheme apply the characteristic based limiter
   if (method=='tvdmu'.or.method=='tvdmu1') call tvdlimit2(method,qdt,ixImin1,&
      ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,idims,wLC,wRC,wnew,x,fC,dx1,dx2)

end do ! Next idims

if ((method=='tvdmu').and.useprimitive) then  
    patchw(ixImin1:ixImax1,ixImin2:ixImax2)=.false.
    call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
       ixImax2,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,&
      ixImin1,ixImin2,ixImax1,ixImax2,'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImin2,ixImax1,&
   ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImin2,&
   ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
   .false.)

end subroutine tvdmusclf
!=============================================================================

subroutine glmSolve(wLC,wRC,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,idir)
include 'amrvacdef.f'
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
    intent(inout) :: wLC, wRC
integer, intent(in)                                  :: ixImin1,ixImin2,&
   ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)                  &
    :: dB, dPsi
!-----------------------------------------------------------------------------

! This implements eq. (42) in Dedner et al. 2002 JcP 175
! Gives the Riemann solution on the interface 
! for the normal B component and Psi in the GLM-MHD system.

! 23/04/2013 Oliver Porth

 !BL(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir)
 !BR(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir)

 !PsiL(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
 !PsiR(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)

dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   b0_+idir) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir)
dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   psi_) - wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)

wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir)   = half * (wRC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b0_+idir) + wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_&
   +idir)) - half/cmax_global * dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)       = half * (wRC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,psi_) + wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   psi_))         - half*cmax_global * dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idir) = wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b0_+idir)
wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLC(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,psi_)

end subroutine glmSolve
!=============================================================================

