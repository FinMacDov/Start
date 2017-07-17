!============================================================================
subroutine eigenvectors(idims,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,leftm,righm,lamda)
! calculate eigenvectors and eigenvalues for flux Jacobian matrix of x_idims
! (Jiang, Chaowei et al. 2011, ApJ, 727, 101)
! Chun Xia 
! 13 May 2016
include 'amrvacdef.f'

integer, intent(in) :: idims,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
double precision :: leftm(ixImin1:ixImax1,ixImin2:ixImax2,1:nwwave,1:nwprim)
double precision :: righm(ixImin1:ixImax1,ixImin2:ixImax2,1:nwprim,1:nwwave)
double precision :: lamda(ixImin1:ixImax1,ixImin2:ixImax2,1:nwwave),&
    betad(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: a2,betau,&
   signbi,va,vf,vs,af,as,tmp
double precision :: sqrthalf
integer :: idirs, ii,is(ndir-1),si(ndir-1)

!----------------------------------------------------------------------------
sqrthalf=dsqrt(0.5d0)

a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=eqpar(gamma_)*w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,p_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)


betau(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt((w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)**2&
   +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2 )-w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b0_+idims)**2) 
where(betau(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==0)
  betau(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrthalf
end where
ii=0
do idirs=1,ndir
  if(idirs==idims) cycle 
  ii=ii+1
  is(ii)=idirs
  betad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirs)=w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,b0_+idirs)/betau(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end do
si(1)=is(2)
si(2)=is(1)
signbi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sign(1.d0,w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b0_+idims))
va(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dabs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   b0_+idims))/dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   + (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)**2+w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,b2_)**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b3_)**2 )&
   /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(vf(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)**2-4.d0*a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   *va(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2)
vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(0.5d0*(vf(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(0.5d0*(vf(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
af(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt((a2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)-vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2)&
   /(vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2-vs(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)**2))
as(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt((vf(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)**2-a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))&
   /(vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2-vs(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)**2))
leftm=0.d0
righm=0.d0

! entropy wave
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,entroW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_,entroW_)=1.d0
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,entroW_,rho_)=1.d0
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,entroW_,  p_)=-1.d0/a2(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)

! magnetic flux wave
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,diverW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idims,diverW_)=1.d0
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,diverW_,b0_+idims)=1.d0
! alfven wave down-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)+va(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),alfvRW_)= sqrthalf&
   *betad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,si(1))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),alfvRW_)=-sqrthalf&
   *betad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,si(2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),alfvRW_)=-sqrthalf&
   *betad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,si(1))*dsqrt(w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_))*signbi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),alfvRW_)= sqrthalf&
   *betad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,si(2))*dsqrt(w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_))*signbi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_,rho_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),alfvRW_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_,rho_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),alfvRW_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_, b0_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),alfvRW_)/w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_, b0_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),alfvRW_)/w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)
! alfven wave up-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvLW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)-va(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),alfvLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),alfvRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),alfvLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),alfvRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),alfvLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),alfvRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),alfvLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),alfvRW_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvLW_,rho_+is(1))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_,rho_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvLW_,rho_+is(2))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_,rho_+is(2))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvLW_, b0_+is(1))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_, b0_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvLW_, b0_+is(2))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,alfvRW_, b0_+is(2))
! fast wave down-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)+vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,fastRW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*af(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,fastRW_)=af(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)=-as(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*signbi&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),fastRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(1))

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,fastRW_)=af(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*eqpar(gamma_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)=as(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),fastRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+idims)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,fastRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),fastRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,  p_      )=af(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)/(2.d0*a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_, b0_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   rho_))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_, b0_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),fastRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   rho_))
! fast wave up-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)-vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,fastLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,fastRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,fastLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,fastRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),fastRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),fastLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),fastRW_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,fastLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,fastRW_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),fastRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),fastLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),fastRW_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_,rho_      )= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_      )
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_,rho_+idims)=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+idims)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_,rho_+is(1))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_,rho_+is(2))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,rho_+is(2))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_,  p_      )= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_,  p_      )

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_, b0_+is(1))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_, b0_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastLW_, b0_+is(2))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,fastRW_, b0_+is(2))
! slow wave down-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)+vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,slowRW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,rho_)*as(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,slowRW_)=as(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)=af(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*vf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*signbi&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),slowRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(1))

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,slowRW_)=as(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*eqpar(gamma_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)=-af(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)*dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),slowRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(2))
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)*betad(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+idims)=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,slowRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),slowRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,  p_      )=as(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)/(2.d0*a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
   *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_, b0_+is(1))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   rho_))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_, b0_+is(2))=righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),slowRW_)/(2.d0&
   *a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   rho_))
! slow wave up-stream
lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_)=w(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2,v0_+idims)-vs(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,slowLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_      ,slowRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,slowLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+idims,slowRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(1),slowRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),slowLW_)=-righm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_+is(2),slowRW_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,slowLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,  p_      ,slowRW_)

righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(1),slowRW_)
righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),slowLW_)= &
   righm(ixOmin1:ixOmax1,ixOmin2:ixOmax2, b0_+is(2),slowRW_)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_,rho_      )= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_      )
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_,rho_+idims)=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+idims)
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_,rho_+is(1))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_,rho_+is(2))=-leftm&
   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,rho_+is(2))

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_,  p_      )= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_,  p_      )

leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_, b0_+is(1))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_, b0_+is(1))
leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowLW_, b0_+is(2))= &
   leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,slowRW_, b0_+is(2))

end subroutine eigenvectors
!============================================================================
subroutine characteristic_project(idims,iside,ixImin1,ixImin2,ixImax1,ixImax2,&
   ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,dxndim,qdt)
! prepare time dependent boundary conditions using projected-characteristics
! method (Jiang, Chaowei et al. 2011, ApJ, 727, 101)
! Chun Xia
! 13 May 2016
include 'amrvacdef.f'

integer, intent(in) :: idims, iside, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),x(ixImin1:ixImax1,&
   ixImin2:ixImax2,1:ndim),dxndim(1:ndim),qdt
double precision :: jcobi(ixImin1:ixImax1,ixImin2:ixImax2,1:nwprim,1:nwprim)
double precision :: leftm(ixImin1:ixImax1,ixImin2:ixImax2,1:nwwave,1:nwprim)
double precision :: righm(ixImin1:ixImax1,ixImin2:ixImax2,1:nwprim,1:nwwave)
double precision :: ritmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nwprim,1:nwwave)
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwprim) :: Sx,&
    dw
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwwave) :: lamda, RHS
double precision :: invdx(ndir), dxall(ndir), sqrthalf
double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
integer :: idirs, ix1,ix2, ixRmin1,ixRmin2,ixRmax1,ixRmax2,iwave,iprim,ii,&
   is(ndir-1)
logical :: patchw(ixImin1:ixImax1,ixImin2:ixImax2)
!----------------------------------------------------------------------------
dxall(1:ndim)=dxndim(1:ndim)
if(ndir>ndim) dxall(ndim+1:ndir)=0.d0
invdx(1:ndim)=1.d0/dxall(1:ndim)
!print*,'cp ok 0'
call primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,&
   ixImax2,w,x)
!print*,'cp ok 1'
ii=0
do idirs=1,ndir
  if(idirs==idims) cycle
  ii=ii+1
  is(ii)=idirs
end do
Sx=0.d0
do ii=1,ndir-1
  idirs=is(ii)
  if(dxall(idirs)==0) cycle
  call eigenvectors(idirs,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,leftm,righm,lamda)
  ! negative part of flux Jacobian matrix
  ixRmin1=ixOmin1+kr(idirs,1);ixRmin2=ixOmin2+kr(idirs,2)
  ixRmax1=ixOmax1+kr(idirs,1);ixRmax2=ixOmax2+kr(idirs,2);
  dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)=(w(ixRmin1:ixRmax1,&
     ixRmin2:ixRmax2,1:nwprim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim))&
     *invdx(idirs)
  lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwwave)=0.5d0*(lamda&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwwave)-dabs(lamda(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:nwwave)))
  do iwave=1,nwwave
    do iprim=1,nwprim
      ritmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iprim,iwave)=lamda&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)*righm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iprim,iwave)
    end do
  end do
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
    jcobi(ix1,ix2,1:nwprim,1:nwprim)=matmul(ritmp(ix1,ix2,1:nwprim,1:nwwave),&
       leftm(ix1,ix2,1:nwwave,1:nwprim))
  end do
  end do
  do iprim=1,nwprim
    Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iprim)=Sx(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,iprim)-sum(jcobi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwprim,iprim)*dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim),2+1)
  end do
  ! positive part of flux Jacobian matrix
  ixRmin1=ixOmin1-kr(idirs,1);ixRmin2=ixOmin2-kr(idirs,2)
  ixRmax1=ixOmax1-kr(idirs,1);ixRmax2=ixOmax2-kr(idirs,2);
  dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)=(w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:nwprim)-w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nwprim))&
     *invdx(idirs)
  lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwwave)=0.5d0*(lamda&
     (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwwave)+dabs(lamda(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:nwwave)))
  do iwave=1,nwwave
    do iprim=1,nwprim
      ritmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iprim,iwave)=lamda&
         (ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)*righm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iprim,iwave)
    end do
  end do
  do ix2=ixOmin2,ixOmax2
  do ix1=ixOmin1,ixOmax1
    jcobi(ix1,ix2,1:nwprim,1:nwprim)=matmul(ritmp(ix1,ix2,1:nwprim,1:nwwave),&
       leftm(ix1,ix2,1:nwwave,1:nwprim))
  end do
  end do
  do iprim=1,nwprim
    Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iprim)=Sx(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,iprim)+sum(jcobi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nwprim,iprim)*dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim),2+1)
  end do
end do
if(iside==1) then
  ixRmin1=ixOmin1;ixRmin2=ixOmin2;ixRmax1=ixOmax1;ixRmax2=ixOmax2;
  ixRmin1=ixOmin1+kr(idims,1)
  ixRmin2=ixOmin2+kr(idims,2)
  call getdivb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,divb)
  select case (idims)
  case (1)
     divb(ixOmin1,ixOmin2:ixOmax2) = divb(ixOmin1+1,ixOmin2:ixOmax2) 
  case (2)
     divb(ixOmin1:ixOmax1,ixOmin2) = divb(ixOmin1:ixOmax1,ixOmin2+1) 
  end select
else
  ixRmin1=ixOmin1;ixRmin2=ixOmin2;ixRmax1=ixOmax1;ixRmax2=ixOmax2;
  ixRmax1=ixOmax1-kr(idims,1)
  ixRmax2=ixOmax2-kr(idims,2)
  call getdivb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,ixRmax1,&
     ixRmax2,divb)
  select case (idims)
  case (1)
     divb(ixOmax1,ixOmin2:ixOmax2) = divb(ixOmax1-1,ixOmin2:ixOmax2)
  case (2)
     divb(ixOmin1:ixOmax1,ixOmax2) = divb(ixOmin1:ixOmax1,ixOmax2-1)
  end select
end if
do idirs=1,ndir
   Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b0_+idirs)=Sx(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,b0_+idirs)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v0_&
      +idirs)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end do

call eigenvectors(idims,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,leftm,righm,lamda)
if(iside==1) then
  ixRmin1=ixOmin1+kr(idims,1);ixRmin2=ixOmin2+kr(idims,2)
  ixRmax1=ixOmax1+kr(idims,1);ixRmax2=ixOmax2+kr(idims,2);
  dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)=(w(ixRmin1:ixRmax1,&
     ixRmin2:ixRmax2,1:nwprim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim))&
     *invdx(idims)
else
  ixRmin1=ixOmin1-kr(idims,1);ixRmin2=ixOmin2-kr(idims,2)
  ixRmax1=ixOmax1-kr(idims,1);ixRmax2=ixOmax2-kr(idims,2);
  dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)=(w(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:nwprim)-w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nwprim))&
     *invdx(idims)
end if
do iwave=1,nwwave
  RHS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)=-lamda(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,iwave)*sum(leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave,&
     1:nwprim)*dw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim),2+1)&
     +sum(leftm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave,1:nwprim)&
     *Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim),2+1)
  where(lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)>=0.d0 .and. iside==1 &
     .or. lamda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)<=0.d0 .and. iside==2)
    RHS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwave)=0.d0
  end where
end do
do iprim=1,nwprim
  Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iprim)=sum(righm(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,iprim,1:nwwave)*RHS(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     1:nwwave),2+1)
end do
w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   1:nwprim)+qdt*Sx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwprim)
patchw=.false.
call conserve(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,ixImax1,ixImax2,&
   w,x,patchw)

end subroutine characteristic_project
!============================================================================
