!=============================================================================
subroutine refine_grid(child_igrid,child_ipe,igrid,ipe,active)

include 'amrvacdef.f'

integer, dimension(2,2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic1,ic2
!-----------------------------------------------------------------------------

! allocate solution space for new children
do ic2=1,2
do ic1=1,2
   call alloc_node(child_igrid(ic1,ic2))
end do
end do

if ((time_advance .and. active).or.convert.or.firstprocess) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   do ic2=1,2
   do ic1=1,2
      call initial_condition(child_igrid(ic1,ic2))
   end do
   end do
end if

! remove solution space of igrid
call dealloc_node(igrid)

end subroutine refine_grid
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)

include 'amrvacdef.f'

integer, dimension(2,2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe

integer :: ixmin1,ixmin2,ixmax1,ixmax2, ichild, ixComin1,ixComin2,ixComax1,&
   ixComax2, ic1,ic2
double precision :: dxCo1,dxCo2, xComin1,xComin2, dxFi1,dxFi2, xFimin1,xFimin2
!-----------------------------------------------------------------------------
if (typegridfill=="linear") then
   dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
   if (amrentropy) then

      ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

      call e_to_rhos(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
         pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then

      ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;

      call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
         pw(igrid)%w,px(igrid)%x)
   end if

   xComin1=rnode(rpxmin1_,igrid)
   xComin2=rnode(rpxmin2_,igrid)
   dxCo1=rnode(rpdx1_,igrid)
   dxCo2=rnode(rpdx2_,igrid)

end if

do ic2=1,2
do ic1=1,2
   ichild=child_igrid(ic1,ic2)

   ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
   ixComin2=ixMlo2+(ic2-1)*(ixMhi2-ixMlo2+1)/2
   ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2
   ixComax2=ixMhi2+(ic2-2)*(ixMhi2-ixMlo2+1)/2

   if (typegridfill=="linear") then
      xFimin1=rnode(rpxmin1_,ichild)
      xFimin2=rnode(rpxmin2_,ichild)
      dxFi1=rnode(rpdx1_,ichild)
      dxFi2=rnode(rpdx2_,ichild)

      call prolong_2nd(pw(igrid)%w,px(igrid)%x,ixComin1,ixComin2,ixComax1,&
         ixComax2,pw(ichild)%w,px(ichild)%x, &
                   dxCo1,dxCo2,xComin1,xComin2,dxFi1,dxFi2,xFimin1,xFimin2,&
                      ichild)

   else
      call prolong_1st(pw(igrid)%w,ixComin1,ixComin2,ixComax1,ixComax2,&
         pw(ichild)%w,px(ichild)%x)
   end if
end do
end do

if (typegridfill=="linear") then
   if (amrentropy) then
      call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
         pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixmin1,ixmin2,ixmax1,ixmax2,&
         pw(igrid)%w,px(igrid)%x,patchfalse)
   end if
end if

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2ab(wCo,xCo,ixComin1,ixComin2,ixComax1,ixComax2,wFi,xFi,&
   dxCo1,dxCo2,xComin1,xComin2,dxFi1,dxFi2,xFimin1,xFimin2,igridFi)
! interpolate children blocks including ghost cells

include 'amrvacdef.f'

integer, intent(in) :: ixComin1,ixComin2,ixComax1,ixComax2, igridFi
double precision, intent(in) :: dxCo1,dxCo2, xComin1,xComin2, dxFi1,dxFi2,&
    xFimin1,xFimin2
double precision, intent(in) :: wCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
    xCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim), xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:ndim)
double precision, intent(inout) :: wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer :: ixCo1,ixCo2, jxCo1,jxCo2, hxCo1,hxCo2, ixFi1,ixFi2, ix1,ix2, idim,&
    iw
integer :: ixFimin1,ixFimin2,ixFimax1,ixFimax2, ixCgmin1,ixCgmin2,ixCgmax1,&
   ixCgmax2, el
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo1,xCo2, xFi1,xFi2, eta1,eta2, invdxCo1,invdxCo2

!-----------------------------------------------------------------------------

invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;
el=ceiling(real(dixB)/2.)
ixCgmin1=ixComin1-el
ixCgmin2=ixComin2-el
ixCgmax1=ixComax1+el
ixCgmax2=ixComax2+el
do ixCo2 = ixCgmin2,ixCgmax2
   ! cell-centered coordinates of coarse grid point
   xCo2=xComin2+(dble(ixCo2-dixB)-half)*dxCo2

   ixFi2=2*(ixCo2-ixComin2)+ixMlo2
do ixCo1 = ixCgmin1,ixCgmax1
   ! cell-centered coordinates of coarse grid point
   xCo1=xComin1+(dble(ixCo1-dixB)-half)*dxCo1

   ixFi1=2*(ixCo1-ixComin1)+ixMlo1


   do idim=1,ndim
      hxCo1=ixCo1-kr(1,idim)
      hxCo2=ixCo2-kr(2,idim)
      jxCo1=ixCo1+kr(1,idim)
      jxCo2=ixCo2+kr(2,idim)

      do iw=1,nw
         slopeL=wCo(ixCo1,ixCo2,iw)-wCo(hxCo1,hxCo2,iw)
         slopeR=wCo(jxCo1,jxCo2,iw)-wCo(ixCo1,ixCo2,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), signR&
              *slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), mcbeta&
              *signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, (dabs(slopeR)&
              +two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,signC&
              *slopeR))
         end select
      end do
   end do
   do ix2=ixFi2,ixFi2+1
      if (ixFi2==0) cycle
      ! cell-centered coordinates of fine grid point
      xFi2=xFimin2+(dble(ix2-dixB)-half)*dxFi2
   do ix1=ixFi1,ixFi1+1
      if (ixFi1==0) cycle
      ! cell-centered coordinates of fine grid point
      xFi1=xFimin1+(dble(ix1-dixB)-half)*dxFi1


      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2;
      else
         eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2) &
            /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1,ix2))) 
         eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2) &
            /sum(pgeo(igridFi)%dvolume(ix1,ixFi2:ixFi2+1))) 

      end if

      wFi(ix1,ix2,1:nw) = wCo(ixCo1,ixCo2,1:nw) + (slope(1:nw,1)*eta1)&
         +(slope(1:nw,2)*eta2)
   end do
   end do
end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,wFi,&
      xFi)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,wFi,&
      xFi,patchfalse)
end if

end subroutine prolong_2ab
!=============================================================================
subroutine prolong_2nd(wCo,xCo,ixComin1,ixComin2,ixComax1,ixComax2,wFi,xFi,&
   dxCo1,dxCo2,xComin1,xComin2,dxFi1,dxFi2,xFimin1,xFimin2,igridFi)

include 'amrvacdef.f'

integer, intent(in) :: ixComin1,ixComin2,ixComax1,ixComax2, igridFi
double precision, intent(in) :: dxCo1,dxCo2, xComin1,xComin2, dxFi1,dxFi2,&
    xFimin1,xFimin2
double precision, intent(in) :: wCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
    xCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim), xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:ndim)
double precision, intent(inout) :: wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer :: ixCo1,ixCo2, jxCo1,jxCo2, hxCo1,hxCo2, ixFi1,ixFi2, ix1,ix2, idim,&
    iw
integer :: ixFimin1,ixFimin2,ixFimax1,ixFimax2
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo1,xCo2, xFi1,xFi2, eta1,eta2, invdxCo1,invdxCo2

!-----------------------------------------------------------------------------

invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;
do ixCo2 = ixComin2,ixComax2
   ! cell-centered coordinates of coarse grid point
   xCo2=xComin2+(dble(ixCo2-dixB)-half)*dxCo2

   ixFi2=2*(ixCo2-ixComin2)+ixMlo2
do ixCo1 = ixComin1,ixComax1
   ! cell-centered coordinates of coarse grid point
   xCo1=xComin1+(dble(ixCo1-dixB)-half)*dxCo1

   ixFi1=2*(ixCo1-ixComin1)+ixMlo1


   do idim=1,ndim
      hxCo1=ixCo1-kr(1,idim)
      hxCo2=ixCo2-kr(2,idim)
      jxCo1=ixCo1+kr(1,idim)
      jxCo2=ixCo2+kr(2,idim)

      do iw=1,nw
         slopeL=wCo(ixCo1,ixCo2,iw)-wCo(hxCo1,hxCo2,iw)
         slopeR=wCo(jxCo1,jxCo2,iw)-wCo(ixCo1,ixCo2,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), signR&
              *slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), mcbeta&
              *signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, (dabs(slopeR)&
              +two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,signC&
              *slopeR))
         end select
      end do
   end do
   do ix2=ixFi2,ixFi2+1
      ! cell-centered coordinates of fine grid point
      xFi2=xFimin2+(dble(ix2-dixB)-half)*dxFi2
   do ix1=ixFi1,ixFi1+1
      ! cell-centered coordinates of fine grid point
      xFi1=xFimin1+(dble(ix1-dixB)-half)*dxFi1


      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2;
      else
         eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2) &
            /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1,ix2))) 
         eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeo(igridFi)%dvolume(ix1,ix2) &
            /sum(pgeo(igridFi)%dvolume(ix1,ixFi2:ixFi2+1))) 

      end if

      wFi(ix1,ix2,1:nw) = wCo(ixCo1,ixCo2,1:nw) + (slope(1:nw,1)*eta1)&
         +(slope(1:nw,2)*eta2)
   end do
   end do
end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,wFi,&
      xFi)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,wFi,&
      xFi,patchfalse)
end if

end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(wCo,ixComin1,ixComin2,ixComax1,ixComax2,wFi,xFi)

include 'amrvacdef.f'

integer, intent(in) :: ixComin1,ixComin2,ixComax1,ixComax2
double precision, intent(in) :: wCo(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
    xFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim)
double precision, intent(out) :: wFi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

integer :: ixCo1,ixCo2, ixFi1,ixFi2, iw
integer :: ixFimin1,ixFimin2,ixFimax1,ixFimax2
!-----------------------------------------------------------------------------
do ixCo2 = ixComin2,ixComax2
   ixFi2=2*(ixCo2-ixComin2)+ixMlo2
do ixCo1 = ixComin1,ixComax1
   ixFi1=2*(ixCo1-ixComin1)+ixMlo1
   forall(iw=1:nw) wFi(ixFi1:ixFi1+1,ixFi2:ixFi2+1,iw)=wCo(ixCo1,ixCo2,iw)
end do
end do

end subroutine prolong_1st
!=============================================================================
