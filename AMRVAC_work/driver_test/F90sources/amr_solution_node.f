!=============================================================================
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
include 'amrvacdef.f'

! getnode = get first available igrid on processor ipe

integer, intent(in) :: ipe

integer :: igrid, igrid_available
!----------------------------------------------------------------------------
igrid_available=0

do igrid=1,ngridshi
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available==0) then
   write(unitterm,*) " out of nodal space - allowed ",ngridshi," grids"
   call mpistop("")
else
   getnode=igrid_available
   igrid_inuse(igrid,ipe)=.true.
end if

if (ipe==mype) then
   ! initialize nodal block
   node(1:nodehi,getnode) = 0
   rnode(1:rnodehi,getnode) = zero
end if

end function getnode
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: level, ig1,ig2, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
    ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmax1,ixCoCoGmax2, ix
double precision :: rXmin1,rXmin2, dx1,dx2
!-----------------------------------------------------------------------------
ixCoGmin1=1;ixCoGmin2=1;
ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;

! initialize solution space
allocate(pwold(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw), &
   pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw), pwCoarse(igrid)%w&
   (ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw))
pwold(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=0.d0
pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=0.d0
pwCoarse(igrid)%w(ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:nw)=0.d0

if(residmin>smalldouble) allocate(pwres(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:nwflux))
if (errorestimate==1) then
   ixCoCoGmin1=1;ixCoCoGmin2=1;
   ixCoCoGmax1=ixCoGmax1/2+dixB;ixCoCoGmax2=ixCoGmax2/2+dixB;
   allocate(pwCoCo(igrid)%w(ixCoCoGmin1:ixCoCoGmax1,ixCoCoGmin2:ixCoCoGmax2,&
      1:nw))
end if

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig1=igrid_to_node(igrid,mype)%node%ig1;ig2=igrid_to_node(igrid,mype)%node%ig2;

node(plevel_,igrid)=level
node(pig1_,igrid)=ig1
node(pig2_,igrid)=ig2

! set dx information
rnode(rpdx1_,igrid)=dx(1,level)
rnode(rpdx2_,igrid)=dx(2,level)

! determine the minimal and maximal corners
rnode(rpxmin1_,igrid)=xprobmin1+dble(ig1-1)*dg1(level)
rnode(rpxmin2_,igrid)=xprobmin2+dble(ig2-1)*dg2(level)
rnode(rpxmax1_,igrid)=xprobmax1-dble(ng1(level)-ig1)*dg1(level)
rnode(rpxmax2_,igrid)=xprobmax2-dble(ng2(level)-ig2)*dg2(level)


allocate(px(igrid)%x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndim),pxCoarse(igrid)%x&
   (ixCoGmin1:ixCoGmax1,ixCoGmin2:ixCoGmax2,1:ndim))
dx1=rnode(rpdx1_,igrid)
dx2=rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-dixB*dx1
rXmin2=rnode(rpxmin2_,igrid)-dixB*dx2
do ix=ixGlo1,ixGhi1
   px(igrid)%x(ix,ixGlo2:ixGhi2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixGlo2,ixGhi2
   px(igrid)%x(ixGlo1:ixGhi1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do
dx1=2.0d0*rnode(rpdx1_,igrid)
dx2=2.0d0*rnode(rpdx2_,igrid)
rXmin1=rnode(rpxmin1_,igrid)-dixB*dx1
rXmin2=rnode(rpxmin2_,igrid)-dixB*dx2
do ix=ixCoGmin1,ixCoGmax1
   pxCoarse(igrid)%x(ix,ixCoGmin2:ixCoGmax2,1)=rXmin1+(dble(ix)-half)*dx1
end do
do ix=ixCoGmin2,ixCoGmax2
   pxCoarse(igrid)%x(ixCoGmin1:ixCoGmax1,ix,2)=rXmin2+(dble(ix)-half)*dx2
end do



! find the blocks on the boundaries
if (rnode(rpxmin1_,igrid)-dx(1,level)<xprobmin1.or.rnode(rpxmin2_,igrid)&
   -dx(2,level)<xprobmin2.or.rnode(rpxmax1_,igrid)+dx(1,level)>xprobmax1&
   .or.rnode(rpxmax2_,igrid)+dx(2,level)>xprobmax2) then
   phyboundblock(igrid)=.true.
else
   phyboundblock(igrid)=.false.
end if

if (.not.slab) call getgridgeo(igrid)

if (B0field) then
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
end if

end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

deallocate(pw(igrid)%w,pwold(igrid)%w,pwCoarse(igrid)%w)
if(residmin>smalldouble) deallocate(pwres(igrid)%w)
if (errorestimate==1) deallocate(pwCoCo(igrid)%w)
deallocate(px(igrid)%x,pxCoarse(igrid)%x)

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
