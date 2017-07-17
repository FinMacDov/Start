!=============================================================================
!> Initialize (and allocate) simulation and grid variables
!> @todo Explain which ones are not initialized here
subroutine initialize_vars
use mod_forest
include 'amrvacdef.f'

integer :: igrid, level, ipe, ig1,ig2
logical :: ok
!-----------------------------------------------------------------------------

! set time, time counter
if(.not.treset)t=zero
if(.not.itreset)it=0
dt=zero
dtimpl=zero
itmin=0

if(.not.time_accurate.or.residmin>smalldouble) then
  residual=one
endif 

! set all dt to zero
dt_grid(1:ngridshi)=zero

! check resolution
if (mod(ixGhi1,2)/=0.or.mod(ixGhi2,2)/=0) then
   call mpistop("mesh widths must give even number grid points")
end if
ixMlo1=ixGlo1+dixB;ixMlo2=ixGlo2+dixB;ixMhi1=ixGhi1-dixB;ixMhi2=ixGhi2-dixB;
if (errorestimate==1) then
   if (mod(ixMhi1-ixMlo1+1,4)/=0.or.mod(ixMhi2-ixMlo2+1,4)/=0) then
      call mpistop("mesh widths must be divisible by 4 for Richardson")
   end if
end if

if (nbufferx1>(ixMhi1-ixMlo1+1).or.nbufferx2>(ixMhi2-ixMlo2+1)) then
   write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
   write(unitterm,*) "Decrease nbufferx or increase mesh size"
   call mpistop("")
end if

! initialize dx arrays on finer (>1) levels
do level=2,mxnest
   dx(1,level) = dx(1,level-1) * half
   dx(2,level) = dx(2,level-1) * half  ! refine ratio 2
end do

! domain decomposition
! physical extent of a grid block at level 1, per dimension
dg1(1)=dx(1,1)*dble(ixGhi1-2*dixB)
dg2(1)=dx(2,1)*dble(ixGhi2-2*dixB)
! number of grid blocks at level 1 in simulation domain, per dimension
ng1(1)=nint((xprobmax1-xprobmin1)/dg1(1))
ng2(1)=nint((xprobmax2-xprobmin2)/dg2(1))
! total number of grid blocks at level 1
nglev1=ng1(1)*ng2(1)

do level=2,mxnest
   dg1(level)=half*dg1(level-1);dg2(level)=half*dg2(level-1);
   ng1(level)=ng1(level-1)*2;ng2(level)=ng2(level-1)*2;
end do

! check that specified stepsize correctly divides domain
ok=((abs(dble(ng1(1))*dg1(1)-(xprobmax1-xprobmin1))<=smalldouble)&
   .and.(abs(dble(ng2(1))*dg2(1)-(xprobmax2-xprobmin2))<=smalldouble))
if (.not.ok) then
   write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
   call mpistop("domain cannot be divided by meshes of given gridsize")
end if


poleB=.false.
if (.not.slab) call set_pole

do igrid=1,ngridshi
   nullify(pwold(igrid)%w,pw(igrid)%w,pw1(igrid)%w, pwCoarse(igrid)%w,&
      pwCoCo(igrid)%w)
   nullify(px(igrid)%x,pxCoarse(igrid)%x)
   nullify(pgeo(igrid)%surfaceC1,pgeo(igrid)%surfaceC2,pgeo(igrid)%surface1,&
      pgeo(igrid)%surface2, pgeo(igrid)%dvolume,pgeo(igrid)%dx)
   nullify(pgeoCoarse(igrid)%surfaceC1,pgeoCoarse(igrid)%surfaceC2,&
      pgeoCoarse(igrid)%surface1,pgeoCoarse(igrid)%surface2,&
       pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%dx)
   if (B0field) then
    nullify(pB0_cell(igrid)%w,pB0_face1(igrid)%w,pB0_face2(igrid)%w)
   end if
   if (nstep>2) then
     nullify(pw2(igrid)%w)
   end if
   if (nstep>3) then
     nullify(pw3(igrid)%w)
   end if
   if (nstep>4) then
      nullify(pw4(igrid)%w)
   end if
   if (residmin>smalldouble) then
     nullify(pwres(igrid)%w)
  end if
end do

! on each processor, create for later use a default patch array
allocate(patchfalse(ixGlo1:ixGhi1,ixGlo2:ixGhi2))
patchfalse(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=.false.

! initialize connectivity data
igridstail=0

! allocate memory for forest data structures
allocate(level_head(mxnest),level_tail(mxnest))
do level=1,mxnest
   nullify(level_head(level)%node,level_tail(level)%node)
end do

allocate(igrid_to_node(ngridshi,0:npe-1))
do ipe=0,npe-1
   do igrid=1,ngridshi
      nullify(igrid_to_node(igrid,ipe)%node)
   end do
end do

allocate(sfc(1:3,ngridshi*npe))

allocate(igrid_to_sfc(ngridshi))

sfc=0
allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

allocate(nleafs_level(1:nlevelshi))

allocate(coarsen(ngridshi,0:npe-1),refine(ngridshi,0:npe-1))
coarsen=.false.
refine=.false.
if (nbufferx1/=0.or.nbufferx2/=0) then
   allocate(buffer(ngridshi,0:npe-1))
   buffer=.false.
end if
allocate(igrid_inuse(ngridshi,0:npe-1))
igrid_inuse=.false.

allocate(tree_root(1:ng1(1),1:ng2(1)))
do ig2=1,ng2(1)
do ig1=1,ng1(1)
   nullify(tree_root(ig1,ig2)%node)
end do
end do



! default the physical scaling parameters:
UNIT_LENGTH   = ONE
UNIT_DENSITY  = ONE
UNIT_VELOCITY = ONE

end subroutine initialize_vars
!=============================================================================
