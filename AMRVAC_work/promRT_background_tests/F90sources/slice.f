!=============================================================================
subroutine write_slice
include 'amrvacdef.f'
! Writes a D-1 slice .dat-file with proper Morton order 
! by Oliver Porth
! 22.Nov 2011
integer :: islice
logical, save :: firstslice=.true.
!-----------------------------------------------------------------------------
if (firstslice) then
   slice=slicenext
   firstslice=.false.
end if

do islice=1,nslices
   call put_slice(slicedir(islice),slicecoord(islice))
end do

slice=slice+1
end subroutine write_slice
!=============================================================================
subroutine put_slice(dir,xslice)
use mod_forest, only: Morton_sub_start, Morton_sub_stop
include 'amrvacdef.f'
! Writes a D-1 slice .dat-file with proper Morton order 
! Can also write csv files by setting sliceascii=.true.
! For ONED simulations, the output will be appended to one csv-file per slice
! In the latter two cases, the slices are sensitive to the saveprim switch
! Thus csv-files with primitive variables are obtained.  
! by Oliver Porth
! 22.Nov 2011
integer, intent(in) :: dir
double precision, intent(in) :: xslice
! .. local ..
integer :: Njgrid, jgrid
integer, dimension(ndim-1) :: ixsubGlo, ixsubGhi
integer, dimension(ndim-1) :: ixsubMlo, ixsubMhi
integer :: size_subblock_io, nx1,nx2, slice_fh, nwexpand
integer, dimension(ndim) :: sizes, subsizes, start
double precision,dimension(0:nw+nwauxio)          :: normconv 

 integer(kind=MPI_ADDRESS_KIND):: size_double, lb
!-----------------------------------------------------------------------------
! Preamble: 
call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,size_double,ierrmpi)
nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;
slice_fh=unitslice

if (sliceascii.or.ndim==1) then
   nwexpand = nwauxio
else
   nwexpand = 0
end if

! Do a last consistency check:
select case(dir)
case(1)
   if(xslice<xprobmin1.or.xslice>xprobmax1) &
   call mpistop("slice out of bounds")

case(2)
   if(xslice<xprobmin2.or.xslice>xprobmax2) &
   call mpistop("slice out of bounds")

end select

! Traverse the forest and fill nodes:
call select_slice(dir,xslice,.false.,slice_fh,normconv)

! Create the MPI-datatype and select indices:


select case(dir)
case (1)
   ixsubGlo(1) = ixGlo2; ixsubGhi(1) = ixGhi2;
   sizes(1) = ixGhi2
   subsizes(1)=nx2
   start(1)=ixMlo2-1
   size_subblock_io=nx2*(nw+nwexpand)*size_double
case (2)
   ixsubGlo(1) = ixGlo1; ixsubGhi(1) = ixGhi1;
   sizes(1) = ixGhi1
   subsizes(1)=nx1
   start(1)=ixMlo1-1
   size_subblock_io=nx1*(nw+nwexpand)*size_double
case default
   call mpistop("slice direction not clear in put_slice")
end select




ixsubMlo(2-1) = ixsubGlo(2-1)+dixB;
ixsubMhi(2-1) = ixsubGhi(2-1)-dixB;


sizes(ndim)=nw+nwexpand
subsizes(ndim)=nw+nwexpand
start(ndim)=0

call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_subblock_io,ierrmpi)
call MPI_TYPE_COMMIT(type_subblock_io,ierrmpi)

sizes(ndim)=2
subsizes(ndim)=2
start(ndim)=0
call MPI_TYPE_CREATE_SUBARRAY(ndim,sizes,subsizes,start, MPI_ORDER_FORTRAN,&
   MPI_DOUBLE_PRECISION, type_subblock_x_io,ierrmpi)
call MPI_TYPE_COMMIT(type_subblock_x_io,ierrmpi)

! local number of sub-grids:
Njgrid=Morton_sub_stop(mype)-Morton_sub_start(mype)+1

! Now output using various schemes: 
if (ndim==1) then 
   call put_slice_zerod
else if (sliceascii) then 
   call put_slice_csv
else
   call put_slice_dat
end if

! If we need the subnodes later, remove deallocation here:
do jgrid=1,Njgrid
   call dealloc_subnode(jgrid)
end do

call MPI_TYPE_FREE(type_subblock_io,ierrmpi)
call MPI_TYPE_FREE(type_subblock_x_io,ierrmpi)


contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine put_slice_csv

character(len=1024) :: filename, xlabel
character(len=79)   :: xxlabel
logical             :: fileopen
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
integer :: iw, ipe, itag
character(len=1024) ::line
integer :: status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------
if (mype==0) then
 inquire(slice_fh,opened=fileopen)
 if(.not.fileopen)then
      ! generate filename: 
      write(xlabel,"(D9.2)")xslice
      xxlabel=trim(xlabel)
      if(xslice>=zero)then
         write(xxlabel(1:1),"(a)") "p"
      else
         write(xxlabel(1:1),"(a)") "n"
      endif
      write(filename,"(a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_x'//trim(xxlabel)//'_n',slice,'.csv'
      open(slice_fh,file=filename,status='unknown',form='formatted')
   end if
   ! get and write the header: 
   call getheadernames(wnamei,xandwnamei,outfilehead)
   line=''
   do iw=1,ndim+nw+nwauxio-1
      line = trim(line)//trim(xandwnamei(iw))//', '
   end do
   line = trim(line)//trim(xandwnamei(ndim+nw+nwauxio))
   write(slice_fh,'(a)')trim(line)
   ! create a recv buffer using allocate, will be deallocated at the end of the routine:
   call alloc_subnode(Njgrid+1,dir,nwauxio)

   ! write to file:
   do jgrid=1, Njgrid
      call put_slice_line(jgrid,slice_fh)
   end do
end if

! Also communicate the normconv array since processor zero might not have it yet:
if (npe>1) then
   do ipe=1,npe-1
      do jgrid=1,Morton_sub_stop(ipe)-Morton_sub_start(ipe)+1
         itag=Morton_sub_start(ipe)+jgrid-1
         if (ipe == mype ) then 
            call MPI_SEND(px_sub(jgrid)%x,1,type_subblock_x_io,0,itag,icomm,&
               ierrmpi)
            call MPI_SEND(pw_sub(jgrid)%w,1,type_subblock_io,0,itag,icomm,&
               ierrmpi)
            call MPI_SEND(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,0,itag,&
               icomm,ierrmpi)
         end if
         if (mype == 0) then
            call MPI_RECV(px_sub(Njgrid+1)%x,1,type_subblock_x_io,ipe,itag,&
               icomm,status,ierrmpi)
            call MPI_RECV(pw_sub(Njgrid+1)%w,1,type_subblock_io,ipe,itag,&
               icomm,status,ierrmpi)
            call MPI_RECV(normconv,nw+nwauxio+1,MPI_DOUBLE_PRECISION,ipe,itag,&
               icomm,status,ierrmpi)
            call put_slice_line(Njgrid+1,slice_fh)
         end if
      end do
end do
endif

if (mype==0) then
   close(slice_fh)
   call dealloc_subnode(Njgrid+1)
end if

end subroutine put_slice_csv
!=============================================================================
subroutine put_slice_line(jout,file_handle)
integer, intent(in) :: jout, file_handle
! .. local ..
character(len=1024) ::line, data
integer :: ix1,ix2,idir,iw
double precision, parameter :: minvalue = 1.0d-99
double precision            :: roundoff
!-----------------------------------------------------------------------------


do ix1=ixsubMlo(1),ixsubMhi(1)

      ! Format the line:
      line = ''
      do idir=1,ndim


         write(data,"(es14.6)")roundoff(px_sub(jout)%x(ix1,idir),minvalue)



         line = trim(line)//trim(data)//', '
      end do
      do iw = 1,nw+nwauxio-1


         write(data,"(es14.6)")roundoff(pw_sub(jout)%w(ix1,iw)&
            *normconv(iw),minvalue)


         line = trim(line)//trim(data)//', '
      end do


      write(data,"(es14.6)")roundoff(pw_sub(jout)%w(ix1,nw+nwauxio)&
         *normconv(nw+nwauxio),minvalue)

      line = trim(line)//trim(data)
      write(file_handle,'(a)')trim(line)

end do



end subroutine put_slice_line
!=============================================================================
subroutine put_slice_dat

integer, dimension(ngridshi) :: iorequest
integer, dimension(MPI_STATUS_SIZE,ngridshi) :: iostatus
integer(kind=MPI_OFFSET_KIND) :: offset
integer :: nsubleafs
character(len=1024) :: filename, xlabel
character(len=79)   :: xxlabel
integer :: amode, status(MPI_STATUS_SIZE), iwrite
!-----------------------------------------------------------------------------

nsubleafs=Morton_sub_stop(npe-1)
! generate filename
write(xlabel,"(D9.2)")xslice
xxlabel=trim(xlabel)
if(xslice>=zero)then
   write(xxlabel(1:1),"(a)") "p"
else
   write(xxlabel(1:1),"(a)") "n"
endif
write(filename,"(a,i1.1,a,i4.4,a)") TRIM(filenameout)//'_d',dir,'_x'//trim(xxlabel)//'_n',slice,'.dat'

if(mype==0) then
   open(unit=slice_fh,file=filename,status='replace')
   close(unit=slice_fh)
end if

amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
call MPI_FILE_OPEN(icomm,filename,amode,MPI_INFO_NULL,slice_fh,ierrmpi)
iorequest=MPI_REQUEST_NULL
iwrite=0

do jgrid=1,Njgrid
   iwrite=iwrite+1
   offset=int(size_subblock_io,kind=MPI_OFFSET_KIND) *int(Morton_sub_start&
      (mype)+jgrid-2,kind=MPI_OFFSET_KIND)
   call MPI_FILE_IWRITE_AT(slice_fh,offset,pw_sub(jgrid)%w,1,type_subblock_io,&
       iorequest(iwrite),ierrmpi)
end do

if (iwrite>0) call MPI_WAITALL(iwrite,iorequest,iostatus,ierrmpi)
call MPI_BARRIER(icomm, ierrmpi)
call MPI_FILE_CLOSE(slice_fh,ierrmpi)

if (mype==0) then
   amode=ior(MPI_MODE_APPEND,MPI_MODE_WRONLY)
   call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL, slice_fh,&
      ierrmpi)

   call select_slice(dir,xslice,.true.,slice_fh,normconv)

   call MPI_FILE_WRITE(slice_fh,subsizes(2-1),1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,eqpar,neqpar+nspecialpar, &
      MPI_DOUBLE_PRECISION,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,nsubleafs,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,levmax_sub,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,ndim-1,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,ndir,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,nw,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,neqpar+nspecialpar,1,MPI_INTEGER,status,&
      ierrmpi)
   call MPI_FILE_WRITE(slice_fh,it,1,MPI_INTEGER,status,ierrmpi)
   call MPI_FILE_WRITE(slice_fh,t,1,MPI_DOUBLE_PRECISION,status,ierrmpi)

   call MPI_FILE_CLOSE(slice_fh,ierrmpi)
end if

end subroutine put_slice_dat
!=============================================================================
subroutine put_slice_zerod

integer::  iw
character(len=10) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead
logical, save :: opened=.false.
character(len=1024) ::line, data
character(len=1024) :: filename, xlabel
character(len=79)   :: xxlabel
integer :: amode, iwrite, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------

end subroutine put_slice_zerod

end subroutine put_slice
!=============================================================================
double precision function roundoff(val,minval)
implicit none
double precision,intent(in)         :: val, minval
if (abs(val) .gt. minval) then 
   roundoff = val
else 
   roundoff = 0.0d0
end if
end function roundoff
!=============================================================================
subroutine select_slice(dir,xslice,writeonly,file_handle,normconv)
use mod_forest, only: tree_node_ptr, tree_root, Morton_sub_start,&
    Morton_sub_stop
include 'amrvacdef.f'
integer, intent(in) :: dir
double precision, intent(in) :: xslice
integer, intent(in) :: file_handle
logical, intent(in) :: writeonly
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
! .. local ..
integer :: ig1,ig2, jgrid, slice_fh, ipe, mylevmax
integer, dimension(nlevelshi) :: igslice
!-----------------------------------------------------------------------------
jgrid = 0
mylevmax = 0

! Find the global slice index for every level:
call get_igslice(dir,xslice,igslice)

! Traverse forest to find grids indicating slice:


select case(dir)
case (1)
   ig1 = igslice(1)
   do ig2=1,ng2(1)
         call traverse_slice(tree_root(ig1,ig2))
   end do
case (2)
   ig2 = igslice(1)
      do ig1=1,ng1(1)
         call traverse_slice(tree_root(ig1,ig2))
      end do
case default
   call mpistop("slice direction not clear in select_slice")
end select



if (.not.writeonly) then
! Synchronize the levmax_sub for output (only rank 0 needs it): 
   levmax_sub = mylevmax
   call MPI_ALLREDUCE(MPI_IN_PLACE,levmax_sub,1,MPI_INTEGER,MPI_MAX,icomm,&
      ierrmpi)

! Communicate the subgrid indices according to new Morton sub-sfc:
   Morton_sub_start(:) = 1
   do ipe=0,npe-1
      call MPI_GATHER(jgrid,1,MPI_INTEGER,Morton_sub_stop,1,MPI_INTEGER,ipe,&
         icomm,ierrmpi)
   end do

   do ipe = 0, npe-2
      Morton_sub_start(ipe+1) = Morton_sub_stop(ipe)+Morton_sub_start(ipe+1)
      Morton_sub_stop(ipe+1)  = Morton_sub_stop(ipe)+Morton_sub_stop(ipe+1)
   end do
end if

contains
!=============================================================================
! internal procedures
!=============================================================================
recursive subroutine traverse_slice(tree)
implicit none
type(tree_node_ptr) :: tree
integer :: ic1,ic2
integer, dimension(MPI_STATUS_SIZE) :: status
!-----------------------------------------------------------------------------
if (writeonly) then
call MPI_FILE_WRITE(file_handle,tree%node%leaf,1,MPI_LOGICAL, status,ierrmpi)
end if

if (tree%node%leaf) then
   if (tree%node%ipe == mype.and..not.writeonly) then
      mylevmax = max(mylevmax,tree%node%level)
      call fill_subnode(tree%node%igrid,tree%node%active,jgrid,dir,xslice,&
         normconv)
   end if
   return
end if
! We are out for leaves now, continue for branches

! Get the correct child:
select case (dir)
case (1)
   ic1 = igslice(tree%node%level+1) - 2 * tree%node%ig1 + 2

case (2)
   ic2 = igslice(tree%node%level+1) - 2 * tree%node%ig2 + 2

case default
   call mpistop("slice direction not clear in traverse_slice")
end select

! Recursively descend into the correct child branch:


select case(dir)
case (1)
   do ic2=1,2
      call traverse_slice(tree%node%child(ic1,ic2))
   end do
case (2)
   do ic1=1,2
      call traverse_slice(tree%node%child(ic1,ic2))
   end do
case default
   call mpistop("slice direction not clear in traverse_slice")
end select



end subroutine traverse_slice
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine select_slice
!=============================================================================
subroutine fill_subnode(igrid,active,jgrid,dir,xslice,normconv)
include 'amrvacdef.f'
integer, intent(in) :: igrid, dir
integer, intent(inout) :: jgrid
logical, intent(in)  :: active
double precision, intent(in) :: xslice
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 
! .. local ..
integer :: ixslice, i, nwexpand
logical                :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2)&
   =.false.,mask(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=.false.
double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw+nwauxio)   :: w
!-----------------------------------------------------------------------------
! increase grid-count:
jgrid=jgrid+1
! Allocate subdim solution array:
if (sliceascii.or.ndim==1) then
   nwexpand = nwauxio
else
   nwexpand = 0
end if
call alloc_subnode(jgrid,dir,nwexpand)
call fill_subnode_info(igrid,jgrid,dir)

mask(ixMlo1:ixMhi1,ixMlo2:ixMhi2)=.true.

! Now hunt for the index closest to the slice:


select case (dir)
case (1)
   ixslice = minloc(dabs(xslice-px(igrid)%x(:,ixMlo2,dir)),1,mask(:,ixMlo2))
case (2)
   ixslice = minloc(dabs(xslice-px(igrid)%x(ixMlo1,:,dir)),1,mask(ixMlo1,:))
case default
   call mpistop("slice direction not clear in fill_subnode")
end select



! Make a local copy of the pw(igrid)%w() array and compute nwauxio variables on this: 
w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=pw(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   1:nw)
if(saveprim.and.(sliceascii.or.ndim==1)) then 
   call primitive(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,w,&
      px(igrid)%x)
   normconv(0:nw)=normvar(0:nw)
else
  normconv(0:nw)=one
end if

if(nwexpand>0)then
! auxiliary io variables can be computed and added by user
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))
dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     myB0      => pB0_cell(igrid)
     myB0_face1 => pB0_face1(igrid)
     myB0_face2 => pB0_face2(igrid)
  end if
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one
  call specialvar_output(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
     ixMhi2,w,px(igrid)%x,normconv)
endif


if (B0field) w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b1_:b0_+ndir) = w(ixMlo1:ixMhi1,&
   ixMlo2:ixMhi2,b1_:b0_+ndir) + pB0_cell(igrid)%w(ixMlo1:ixMhi1,&
   ixMlo2:ixMhi2,1:ndir)

if((.not.saveprim) .and. B0field) then
   w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,e_)=w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,e_) &
           +half*( pB0_cell(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b1_-b0_)**2&
              +pB0_cell(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b2_-b0_)**2&
              +pB0_cell(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b3_-b0_)**2 ) &
           + ( w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b1_)*pB0_cell(igrid)%w&
              (ixMlo1:ixMhi1,ixMlo2:ixMhi2,b1_-b0_)+w(ixMlo1:ixMhi1,&
              ixMlo2:ixMhi2,b2_)*pB0_cell(igrid)%w(ixMlo1:ixMhi1,&
              ixMlo2:ixMhi2,b2_-b0_)+w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b3_)&
              *pB0_cell(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,b3_-b0_) )
endif



! Fill the subdimensional solution and position array:



select case (dir)
case (1)
   pw_sub(jgrid)%w(ixGlo2:ixGhi2,1:nw+nwexpand) = w(ixslice,ixGlo2:ixGhi2,&
      1:nw+nwexpand)
   px_sub(jgrid)%x(ixGlo2:ixGhi2,1:ndim) = px(igrid)%x(ixslice,ixGlo2:ixGhi2,&
      1:ndim)
case (2)
   pw_sub(jgrid)%w(ixGlo1:ixGhi1,1:nw+nwexpand) = w(ixGlo1:ixGhi1,ixslice,&
      1:nw+nwexpand)
   px_sub(jgrid)%x(ixGlo1:ixGhi1,1:ndim) = px(igrid)%x(ixGlo1:ixGhi1,ixslice,&
      1:ndim)
case default
   call mpistop("slice direction not clear in fill_subnode")
end select



end subroutine fill_subnode
!=============================================================================
subroutine alloc_subnode(jgrid,dir,nwexpand)
include 'amrvacdef.f'
integer, intent(in) :: jgrid, dir, nwexpand
!-----------------------------------------------------------------------------
! take care, what comes out is not necessarily a right handed system!


select case (dir)
case (1)
allocate(pw_sub(jgrid)%w(ixGlo2:ixGhi2,1:nw+nwexpand),&
     px_sub(jgrid)%x(ixGlo2:ixGhi2,1:ndim))
case (2)
allocate(pw_sub(jgrid)%w(ixGlo1:ixGhi1,1:nw+nwexpand),&
     px_sub(jgrid)%x(ixGlo1:ixGhi1,1:ndim))
case default
   call mpistop("slice direction not clear in alloc_subnode")
end select


end subroutine alloc_subnode
!=============================================================================
subroutine dealloc_subnode(jgrid)
include 'amrvacdef.f'
integer, intent(in) :: jgrid
!-----------------------------------------------------------------------------
if (jgrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_subnode")
end if

deallocate(pw_sub(jgrid)%w,px_sub(jgrid)%x)

! reset the global node info:
node_sub(:,jgrid)=0
rnode_sub(:,jgrid)=zero

end subroutine dealloc_subnode
!=============================================================================
subroutine fill_subnode_info(igrid,jgrid,dir)
include 'amrvacdef.f'
integer, intent(in) :: igrid,jgrid,dir
!-----------------------------------------------------------------------------

node_sub(plevel_,jgrid)=node(plevel_,igrid)


select case(dir)
case (1)
   node_sub(pig1_,jgrid)=node(pig2_,igrid)
   rnode_sub(rpdx1_,jgrid)=rnode(rpdx2_,igrid)
   rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin2_,igrid)
case (2)
   node_sub(pig1_,jgrid)=node(pig1_,igrid)
   rnode_sub(rpdx1_,jgrid)=rnode(rpdx1_,igrid)
   rnode_sub(rpxmin1_,jgrid)=rnode(rpxmin1_,igrid)
case default
   call mpistop("slice direction not clear in fill_subnode_info")
end select



end subroutine fill_subnode_info
!=============================================================================

subroutine get_igslice(dir,x,igslice)
include 'amrvacdef.f'
integer, intent(in) :: dir
double precision, intent(in) :: x
integer, dimension(nlevelshi), intent(out) :: igslice
! .. local ..
integer level
!-----------------------------------------------------------------------------
select case (dir)
case (1)
do level = 1, mxnest
   igslice(level) = int((x-xprobmin1)/dg1(level))+1
! Gets out of domain when x==xprobmax1, not caught by put_slice, so limit:
   if (x>=xprobmax1) igslice(level) =  int((xprobmax1-xprobmin1)/dg1(level))
! This is already caught by control in put_slice, but anyways:
   if (x<=xprobmin1) igslice(level) =  1
end do
case (2)
do level = 1, mxnest
   igslice(level) = int((x-xprobmin2)/dg2(level))+1
! Gets out of domain when x==xprobmax2, not caught by put_slice, so limit:
   if (x>=xprobmax2) igslice(level) =  int((xprobmax2-xprobmin2)/dg2(level))
! This is already caught by control in put_slice, but anyways:
   if (x<=xprobmin2) igslice(level) =  1
end do
case default
   call mpistop("slice direction not clear in get_igslice")
end select
end subroutine get_igslice
!=============================================================================
