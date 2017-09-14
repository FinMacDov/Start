!=============================================================================
subroutine getbc(time,qdt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,pwuse,pwuseCo,&
   pgeoFi,pgeoCo,richardson,nwstart,nwbc)

include 'amrvacdef.f'

double precision, intent(in)               :: time, qdt
integer, intent(in)                        :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
   nwstart,nwbc
type(walloc), dimension(ngridshi)          :: pwuse, pwuseCo
type(geoalloc), target,dimension(ngridshi) :: pgeoFi, pgeoCo
logical, intent(in)                        :: richardson

integer :: ixMmin1,ixMmin2,ixMmax1,ixMmax2, ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2, ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2, idims, iside
integer :: my_neighbor_type, ipole
integer :: iigrid, igrid, ineighbor, ipe_neighbor
integer :: nrecvs, nsends, isizes
integer :: ixRmin1,ixRmin2,ixRmax1,ixRmax2, ixSmin1,ixSmin2,ixSmax1,ixSmax2
integer :: ixBmin1,ixBmin2,ixBmax1,ixBmax2
integer :: kmin1,kmin2,kmax1,kmax2
integer :: i1,i2, n_i1,n_i2, ic1,ic2, inc1,inc2, n_inc1,n_inc2
integer, dimension(-1:1) :: ixS_srl_min1,ixS_srl_min2,ixS_srl_max1,&
   ixS_srl_max2, ixR_srl_min1,ixR_srl_min2,ixR_srl_max1,ixR_srl_max2,&
    ixS_r_min1,ixS_r_min2,ixS_r_max1,ixS_r_max2
integer, dimension(0:3) :: ixR_r_min1,ixR_r_min2,ixR_r_max1,ixR_r_max2,&
    ixS_p_min1,ixS_p_min2,ixS_p_max1,ixS_p_max2, ixR_p_min1,ixR_p_min2,&
   ixR_p_max1,ixR_p_max2, ixS_old_min1,ixS_old_min2,ixS_old_max1,ixS_old_max2,&
    ixR_old_min1,ixR_old_min2,ixR_old_max1,ixR_old_max2
integer, dimension(-1:1,-1:1) :: type_send_srl, type_recv_srl, type_send_r
integer, dimension(0:3,0:3) :: type_recv_r, type_send_p, type_recv_p,&
    type_send_old, type_recv_old
integer :: isend_buf(npwbuf), ipwbuf
type(walloc) :: pwbuf(npwbuf)
logical  :: isphysbound

double precision :: time_bcin

!-----------------------------------------------------------------------------
time_bcin=MPI_WTIME()

call init_bc
if (internalboundary) then 
   call getintbc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,pwuse)
end if

! default : no singular axis
ipole=0

irecv=0
isend=0
isend_buf=0
ipwbuf=1
nrecvs=nrecv_bc_srl+nrecv_bc_r
nsends=nsend_bc_srl+nsend_bc_r
if (richardson) then
   nrecvs=nrecvs+nrecv_bc_p
   nsends=nsends+nsend_bc_p
end if
if (nrecvs>0) then
   allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
   recvrequest=MPI_REQUEST_NULL
end if
if (nsends>0) then
   allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
   sendrequest=MPI_REQUEST_NULL
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
   do i2=-1,1
   do i1=-1,1
      if (i1==0.and.i2==0) cycle

      my_neighbor_type=neighbor_type(i1,i2,igrid)
      select case (my_neighbor_type)
      case (2)
         if (richardson) call bc_recv_old
      case (3)
         call bc_recv_srl
      case (4)
         call bc_recv_restrict
      end select
   end do
   end do
end do

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
   if (any(neighbor_type(:,:,igrid)==2)) then
      if (richardson) then
         dxlevel(1)=two*rnode(rpdx1_,igrid)
         dxlevel(2)=two*rnode(rpdx2_,igrid);
      else
         dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
      end if
      call coarsen_grid(pwuse(igrid)%w,px(igrid)%x,ixGmin1,ixGmin2,ixGmax1,&
         ixGmax2,ixMmin1,ixMmin2,ixMmax1,ixMmax2,pwuseCo(igrid)%w,&
         pxCoarse(igrid)%x, ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixCoMmin1,&
         ixCoMmin2,ixCoMmax1,ixCoMmax2,pgeoFi(igrid),pgeoCo(igrid),&
          coarsenprimitive,.true.)
   end if
   do i2=-1,1
   do i1=-1,1
      if (i1==0.and.i2==0) cycle

      
      my_neighbor_type=neighbor_type(i1,i2,igrid)
      select case (my_neighbor_type)
      case (2)
         call bc_send_restrict
      case (3)
         call bc_send_srl
      case (4)
         if (richardson) call bc_send_old
      end select
   end do
   end do
end do

if (irecv/=nrecvs) then
   call mpistop("number of recvs in phase1 in amr_ghostcells is incorrect")
end if
if (isend/=nsends) then
   call mpistop("number of sends in phase1 in amr_ghostcells is incorrect")
end if

if (irecv>0) then
   call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
   deallocate(recvstatus,recvrequest)
end if
if (isend>0) then
   call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendstatus,sendrequest)
   do ipwbuf=1,npwbuf
      if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
   end do
end if


if (.not.richardson) then
   irecv=0
   isend=0
   isend_buf=0
   ipwbuf=1
   nrecvs=nrecv_bc_p
   nsends=nsend_bc_p
   if (nrecvs>0) then
      allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
      recvrequest=MPI_REQUEST_NULL
   end if
   if (nsends>0) then
      allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
      sendrequest=MPI_REQUEST_NULL
   end if

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
      do i2=-1,1
      do i1=-1,1
         if (i1==0.and.i2==0) cycle

         my_neighbor_type=neighbor_type(i1,i2,igrid)
         if (my_neighbor_type==2) call bc_recv_prolong
      end do
      end do
   end do
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
     if (any(neighbor_type(:,:,igrid)==4)) then
      do i2=-1,1
      do i1=-1,1
         if (i1==0.and.i2==0) cycle

         
         my_neighbor_type=neighbor_type(i1,i2,igrid)
         if (my_neighbor_type==4) call bc_send_prolong
      end do
      end do
     end if
   end do


   if (irecv/=nrecvs) then
      call mpistop("number of recvs in phase2 in amr_ghostcells is incorrect")
   end if
   if (isend/=nsends) then
      call mpistop("number of sends in phase2 in amr_ghostcells is incorrect")
   end if

   if (irecv>0) then
      call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
      deallocate(recvstatus,recvrequest)
   end if

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
      if (any(neighbor_type(:,:,igrid)==2)) then
         do i2=-1,1
         do i1=-1,1
            if (i1==0.and.i2==0) cycle
            my_neighbor_type=neighbor_type(i1,i2,igrid)
            if (my_neighbor_type==2) call bc_prolong
         end do
         end do
      end if
   end do

   if (isend>0) then
      call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
      deallocate(sendstatus,sendrequest)
      do ipwbuf=1,npwbuf
         if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
      end do
   end if
end if

if(bcphys) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     saveigrid=igrid
     dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
     do idims=1,ndim
        ! to avoid using as yet unknown corner info in more than 1D, we
        ! fill only interior mesh ranges of the ghost cell ranges at first,
        ! and progressively enlarge the ranges to include corners later
        kmin1=0; kmax1=0;
        
         kmin2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0,-1,igrid)==1)
         kmax2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0, 1,igrid)==1)
        
        ixBmin1=ixGmin1+kmin1*dixB;ixBmin2=ixGmin2+kmin2*dixB;
        ixBmax1=ixGmax1-kmax1*dixB;ixBmax2=ixGmax2-kmax2*dixB;
        do iside=1,2
           i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3);
           if (aperiodB(idims)) then 
              call physbound(i1,i2,igrid,isphysbound)
              if (neighbor_type(i1,i2,igrid)/=1 .and. .not. isphysbound) cycle
           else 
              if (neighbor_type(i1,i2,igrid)/=1) cycle
           end if
           if (richardson) then
              if(.not.slab)mygeo=>pgeoCo(igrid)
              call bc_phys(iside,idims,time,qdt,pwuse(igrid)%w,&
                 pxCoarse(igrid)%x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
                 ixBmin2,ixBmax1,ixBmax2)
           else
              if(.not.slab)mygeo=>pgeoFi(igrid)
              if (B0field) then
                 myB0_cell => pB0_cell(igrid)
                 myB0_face1 => pB0_face1(igrid)
                 myB0_face2 => pB0_face2(igrid)
              end if
              call bc_phys(iside,idims,time,qdt,pwuse(igrid)%w,px(igrid)%x,&
                 ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,ixBmax1,&
                 ixBmax2)
           end if
        end do
     end do
  end do
end if

if (npe>1) call put_bc_comm_types

if (nwaux>0) call fix_auxiliary

time_bc=time_bc+(MPI_WTIME()-time_bcin)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,i2,igrid)
ipe_neighbor=neighbor(2,i1,i2,igrid)


if (ipole==0) then
   n_i1=-i1;n_i2=-i2;
   if (ipe_neighbor==mype) then
      ixSmin1=ixS_srl_min1(i1);ixSmin2=ixS_srl_min2(i2)
      ixSmax1=ixS_srl_max1(i1);ixSmax2=ixS_srl_max2(i2);
      ixRmin1=ixR_srl_min1(n_i1);ixRmin2=ixR_srl_min2(n_i2)
      ixRmax1=ixR_srl_max1(n_i1);ixRmax2=ixR_srl_max2(n_i2);
      pwuse(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,nwstart+1:nwstart&
         +nwbc)=pwuse(igrid)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
         +1:nwstart+nwbc)
   else
      isend=isend+1
      itag=(3**2+4**2)*(ineighbor-1)+(n_i1+1)*3**(1-1)+(n_i2+1)*3**(2-1)
      call MPI_ISEND(pwuse(igrid)%w,1,type_send_srl(i1,i2), ipe_neighbor,itag,&
         icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixSmin1=ixS_srl_min1(i1);ixSmin2=ixS_srl_min2(i2);ixSmax1=ixS_srl_max1(i1)
   ixSmax2=ixS_srl_max2(i2);
   select case (ipole)
   case (1)
      n_i1=i1;n_i2=-i2;
   case (2)
      n_i1=-i1;n_i2=i2;
   end select
   if (ipe_neighbor==mype) then
      ixRmin1=ixR_srl_min1(n_i1);ixRmin2=ixR_srl_min2(n_i2)
      ixRmax1=ixR_srl_max1(n_i1);ixRmax2=ixR_srl_max2(n_i2);
      call pole_copy(pwuse(ineighbor),ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
         pwuse(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
            isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
         +1:nwstart+nwbc))
      call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmin2,ixSmax1,ixSmax2,&
         pwuse(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**2+4**2)*(ineighbor-1)+(n_i1+1)*3**(1-1)+(n_i2+1)*3**(2-1)
      isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*nwbc
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
          ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_srl
!=============================================================================
subroutine bc_send_restrict
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)) return

ineighbor=neighbor(1,i1,i2,igrid)
ipe_neighbor=neighbor(2,i1,i2,igrid)

if (ipole==0) then
   n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;
   if (ipe_neighbor==mype) then
      ixSmin1=ixS_r_min1(i1);ixSmin2=ixS_r_min2(i2);ixSmax1=ixS_r_max1(i1)
      ixSmax2=ixS_r_max2(i2);
      ixRmin1=ixR_r_min1(n_inc1);ixRmin2=ixR_r_min2(n_inc2)
      ixRmax1=ixR_r_max1(n_inc1);ixRmax2=ixR_r_max2(n_inc2);
      pwuse(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,nwstart+1:nwstart&
         +nwbc)=pwuseCo(igrid)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
         +1:nwstart+nwbc)
   else
      isend=isend+1
      itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
      call MPI_ISEND(pwuseCo(igrid)%w,1,type_send_r(i1,i2), ipe_neighbor,itag,&
         icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixSmin1=ixS_r_min1(i1);ixSmin2=ixS_r_min2(i2);ixSmax1=ixS_r_max1(i1)
   ixSmax2=ixS_r_max2(i2);
   select case (ipole)
   case (1)
      n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;
   case (2)
      n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);
   end select
   if (ipe_neighbor==mype) then
      ixRmin1=ixR_r_min1(n_inc1);ixRmin2=ixR_r_min2(n_inc2)
      ixRmax1=ixR_r_max1(n_inc1);ixRmax2=ixR_r_max2(n_inc2);
      call pole_copy(pwuse(ineighbor),ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
         pwuseCo(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
            isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
         +1:nwstart+nwbc))
      call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmin2,ixSmax1,ixSmax2,&
         pwuseCo(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
      isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*nwbc
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
          ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_restrict
!=============================================================================
subroutine bc_send_prolong
integer :: ii1,ii2
!-----------------------------------------------------------------------------
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
   inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1

   ixSmin1=ixS_p_min1(inc1);ixSmin2=ixS_p_min2(inc2);ixSmax1=ixS_p_max1(inc1)
   ixSmax2=ixS_p_max2(inc2);

   if(bcphys) then
     do idims=1,ndim
        do iside=1,2
           ii1=kr(1,idims)*(2*iside-3);ii2=kr(2,idims)*(2*iside-3);

           if (neighbor_type(ii1,ii2,igrid)/=1) cycle

           if ((  (iside==1.and.idims==1.and.ixSmin1<ixMlo1)&
              .or. (iside==1.and.idims==2.and.ixSmin2<ixMlo2)) &
              .or.( (iside==2.and.idims==1.and.ixSmax1>ixMhi1)&
              .or. (iside==2.and.idims==2.and.ixSmax2>ixMhi2)))then
            ixBmin1=merge(ixGmin1,ixSmin1,idims==1)
            ixBmin2=merge(ixGmin2,ixSmin2,idims==2);
            ixBmax1=merge(ixGmax1,ixSmax1,idims==1)
            ixBmax2=merge(ixGmax2,ixSmax2,idims==2);
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
            kmin1=0; kmax1=0;
           
            kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,&
               -1,  igrid)==1)
            kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,  igrid)&
               ==1)
            if(neighbor_type(0,-1,igrid)==1.and.(neighbor_type(1,0,igrid)==1&
               .or. neighbor_type(-1,0,igrid)==1) .and. i2== 1) kmin2=0
            if(neighbor_type(0, 1,igrid)==1.and.(neighbor_type(1,0,igrid)==1&
               .or. neighbor_type(-1,0,igrid)==1) .and. i2==-1) kmax2=0
           
            ixBmin1=ixBmin1+kmin1;ixBmin2=ixBmin2+kmin2;
            ixBmax1=ixBmax1-kmax1;ixBmax2=ixBmax2-kmax2;

            if(.not.slab)mygeo=>pgeoFi(igrid)
            if (B0field) then
              myB0_cell => pB0_cell(igrid)
              myB0_face1 => pB0_face1(igrid)
              myB0_face2 => pB0_face2(igrid)
            end if

            call bc_phys(iside,idims,time,qdt,pwuse(igrid)%w, px(igrid)%x,&
               ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,ixBmax1,&
               ixBmax2)
           end if
        end do
     end do
   end if

   ineighbor=neighbor_child(1,inc1,inc2,igrid)
   ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)

   if (ipole==0) then
      n_i1=-i1;n_i2=-i2;
      n_inc1=ic1+n_i1;n_inc2=ic2+n_i2;
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_p_min1(n_inc1);ixRmin2=ixR_p_min2(n_inc2)
         ixRmax1=ixR_p_max1(n_inc1);ixRmax2=ixR_p_max2(n_inc2);
         pwuseCo(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,nwstart&
            +1:nwstart+nwbc) =pwuse(igrid)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
            nwstart+1:nwstart+nwbc)
      else
         isend=isend+1
         itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
         call MPI_ISEND(pwuse(igrid)%w,1,type_send_p(inc1,inc2), ipe_neighbor,&
            itag,icomm,sendrequest(isend),ierrmpi)
      end if
   else
      select case (ipole)
      case (1)
         n_inc1=inc1;n_inc2=ic2-i2;
      case (2)
         n_inc1=ic1-i1;n_inc2=inc2;
      end select
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_p_min1(n_inc1);ixRmin2=ixR_p_min2(n_inc2)
         ixRmax1=ixR_p_max1(n_inc1);ixRmax2=ixR_p_max2(n_inc2);
         call pole_copy(pwuseCo(ineighbor),ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
            pwuse(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
               isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
            +1:nwstart+nwbc))
         call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmin2,ixSmax1,ixSmax2,&
            pwuse(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
         isend=isend+1
         isend_buf(ipwbuf)=isend
         itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
         isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*nwbc
         call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
         ipwbuf=1+modulo(ipwbuf,npwbuf)
      end if
   end if
end do
end do

end subroutine bc_send_prolong
!=============================================================================
subroutine bc_send_old
!-----------------------------------------------------------------------------
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
   inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1

   ineighbor=neighbor_child(1,inc1,inc2,igrid)
   ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)

   if (ipole==0) then
      n_i1=-i1;n_i2=-i2;
      n_inc1=ic1+n_i1;n_inc2=ic2+n_i2;
      if (ipe_neighbor==mype) then
         ixSmin1=ixS_old_min1(inc1);ixSmin2=ixS_old_min2(inc2)
         ixSmax1=ixS_old_max1(inc1);ixSmax2=ixS_old_max2(inc2);
         ixRmin1=ixR_old_min1(n_inc1);ixRmin2=ixR_old_min2(n_inc2)
         ixRmax1=ixR_old_max1(n_inc1);ixRmax2=ixR_old_max2(n_inc2);
         pwuse(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,nwstart+1:nwstart&
            +nwbc) =pwold(igrid)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
            +1:nwstart+nwbc)
      else
         isend=isend+1
         itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
         call MPI_ISEND(pwold(igrid)%w,1,type_send_old(inc1,inc2),&
             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      end if
   else
      ixSmin1=ixS_old_min1(inc1);ixSmin2=ixS_old_min2(inc2)
      ixSmax1=ixS_old_max1(inc1);ixSmax2=ixS_old_max2(inc2);
      select case (ipole)
      case (1)
         n_inc1=inc1;n_inc2=ic2-i2;
      case (2)
         n_inc1=ic1-i1;n_inc2=inc2;
      end select
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_old_min1(n_inc1);ixRmin2=ixR_old_min2(n_inc2)
         ixRmax1=ixR_old_max1(n_inc1);ixRmax2=ixR_old_max2(n_inc2);
         call pole_copy(pwuse(ineighbor),ixRmin1,ixRmin2,ixRmax1,ixRmax2,&
            pwold(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
               isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,nwstart&
            +1:nwstart+nwbc))
         call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmin2,ixSmax1,ixSmax2,&
            pwold(igrid),ixSmin1,ixSmin2,ixSmax1,ixSmax2)
         isend=isend+1
         isend_buf(ipwbuf)=isend
         itag=(3**2+4**2)*(ineighbor-1)+3**2+n_inc1*4**(1-1)+n_inc2*4**(2-1)
         isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*nwbc
         call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
         ipwbuf=1+modulo(ipwbuf,npwbuf)
      end if
   end if

end do
end do

end subroutine bc_send_old
!=============================================================================
subroutine bc_recv_srl
!-----------------------------------------------------------------------------
ipe_neighbor=neighbor(2,i1,i2,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   itag=(3**2+4**2)*(igrid-1)+(i1+1)*3**(1-1)+(i2+1)*3**(2-1)
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_srl(i1,i2), ipe_neighbor,itag,&
      icomm,recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_srl
!=============================================================================
subroutine bc_recv_restrict
!-----------------------------------------------------------------------------
do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
   inc2=2*i2+ic2
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1
   ipe_neighbor=neighbor_child(2,inc1,inc2,igrid)
   if (ipe_neighbor/=mype) then
      irecv=irecv+1
      itag=(3**2+4**2)*(igrid-1)+3**2+inc1*4**(1-1)+inc2*4**(2-1)
      call MPI_IRECV(pwuse(igrid)%w,1,type_recv_r(inc1,inc2), ipe_neighbor,&
         itag,icomm,recvrequest(irecv),ierrmpi)
   end if
end do
end do

end subroutine bc_recv_restrict
!=============================================================================
subroutine bc_recv_prolong
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)) return

ipe_neighbor=neighbor(2,i1,i2,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc1=ic1+i1;inc2=ic2+i2;
   itag=(3**2+4**2)*(igrid-1)+3**2+inc1*4**(1-1)+inc2*4**(2-1)
   call MPI_IRECV(pwuseCo(igrid)%w,1,type_recv_p(inc1,inc2), ipe_neighbor,&
      itag,icomm,recvrequest(irecv),ierrmpi)  
end if

end subroutine bc_recv_prolong
!=============================================================================
subroutine bc_recv_old
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);ic2=1+modulo(node(pig2_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)) return

ipe_neighbor=neighbor(2,i1,i2,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc1=ic1+i1;inc2=ic2+i2;
   itag=(3**2+4**2)*(igrid-1)+3**2+inc1*4**(1-1)+inc2*4**(2-1)
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_old(inc1,inc2), ipe_neighbor,&
      itag,icomm,recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_old
!=============================================================================
subroutine bc_prolong

integer :: ixFimin1,ixFimin2,ixFimax1,ixFimax2,ixComin1,ixComin2,ixComax1,&
   ixComax2,ii1,ii2
double precision :: dxFi1,dxFi2, dxCo1,dxCo2, xFimin1,xFimin2, xComin1,&
   xComin2, invdxCo1,invdxCo2
!-----------------------------------------------------------------------------
ixFimin1=ixR_srl_min1(i1);ixFimin2=ixR_srl_min2(i2);ixFimax1=ixR_srl_max1(i1)
ixFimax2=ixR_srl_max2(i2);

dxFi1=rnode(rpdx1_,igrid);dxFi2=rnode(rpdx2_,igrid);
dxCo1=two*dxFi1;dxCo2=two*dxFi2;
invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;

xFimin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxFi1
xFimin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxFi2;
xComin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxCo1
xComin2=rnode(rpxmin2_,igrid)-dble(dixB)*dxCo2;


! moved the physical boundary filling here, to only fill the
! part needed

ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+1-1
ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-xComin2)*invdxCo2)+1-1;
ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+1
ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-xComin2)*invdxCo2)+1+1;

if(bcphys) then
  do idims=1,ndim
     do iside=1,2
        ii1=kr(1,idims)*(2*iside-3);ii2=kr(2,idims)*(2*iside-3);
  
        if (neighbor_type(ii1,ii2,igrid)/=1) cycle
  
        if  (( (iside==1.and.idims==1.and.ixComin1<ixCoGmin1&
           +dixB).or.(iside==1.and.idims==2.and.ixComin2<ixCoGmin2&
           +dixB) ) .or.( (iside==2.and.idims==1.and.ixComax1>ixCoGmax1&
           -dixB).or. (iside==2.and.idims==2.and.ixComax2>ixCoGmax2&
           -dixB)))then
          ixBmin1=merge(ixCoGmin1,ixComin1,idims==1)
          ixBmin2=merge(ixCoGmin2,ixComin2,idims==2);
          ixBmax1=merge(ixCoGmax1,ixComax1,idims==1)
          ixBmax2=merge(ixCoGmax2,ixComax2,idims==2);
          if(.not.slab)mygeo=>pgeoCo(igrid)
  
          call bc_phys(iside,idims,time,0.d0,pwuseCo(igrid)%w,&
              pxCoarse(igrid)%x,ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
             ixBmin1,ixBmin2,ixBmax1,ixBmax2)
        end if
     end do
  end do
end if

if (amrentropy) then
   call e_to_rhos(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,ixComin2,&
      ixComax1,ixComax2,pwuseCo(igrid)%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
   call primitive(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,ixComin2,&
      ixComax1,ixComax2,pwuseCo(igrid)%w,pxCoarse(igrid)%x)
end if

select case (typeghostfill)
case ("linear")
   call interpolation_linear(pwuse(igrid),ixFimin1,ixFimin2,ixFimax1,ixFimax2,&
      dxFi1,dxFi2,xFimin1,xFimin2, pwuseCo(igrid),dxCo1,dxCo2,invdxCo1,&
      invdxCo2,xComin1,xComin2)
case ("copy")
   call interpolation_copy(pwuse(igrid),ixFimin1,ixFimin2,ixFimax1,ixFimax2,&
      dxFi1,dxFi2,xFimin1,xFimin2, pwuseCo(igrid),dxCo1,dxCo2,invdxCo1,&
      invdxCo2,xComin1,xComin2)
case ("unlimit")
   call interpolation_unlimit(pwuse(igrid),ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,dxFi1,dxFi2,xFimin1,xFimin2, pwuseCo(igrid),dxCo1,dxCo2,&
      invdxCo1,invdxCo2,xComin1,xComin2)
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select

if (amrentropy) then
    call rhos_to_e(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,ixComin2,&
       ixComax1,ixComax2,pwuseCo(igrid)%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
    call conserve(ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,ixComin1,ixComin2,&
       ixComax1,ixComax2,pwuseCo(igrid)%w,pxCoarse(igrid)%x,patchfalse)
end if

end subroutine bc_prolong
!=============================================================================
subroutine interpolation_linear(pwFi,ixFimin1,ixFimin2,ixFimax1,ixFimax2,&
   dxFi1,dxFi2,xFimin1,xFimin2, pwCo,dxCo1,dxCo2,invdxCo1,invdxCo2,xComin1,&
   xComin2)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimax1,ixFimax2
double precision, intent(in) :: dxFi1,dxFi2, xFimin1,xFimin2,dxCo1,dxCo2,&
    invdxCo1,invdxCo2, xComin1,xComin2
type(walloc) :: pwCo, pwFi

integer :: ixCo1,ixCo2, jxCo1,jxCo2, hxCo1,hxCo2, ixFi1,ixFi2, ix1,ix2, iw,&
    idims
double precision :: xCo1,xCo2, xFi1,xFi2, eta1,eta2
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nwstart+1:nwstart+nwbc,ndim)
!-----------------------------------------------------------------------------
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1

   ! cell-centered coordinate for coarse cell
   xCo2=xComin2+(dble(ixCo2)-half)*dxCo2
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;ix2=2*int((ixFi2+ixMlo2)/2)-ixMlo2;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi(igrid)%dvolume(ixFi1,ixFi2) &
         /sum(pgeoFi(igrid)%dvolume(ix1:ix1+1,ixFi2))) 
      eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeoFi(igrid)%dvolume(ixFi1,ixFi2) &
         /sum(pgeoFi(igrid)%dvolume(ixFi1,ix2:ix2+1))) 

   end if

   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      hxCo2=ixCo2-kr(2,idims)
      jxCo1=ixCo1+kr(1,idims)
      jxCo2=ixCo2+kr(2,idims)

      do iw=nwstart+1,nwstart+nwbc
         slopeL=pwCo%w(ixCo1,ixCo2,iw)-pwCo%w(hxCo1,hxCo2,iw)
         slopeR=pwCo%w(jxCo1,jxCo2,iw)-pwCo%w(ixCo1,ixCo2,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idims)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), signR&
              *slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idims)=signR*max(zero,min(mcbeta*dabs(slopeR), mcbeta&
              *signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, (dabs(slopeR)&
              +two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idims)=signC*max(zero,min(dabs(slopeC), signC&
              *slopeL,signC*slopeR))
         end select
      end do
   end do

   ! Interpolate from coarse cell using limited slopes
   pwFi%w(ixFi1,ixFi2,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo1,ixCo2,nwstart&
      +1:nwstart+nwbc)+(slope(nwstart+1:nwstart+nwbc,1)*eta1)+(slope(nwstart&
      +1:nwstart+nwbc,2)*eta2)

end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_linear
!=============================================================================
subroutine interpolation_copy(pwFi,ixFimin1,ixFimin2,ixFimax1,ixFimax2,dxFi1,&
   dxFi2,xFimin1,xFimin2, pwCo,dxCo1,dxCo2,invdxCo1,invdxCo2,xComin1,xComin2)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimax1,ixFimax2
double precision, intent(in) :: dxFi1,dxFi2, xFimin1,xFimin2,dxCo1,dxCo2,&
    invdxCo1,invdxCo2, xComin1,xComin2
type(walloc) :: pwCo, pwFi

integer :: ixCo1,ixCo2, ixFi1,ixFi2
double precision :: xFi1,xFi2
!-----------------------------------------------------------------------------
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! Copy from coarse cell
   pwFi%w(ixFi1,ixFi2,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo1,ixCo2,nwstart&
      +1:nwstart+nwbc)

end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_copy
!=============================================================================
subroutine interpolation_unlimit(pwFi,ixFimin1,ixFimin2,ixFimax1,ixFimax2,&
   dxFi1,dxFi2,xFimin1,xFimin2, pwCo,dxCo1,dxCo2,invdxCo1,invdxCo2,xComin1,&
   xComin2)

integer, intent(in) :: ixFimin1,ixFimin2,ixFimax1,ixFimax2
double precision, intent(in) :: dxFi1,dxFi2, xFimin1,xFimin2, dxCo1,dxCo2,&
   invdxCo1,invdxCo2, xComin1,xComin2
type(walloc) :: pwCo, pwFi

integer :: ixCo1,ixCo2, jxCo1,jxCo2, hxCo1,hxCo2, ixFi1,ixFi2, ix1,ix2, idims
double precision :: xCo1,xCo2, xFi1,xFi2, eta1,eta2
double precision :: slope(nwstart+1:nwstart+nwbc,ndim)
!-----------------------------------------------------------------------------
do ixFi2 = ixFimin2,ixFimax2
   ! cell-centered coordinates of fine grid point
   xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

   ! indices of coarse cell which contains the fine cell
   ixCo2=int((xFi2-xComin2)*invdxCo2)+1

   ! cell-centered coordinate for coarse cell
   xCo2=xComin2+(dble(ixCo2)-half)*dxCo2
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;ix2=2*int((ixFi2+ixMlo2)/2)-ixMlo2;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi(igrid)%dvolume(ixFi1,ixFi2) &
         /sum(pgeoFi(igrid)%dvolume(ix1:ix1+1,ixFi2))) 
      eta2=(xFi2-xCo2)*invdxCo2 *two*(one-pgeoFi(igrid)%dvolume(ixFi1,ixFi2) &
         /sum(pgeoFi(igrid)%dvolume(ixFi1,ix2:ix2+1))) 
   end if

   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      hxCo2=ixCo2-kr(2,idims)
      jxCo1=ixCo1+kr(1,idims)
      jxCo2=ixCo2+kr(2,idims)

      ! get centered slope
      slope(nwstart+1:nwstart+nwbc,idims)=half*(pwCo%w(jxCo1,jxCo2,nwstart&
         +1:nwstart+nwbc)-pwCo%w(hxCo1,hxCo2,nwstart+1:nwstart+nwbc))
   end do

   ! Interpolate from coarse cell using centered slopes
   pwFi%w(ixFi1,ixFi2,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo1,ixCo2,nwstart&
      +1:nwstart+nwbc)+(slope(nwstart+1:nwstart+nwbc,1)*eta1)+(slope(nwstart&
      +1:nwstart+nwbc,2)*eta2)
end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixFimin1,ixFimin2,ixFimax1,&
      ixFimax2,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_unlimit
!=============================================================================
subroutine init_bc

integer :: dixBCo, interpolation_order
integer :: ixoldGmin1,ixoldGmin2,ixoldGmax1,ixoldGmax2, ixoldMmin1,ixoldMmin2,&
   ixoldMmax1,ixoldMmax2, nx1,nx2, nxCo1,nxCo2
!-----------------------------------------------------------------------------
ixMmin1=ixGmin1+dixB;ixMmin2=ixGmin2+dixB;ixMmax1=ixGmax1-dixB
ixMmax2=ixGmax2-dixB;
ixCoGmin1=1;ixCoGmin2=1;
ixCoGmax1=ixGmax1/2+dixB;ixCoGmax2=ixGmax2/2+dixB;
ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmax1=ixCoGmax1-dixB
ixCoMmax2=ixCoGmax2-dixB;

if (richardson) then
   ixoldGmin1=1;ixoldGmin2=1; ixoldGmax1=2*ixMmax1;ixoldGmax2=2*ixMmax2;
   ixoldMmin1=ixoldGmin1+dixB;ixoldMmin2=ixoldGmin2+dixB
   ixoldMmax1=ixoldGmax1-dixB;ixoldMmax2=ixoldGmax2-dixB;
end if

nx1=ixMmax1-ixMmin1+1;nx2=ixMmax2-ixMmin2+1;
nxCo1=nx1/2;nxCo2=nx2/2;

select case (typeghostfill)
case ("copy")
   interpolation_order=1
case ("linear","unlimit")
   interpolation_order=2
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select
dixBCo=int((dixB+1)/2)

if (dixBCo+interpolation_order-1>dixB) then
   call mpistop("interpolation order for prolongation in getbc to high")
end if


ixS_srl_min1(-1)=ixMmin1
ixS_srl_min1(0) =ixMmin1
ixS_srl_min1(1) =ixMmax1+1-dixB
ixS_srl_max1(-1)=ixMmin1-1+dixB
ixS_srl_max1(0) =ixMmax1
ixS_srl_max1(1) =ixMmax1

ixR_srl_min1(-1)=1
ixR_srl_min1(0) =ixMmin1
ixR_srl_min1(1) =ixMmax1+1
ixR_srl_max1(-1)=dixB
ixR_srl_max1(0) =ixMmax1
ixR_srl_max1(1) =ixGmax1


ixS_srl_min2(-1)=ixMmin2
ixS_srl_min2(0) =ixMmin2
ixS_srl_min2(1) =ixMmax2+1-dixB
ixS_srl_max2(-1)=ixMmin2-1+dixB
ixS_srl_max2(0) =ixMmax2
ixS_srl_max2(1) =ixMmax2

ixR_srl_min2(-1)=1
ixR_srl_min2(0) =ixMmin2
ixR_srl_min2(1) =ixMmax2+1
ixR_srl_max2(-1)=dixB
ixR_srl_max2(0) =ixMmax2
ixR_srl_max2(1) =ixGmax2


if (levmin/=levmax) then

   ixS_r_min1(-1)=ixCoMmin1
   ixS_r_min1(0) =ixCoMmin1
   ixS_r_min1(1) =ixCoMmax1+1-dixB
   ixS_r_max1(-1)=ixCoMmin1-1+dixB
   ixS_r_max1(0) =ixCoMmax1
   ixS_r_max1(1) =ixCoMmax1

   ixR_r_min1(0)=1
   ixR_r_min1(1)=ixMmin1
   ixR_r_min1(2)=ixMmin1+nxCo1
   ixR_r_min1(3)=ixMmax1+1
   ixR_r_max1(0)=dixB
   ixR_r_max1(1)=ixMmin1-1+nxCo1
   ixR_r_max1(2)=ixMmax1
   ixR_r_max1(3)=ixGmax1

   if (richardson) then
      ixS_old_min1(0)=ixoldMmin1
      ixS_old_min1(1)=ixoldMmin1
      ixS_old_min1(2)=ixoldMmin1+nx1-dixB
      ixS_old_min1(3)=ixoldMmax1+1-dixB
      ixS_old_max1(0)=ixoldMmin1-1+dixB
      ixS_old_max1(1)=ixoldMmin1-1+nx1+dixB
      ixS_old_max1(2)=ixoldMmax1
      ixS_old_max1(3)=ixoldMmax1

      ixR_old_min1(0)=1
      ixR_old_min1(1)=ixMmin1
      ixR_old_min1(2)=1
      ixR_old_min1(3)=ixMmax1+1
      ixR_old_max1(0)=dixB
      ixR_old_max1(1)=ixGmax1
      ixR_old_max1(2)=ixMmax1
      ixR_old_max1(3)=ixGmax1
   else
      ixS_p_min1(0)=ixMmin1-(interpolation_order-1)
      ixS_p_min1(1)=ixMmin1-(interpolation_order-1)
      ixS_p_min1(2)=ixMmin1+nxCo1-dixBCo-(interpolation_order-1)
      ixS_p_min1(3)=ixMmax1+1-dixBCo-(interpolation_order-1)
      ixS_p_max1(0)=ixMmin1-1+dixBCo+(interpolation_order-1)
      ixS_p_max1(1)=ixMmin1-1+nxCo1+dixBCo+(interpolation_order-1)
      ixS_p_max1(2)=ixMmax1+(interpolation_order-1)
      ixS_p_max1(3)=ixMmax1+(interpolation_order-1)

      ixR_p_min1(0)=ixCoMmin1-dixBCo-(interpolation_order-1)
      ixR_p_min1(1)=ixCoMmin1-(interpolation_order-1)
      ixR_p_min1(2)=ixCoMmin1-dixBCo-(interpolation_order-1)
      ixR_p_min1(3)=ixCoMmax1+1-(interpolation_order-1)
      ixR_p_max1(0)=dixB+(interpolation_order-1)
      ixR_p_max1(1)=ixCoMmax1+dixBCo+(interpolation_order-1)
      ixR_p_max1(2)=ixCoMmax1+(interpolation_order-1)
      ixR_p_max1(3)=ixCoMmax1+dixBCo+(interpolation_order-1)
   end if


   ixS_r_min2(-1)=ixCoMmin2
   ixS_r_min2(0) =ixCoMmin2
   ixS_r_min2(1) =ixCoMmax2+1-dixB
   ixS_r_max2(-1)=ixCoMmin2-1+dixB
   ixS_r_max2(0) =ixCoMmax2
   ixS_r_max2(1) =ixCoMmax2

   ixR_r_min2(0)=1
   ixR_r_min2(1)=ixMmin2
   ixR_r_min2(2)=ixMmin2+nxCo2
   ixR_r_min2(3)=ixMmax2+1
   ixR_r_max2(0)=dixB
   ixR_r_max2(1)=ixMmin2-1+nxCo2
   ixR_r_max2(2)=ixMmax2
   ixR_r_max2(3)=ixGmax2

   if (richardson) then
      ixS_old_min2(0)=ixoldMmin2
      ixS_old_min2(1)=ixoldMmin2
      ixS_old_min2(2)=ixoldMmin2+nx2-dixB
      ixS_old_min2(3)=ixoldMmax2+1-dixB
      ixS_old_max2(0)=ixoldMmin2-1+dixB
      ixS_old_max2(1)=ixoldMmin2-1+nx2+dixB
      ixS_old_max2(2)=ixoldMmax2
      ixS_old_max2(3)=ixoldMmax2

      ixR_old_min2(0)=1
      ixR_old_min2(1)=ixMmin2
      ixR_old_min2(2)=1
      ixR_old_min2(3)=ixMmax2+1
      ixR_old_max2(0)=dixB
      ixR_old_max2(1)=ixGmax2
      ixR_old_max2(2)=ixMmax2
      ixR_old_max2(3)=ixGmax2
   else
      ixS_p_min2(0)=ixMmin2-(interpolation_order-1)
      ixS_p_min2(1)=ixMmin2-(interpolation_order-1)
      ixS_p_min2(2)=ixMmin2+nxCo2-dixBCo-(interpolation_order-1)
      ixS_p_min2(3)=ixMmax2+1-dixBCo-(interpolation_order-1)
      ixS_p_max2(0)=ixMmin2-1+dixBCo+(interpolation_order-1)
      ixS_p_max2(1)=ixMmin2-1+nxCo2+dixBCo+(interpolation_order-1)
      ixS_p_max2(2)=ixMmax2+(interpolation_order-1)
      ixS_p_max2(3)=ixMmax2+(interpolation_order-1)

      ixR_p_min2(0)=ixCoMmin2-dixBCo-(interpolation_order-1)
      ixR_p_min2(1)=ixCoMmin2-(interpolation_order-1)
      ixR_p_min2(2)=ixCoMmin2-dixBCo-(interpolation_order-1)
      ixR_p_min2(3)=ixCoMmax2+1-(interpolation_order-1)
      ixR_p_max2(0)=dixB+(interpolation_order-1)
      ixR_p_max2(1)=ixCoMmax2+dixBCo+(interpolation_order-1)
      ixR_p_max2(2)=ixCoMmax2+(interpolation_order-1)
      ixR_p_max2(3)=ixCoMmax2+dixBCo+(interpolation_order-1)
   end if

end if

if (npe>1) then
   do i2=-1,1
   do i1=-1,1
      if (i1==0.and.i2==0) cycle

      call get_bc_comm_type(type_send_srl(i1,i2),ixS_srl_min1(i1),&
         ixS_srl_min2(i2),ixS_srl_max1(i1),ixS_srl_max2(i2),ixGmin1,ixGmin2,&
         ixGmax1,ixGmax2)
      call get_bc_comm_type(type_recv_srl(i1,i2),ixR_srl_min1(i1),&
         ixR_srl_min2(i2),ixR_srl_max1(i1),ixR_srl_max2(i2),ixGmin1,ixGmin2,&
         ixGmax1,ixGmax2)

      if (levmin==levmax) cycle

      call get_bc_comm_type(type_send_r(i1,i2),ixS_r_min1(i1),ixS_r_min2(i2),&
         ixS_r_max1(i1),ixS_r_max2(i2),ixCoGmin1,ixCoGmin2,ixCoGmax1,&
         ixCoGmax2)
      do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
         inc2=2*i2+ic2
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
         inc1=2*i1+ic1
         call get_bc_comm_type(type_recv_r(inc1,inc2),ixR_r_min1(inc1),&
            ixR_r_min2(inc2),ixR_r_max1(inc1),ixR_r_max2(inc2),ixGmin1,&
            ixGmin2,ixGmax1,ixGmax2)
         if (richardson) then
            call get_bc_comm_type(type_send_old(inc1,inc2),&
                ixS_old_min1(inc1),ixS_old_min2(inc2),ixS_old_max1(inc1),&
               ixS_old_max2(inc2),ixoldGmin1,ixoldGmin2,ixoldGmax1,ixoldGmax2)
            call get_bc_comm_type(type_recv_old(inc1,inc2),&
                ixR_old_min1(inc1),ixR_old_min2(inc2),ixR_old_max1(inc1),&
               ixR_old_max2(inc2),ixGmin1,ixGmin2,ixGmax1,ixGmax2)
         else
            call get_bc_comm_type(type_send_p(inc1,inc2),ixS_p_min1(inc1),&
               ixS_p_min2(inc2),ixS_p_max1(inc1),ixS_p_max2(inc2),ixGmin1,&
               ixGmin2,ixGmax1,ixGmax2)
            call get_bc_comm_type(type_recv_p(inc1,inc2),ixR_p_min1(inc1),&
               ixR_p_min2(inc2),ixR_p_max1(inc1),ixR_p_max2(inc2),ixCoGmin1,&
               ixCoGmin2,ixCoGmax1,ixCoGmax2)
         end if
      end do
      end do
   end do
   end do
end if

end subroutine init_bc
!=============================================================================
subroutine get_bc_comm_type(comm_type,ixmin1,ixmin2,ixmax1,ixmax2,ixGmin1,&
   ixGmin2,ixGmax1,ixGmax2)

integer, intent(inout) :: comm_type
integer, intent(in) :: ixmin1,ixmin2,ixmax1,ixmax2, ixGmin1,ixGmin2,ixGmax1,&
   ixGmax2

integer, dimension(ndim+1) :: size, subsize, start
!-----------------------------------------------------------------------------
size(1)=ixGmax1;size(2)=ixGmax2;
size(ndim+1)=nw
subsize(1)=ixmax1-ixmin1+1;subsize(2)=ixmax2-ixmin2+1;
subsize(ndim+1)=nwbc
start(1)=ixmin1-1;start(2)=ixmin2-1;
start(ndim+1)=nwstart

call MPI_TYPE_CREATE_SUBARRAY(ndim+1,size,subsize,start,MPI_ORDER_FORTRAN,&
    MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
call MPI_TYPE_COMMIT(comm_type,ierrmpi)

end subroutine get_bc_comm_type
!=============================================================================
subroutine put_bc_comm_types
!-----------------------------------------------------------------------------
do i2=-1,1
do i1=-1,1
   if (i1==0.and.i2==0) cycle

   call MPI_TYPE_FREE(type_send_srl(i1,i2),ierrmpi)
   call MPI_TYPE_FREE(type_recv_srl(i1,i2),ierrmpi)

   if (levmin==levmax) cycle

   call MPI_TYPE_FREE(type_send_r(i1,i2),ierrmpi)
   do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
      inc2=2*i2+ic2
   do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
      inc1=2*i1+ic1
      call MPI_TYPE_FREE(type_recv_r(inc1,inc2),ierrmpi)
      if (richardson) then
         call MPI_TYPE_FREE(type_send_old(inc1,inc2),ierrmpi)
         call MPI_TYPE_FREE(type_recv_old(inc1,inc2),ierrmpi)
      else
         call MPI_TYPE_FREE(type_send_p(inc1,inc2),ierrmpi)
         call MPI_TYPE_FREE(type_recv_p(inc1,inc2),ierrmpi)
      end if
   end do
   end do
end do
end do

end subroutine put_bc_comm_types
!=============================================================================
subroutine pole_copy(pwrecv,ixRmin1,ixRmin2,ixRmax1,ixRmax2,pwsend,ixSmin1,&
   ixSmin2,ixSmax1,ixSmax2)

integer, intent(in) :: ixRmin1,ixRmin2,ixRmax1,ixRmax2, ixSmin1,ixSmin2,&
   ixSmax1,ixSmax2
type(walloc) :: pwrecv, pwsend

integer :: iw, iB
!-----------------------------------------------------------------------------
select case (ipole)
case (1)
   iside=int((i1+3)/2)
   iB=2*(1-1)+iside
   do iw=nwstart+1,nwstart+nwbc
      select case (typeB(iw,iB))
      case ("symm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw) = pwsend%w&
            (ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw) =-pwsend%w&
            (ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,iw)
      case default
         call mpistop("Boundary condition at pole should be symm or asymm")
      end select
   end do 
case (2)
   iside=int((i2+3)/2)
   iB=2*(2-1)+iside
   do iw=nwstart+1,nwstart+nwbc
      select case (typeB(iw,iB))
      case ("symm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw) = pwsend%w&
            (ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw) =-pwsend%w&
            (ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,iw)
      case default
         call mpistop("Boundary condition at pole should be symm or asymm")
      end select
   end do 
end select

end subroutine pole_copy
!=============================================================================
subroutine fix_auxiliary

integer :: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
   do i2=-1,1
   do i1=-1,1
      if (i1==0.and.i2==0) cycle

      ixmin1=ixR_srl_min1(i1);ixmin2=ixR_srl_min2(i2);ixmax1=ixR_srl_max1(i1)
      ixmax2=ixR_srl_max2(i2);
      if(.not.slab)mygeo=>pgeoFi(igrid)
      call getaux(.true.,pwuse(igrid)%w,px(igrid)%x,ixGmin1,ixGmin2,ixGmax1,&
         ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,"bc")
   end do
   end do
end do

end subroutine fix_auxiliary
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine getbc
!=============================================================================
subroutine physbound(i1,i2,igrid,isphysbound)
use mod_forest
include 'amrvacdef.f'

integer, intent(in)  :: i1,i2, igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: level, ig1,ig2, ign1,ign2
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
ig1 = tree%node%ig1; ig2 = tree%node%ig2;

ign1 = ig1 + i1; ign2 = ig2 + i2;
if (ign1 .gt. ng1(level) .or. ign1 .lt. 1.or.ign2 .gt. ng2(level) &
   .or. ign2 .lt. 1) isphysbound = .true.

end subroutine physbound
!=============================================================================
