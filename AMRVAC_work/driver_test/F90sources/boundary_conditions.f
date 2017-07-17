!=============================================================================
subroutine bc_phys(iside,idims,time,qdt,w,x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
   ixBmin1,ixBmin2,ixBmax1,ixBmax2)

include 'amrvacdef.f'

integer, intent(in) :: iside, idims, ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
   ixBmin2,ixBmax1,ixBmax2
double precision, intent(in) :: time,qdt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
double precision :: wtmp(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nwflux)

integer :: iw, iB, ix1,ix2, ixImin1,ixImin2,ixImax1,ixImax2, ixMmin1,ixMmin2,&
   ixMmax1,ixMmax2
!-----------------------------------------------------------------------------
select case (idims)
case (1)
   if (iside==2) then
      ! maximal boundary
      iB=ismax1
      ixImin1=ixBmax1+1-dixB;ixImin2=ixBmin2;
      ixImax1=ixBmax1;ixImax2=ixBmax2;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = w(ixImin1-1:ixImin1-dixB:&
               -1,ixImin2:ixImax2,iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) =-w(ixImin1-1:ixImin1-dixB:&
               -1,ixImin2:ixImax2,iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,ixImin2:ixImax2,iw) = w(ixImin1-1,ixImin2:ixImax2,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,iw) = max(w(ixImin1&
                     -1,ixImin2:ixImax2,iw),zero)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,iw) = w(ixImin1-1,ixImin2:ixImax2,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,iw) = max(w(ixImin1&
                     -1,ixImin2:ixImax2,iw), &
                                           w(ixImin1-1,ixImin2:ixImax2,iw)&
                                              *ratebdflux)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,iw) = w(ixImin1-1,ixImin2:ixImax2,iw)
              end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = - w(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      iB=ismin1
      ixImin1=ixBmin1;ixImin2=ixBmin2;
      ixImax1=ixBmin1-1+dixB;ixImax2=ixBmax2;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = w(ixImax1+dixB:ixImax1&
               +1:-1,ixImin2:ixImax2,iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) =-w(ixImax1+dixB:ixImax1&
               +1:-1,ixImin2:ixImax2,iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,ixImin2:ixImax2,iw) = w(ixImax1+1,ixImin2:ixImax2,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
               do ix1=ixImin1,ixImax1
                 w(ix1,ixImin2:ixImax2,iw) = min(w(ixImax1+1,ixImin2:ixImax2,&
                    iw),zero)
               end do
            else
               do ix1=ixImin1,ixImax1
                 w(ix1,ixImin2:ixImax2,iw) = w(ixImax1+1,ixImin2:ixImax2,iw)
               end do
            end if
         case("limitinflow")
            if (iw==1+1)then
               do ix1=ixImin1,ixImax1
                 w(ix1,ixImin2:ixImax2,iw) = min(w(ixImax1+1,ixImin2:ixImax2,&
                    iw), &
                                          w(ixImax1+1,ixImin2:ixImax2,iw)&
                                             *ratebdflux)
               end do
            else
               do ix1=ixImin1,ixImax1
                 w(ix1,ixImin2:ixImax2,iw) = w(ixImax1+1,ixImin2:ixImax2,iw)
               end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = - w(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if 
case (2)
   if (iside==2) then
      ! maximal boundary
      iB=ismax2
      ixImin1=ixBmin1;ixImin2=ixBmax2+1-dixB;
      ixImax1=ixBmax1;ixImax2=ixBmax2;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = w(ixImin1:ixImax1,ixImin2&
               -1:ixImin2-dixB:-1,iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) =-w(ixImin1:ixImax1,ixImin2&
               -1:ixImin2-dixB:-1,iw)
         case ("cont")
            do ix2=ixImin2,ixImax2
               w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImin2-1,iw)
            end do
         case("noinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,iw) = max(w(ixImin1:ixImax1,ixImin2&
                     -1,iw),zero)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImin2-1,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,iw) = max(w(ixImin1:ixImax1,ixImin2&
                     -1,iw), &
                                           w(ixImin1:ixImax1,ixImin2&
                                              -1,iw)*ratebdflux)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImin2-1,iw)
              end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = - w(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   else
      ! minimal boundary
      iB=ismin2
      ixImin1=ixBmin1;ixImin2=ixBmin2;
      ixImax1=ixBmax1;ixImax2=ixBmin2-1+dixB;
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux
         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = w(ixImin1:ixImax1,ixImax2&
               +dixB:ixImax2+1:-1,iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) =-w(ixImin1:ixImax1,ixImax2&
               +dixB:ixImax2+1:-1,iw)
         case ("cont")
            do ix2=ixImin2,ixImax2
               w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImax2+1,iw)
            end do
         case("noinflow")
            if (iw==1+2)then
               do ix2=ixImin2,ixImax2
                 w(ixImin1:ixImax1,ix2,iw) = min(w(ixImin1:ixImax1,ixImax2&
                    +1,iw),zero)
               end do
            else
               do ix2=ixImin2,ixImax2
                 w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImax2+1,iw)
               end do
            end if
         case("limitinflow")
            if (iw==1+2)then
               do ix2=ixImin2,ixImax2
                 w(ixImin1:ixImax1,ix2,iw) = min(w(ixImin1:ixImax1,ixImax2&
                    +1,iw), &
                                          w(ixImin1:ixImax1,ixImax2&
                                             +1,iw)*ratebdflux)
               end do
            else
               do ix2=ixImin2,ixImax2
                 w(ixImin1:ixImax1,ix2,iw) = w(ixImin1:ixImax1,ixImax2+1,iw)
               end do
            end if
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("character")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,iw) = - w(ixImin1:ixImax1,&
               ixImin2:ixImax2,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select
      end do
   end if 
end select

! do special case AFTER all normal cases are set
!do iw=1,nwflux+nwaux
! opedit: iw==0 since this breaks fewest of setups.
if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
  call specialbound_usr(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixImin1,ixImin2,&
     ixImax1,ixImax2,0,iB,w,x)
end if

!end do

end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,pwuse)

include 'amrvacdef.f'

double precision, intent(in)   :: time
integer, intent(in)            :: ixGmin1,ixGmin2,ixGmax1,ixGmax2
type(walloc), dimension(ngridshi)          :: pwuse

! .. local ..
integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmax1,ixOmax2,level
!----------------------------------------------------------------------------
ixOmin1=ixGmin1+dixB;ixOmin2=ixGmin2+dixB;ixOmax1=ixGmax1-dixB
ixOmax2=ixGmax2-dixB;

do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
!do iigrid=1,igridstail; igrid=igrids(iigrid);
   dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      myB0_face1 => pB0_face1(igrid)
      myB0_face2 => pB0_face2(igrid)
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   level=node(plevel_,igrid)
   saveigrid=igrid
   call bc_int(level,time,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,pwuse(igrid)%w,px(igrid)%x)
end do

      
end subroutine getintbc
!=============================================================================
