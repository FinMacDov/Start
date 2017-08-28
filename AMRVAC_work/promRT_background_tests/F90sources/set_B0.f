!=============================================================================
subroutine set_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

call set_B0_cell(pB0_cell(igrid)%w,px(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,&
   ixGlo1,ixGlo2,ixGhi1,ixGhi2)

call set_B0_face(igrid,px(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,&
   ixMhi1,ixMhi2)


end subroutine set_B0_grid
!=============================================================================
subroutine set_B0_cell(wB0,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2)

include 'amrvacdef.f'

integer, intent(in):: ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2
double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
!-----------------------------------------------------------------------------
wB0(ixmin1:ixmax1,ixmin2:ixmax2,1:ndir)=zero

! approximate cell-averaged B0 as cell-centered B0
select case (typeaxial)
case ("spherical")
   
   if (dabs(Bdip)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=2.0d0*Bdip*dcos(x(ixmin1:ixmax1,&
         ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=Bdip*dsin(x(ixmin1:ixmax1,&
         ixmin2:ixmax2,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**3
   end if

   if (abs(Bquad)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,1) &
           +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
              2)))/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)&
         +Bquad*dsin(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,2))&
         /x(ixmin1:ixmax1,ixmin2:ixmax2,1)**4
   end if
   if (abs(Boct)>smalldouble) then
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,1)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,1) &
                   +Boct*(10.0d0*dcos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,2))&
                      -2.0d0) &
                        *dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,2))&
                           /x(ixmin1:ixmax1,ixmin2:ixmax2,1)**5
      wB0(ixmin1:ixmax1,ixmin2:ixmax2,2)=wB0(ixmin1:ixmax1,ixmin2:ixmax2,2) &
                   +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
                      ixmin2:ixmax2,2))) &
                        *dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,2))&
                           /x(ixmin1:ixmax1,ixmin2:ixmax2,1)**5
   end if
  
end select

if(dabs(Busr)/=zero) then
   call specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
      ixmax2,x,wB0)
end if

end subroutine set_B0_cell
!=============================================================================
subroutine set_B0_face(igrid,x,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,&
   ixmax1,ixmax2)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,&
   ixmax1,ixmax2
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),dx1,dx2,xmin1,&
   xmin2,xshift1,xshift2
integer :: idims, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ix, idims2
!-----------------------------------------------------------------------------
dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);
xmin1=rnode(rpxmin1_,igrid);xmin2=rnode(rpxmin2_,igrid);


do idims=1,ndim
   ixCmin1=ixmin1-kr(1,idims);ixCmin2=ixmin2-kr(2,idims); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims));
   do idims2=1,ndim
      select case(idims2)
      case(1)
        do ix = ixCmin1,ixCmax1
          xC(ix,ixCmin2:ixCmax2,1)=xmin1+(dble(ix-dixB)-xshift1)*dx1
        end do
      case(2)
        do ix = ixCmin2,ixCmax2
          xC(ixCmin1:ixCmax1,ix,2)=xmin2+(dble(ix-dixB)-xshift2)*dx2
        end do
      end select
   end do


   select case(idims)
   case (1)
      call set_B0_cell(pB0_face1(igrid)%w,xC,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixCmin1,ixCmin2,ixCmax1,ixCmax2) 
   case (2)
      call set_B0_cell(pB0_face2(igrid)%w,xC,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixCmin1,ixCmin2,ixCmax1,ixCmax2) 
   end select
end do

end subroutine set_B0_face
!=============================================================================
subroutine alloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

allocate(pB0_cell(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir))
allocate(pB0_face1(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir))
allocate(pB0_face2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir))

end subroutine alloc_B0_grid
!=============================================================================
subroutine dealloc_B0_grid(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------

deallocate(pB0_cell(igrid)%w)
deallocate(pB0_face1(igrid)%w)
deallocate(pB0_face2(igrid)%w)

end subroutine dealloc_B0_grid
!=============================================================================
