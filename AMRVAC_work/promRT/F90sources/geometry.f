!=============================================================================
subroutine set_pole

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
select case (typeaxial)
case ("spherical") 
case ("cylindrical")
  if (1==-1.and.periodB(1)) then
    if(mod(ng1(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
  if (2==-1.and.periodB(2)) then
    if(mod(ng2(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
end select

end subroutine set_pole
!=============================================================================
subroutine getgridgeo(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: ixmin1,ixmin2,ixmax1,ixmax2, ixCoGmin1,ixCoGmin2,ixCoGmax1,&
   ixCoGmax2, ixCoMmin1,ixCoMmin2,ixCoMmax1,ixCoMmax2, ixComin1,ixComin2,&
   ixComax1,ixComax2, ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmax1,ixCoCoGmax2,&
    ixGextmin1,ixGextmin2,ixGextmax1,ixGextmax2
double precision :: xmin1,xmin2, dx1,dx2
!-----------------------------------------------------------------------------
!ix^L=ixM^LL^LADD1;
ixmin1=ixGlo1+1;ixmin2=ixGlo2+1;ixmax1=ixGhi1-1;ixmax2=ixGhi2-1;
if (2*int(dixB/2)==dixB) then
   ixGextmin1=ixGlo1;ixGextmin2=ixGlo2;ixGextmax1=ixGhi1;ixGextmax2=ixGhi2;
else
   ixGextmin1=ixGlo1-1;ixGextmin2=ixGlo2-1;ixGextmax1=ixGhi1+1
   ixGextmax2=ixGhi2+1;
end if


allocate(pgeo(igrid)%surfaceC1(ixmin1-1:ixmax1,ixmin2:ixmax2),&
   pgeo(igrid)%surfaceC2(ixmin1:ixmax1,ixmin2-1:ixmax2), pgeo(igrid)%surface1&
   (ixmin1-1:ixmax1,ixmin2:ixmax2),pgeo(igrid)%surface2(ixmin1:ixmax1,ixmin2&
   -1:ixmax2), pgeo(igrid)%dvolume(ixGextmin1:ixGextmax1,&
   ixGextmin2:ixGextmax2), pgeo(igrid)%x(ixGextmin1:ixGextmax1,&
   ixGextmin2:ixGextmax2,1:ndim),pgeo(igrid)%dx(ixGextmin1:ixGextmax1,&
   ixGextmin2:ixGextmax2,1:ndim))

dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);
xmin1=rnode(rpxmin1_,igrid);xmin2=rnode(rpxmin2_,igrid);

call fillgeo(pgeo(igrid),ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGextmin1,ixGextmin2,&
   ixGextmax1,ixGextmax2,xmin1,xmin2,dx1,dx2,.false.)

if (errorestimate==1) then
   ixCoGmin1=1;ixCoGmin2=1; ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmin2=ixCoGmin2;ixGextmax1=ixCoGmax1
      ixGextmax2=ixCoGmax2;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmin2=ixCoGmin2-1;ixGextmax1=ixCoGmax1+1
      ixGextmax2=ixCoGmax2+1;
   end if
   ixCoMmin1=ixCoGmin1+dixB;ixCoMmin2=ixCoGmin2+dixB;ixCoMmax1=ixCoGmax1-dixB
   ixCoMmax2=ixCoGmax2-dixB;
   ixComin1=ixCoMmin1-1;ixComin2=ixCoMmin2-1;ixComax1=ixCoMmax1+1
   ixComax2=ixCoMmax2+1;
   ixComin1=ixCoGmin1+1;ixComin2=ixCoGmin2+1;ixComax1=ixCoGmax1-1
   ixComax2=ixCoGmax2-1;

   allocate(pgeoCoarse(igrid)%surfaceC1(ixComin1-1:ixComax1,&
      ixComin2:ixComax2),pgeoCoarse(igrid)%surfaceC2(ixComin1:ixComax1,&
      ixComin2-1:ixComax2), pgeoCoarse(igrid)%surface1(ixComin1&
      -1:ixComax1,ixComin2:ixComax2),pgeoCoarse(igrid)%surface2&
      (ixComin1:ixComax1,ixComin2-1:ixComax2), pgeoCoarse(igrid)%dvolume&
      (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2), pgeoCoarse(igrid)%x&
      (ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1:ndim),&
      pgeoCoarse(igrid)%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
      1:ndim))

   dx1=two*rnode(rpdx1_,igrid);dx2=two*rnode(rpdx2_,igrid);

   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
      ixGextmin1,ixGextmin2,ixGextmax1,ixGextmax2,xmin1,xmin2,dx1,dx2,.false.)

   ixCoCoGmin1=1;ixCoCoGmin2=1; ixCoCoGmax1=ixCoGmax1/2+dixB
   ixCoCoGmax2=ixCoGmax2/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoCoGmin1;ixGextmin2=ixCoCoGmin2;ixGextmax1=ixCoCoGmax1
      ixGextmax2=ixCoCoGmax2;
   else
      ixGextmin1=ixCoCoGmin1-1;ixGextmin2=ixCoCoGmin2-1
      ixGextmax1=ixCoCoGmax1+1;ixGextmax2=ixCoCoGmax2+1;
   end if


   allocate(pgeoCoCo(igrid)%dvolume(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2))

   dx1=4.0d0*rnode(rpdx1_,igrid);dx2=4.0d0*rnode(rpdx2_,igrid);

   call fillgeo(pgeoCoCo(igrid),ixCoCoGmin1,ixCoCoGmin2,ixCoCoGmax1,&
      ixCoCoGmax2,ixGextmin1,ixGextmin2,ixGextmax1,ixGextmax2,xmin1,xmin2,dx1,&
      dx2,.true.)
else
   ixCoGmin1=1;ixCoGmin2=1; ixCoGmax1=ixGhi1/2+dixB;ixCoGmax2=ixGhi2/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmin2=ixCoGmin2;ixGextmax1=ixCoGmax1
      ixGextmax2=ixCoGmax2;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmin2=ixCoGmin2-1;ixGextmax1=ixCoGmax1+1
      ixGextmax2=ixCoGmax2+1;
   end if


   allocate(pgeoCoarse(igrid)%dvolume(ixGextmin1:ixGextmax1,&
      ixGextmin2:ixGextmax2))

   dx1=two*rnode(rpdx1_,igrid);dx2=two*rnode(rpdx2_,igrid);

   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmin2,ixCoGmax1,ixCoGmax2,&
      ixGextmin1,ixGextmin2,ixGextmax1,ixGextmax2,xmin1,xmin2,dx1,dx2,.true.)
end if

end subroutine getgridgeo
!=============================================================================
subroutine putgridgeo(igrid)

  include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (errorestimate==1) then
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surfaceC2,pgeo&
      (igrid)%surface1,pgeo(igrid)%surface2,pgeo(igrid)%dvolume,&
      pgeo(igrid)%dx,pgeo(igrid)%x, pgeoCoarse(igrid)%surfaceC1,&
      pgeoCoarse(igrid)%surfaceC2,pgeoCoarse(igrid)%surface1,&
      pgeoCoarse(igrid)%surface2,pgeoCoarse(igrid)%dvolume,&
      pgeoCoarse(igrid)%dx,pgeoCoCo(igrid)%dvolume)
   
else
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surfaceC2,pgeo&
      (igrid)%surface1,pgeo(igrid)%surface2,pgeo(igrid)%dvolume,&
      pgeo(igrid)%dx,pgeo(igrid)%x,pgeoCoarse(igrid)%dvolume)
end if

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(pgeogrid,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixGextmin1,&
   ixGextmin2,ixGextmax1,ixGextmax2,xmin1,xmin2,dx1,dx2,need_only_volume)

include 'amrvacdef.f'

type(geoalloc) :: pgeogrid
integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGextmin1,ixGextmin2,&
   ixGextmax1,ixGextmax2
double precision, intent(in) :: xmin1,xmin2, dx1,dx2
logical, intent(in) :: need_only_volume

integer :: idims, ix, ixMmin1,ixMmin2,ixMmax1,ixMmax2, ixmin1,ixmin2,ixmax1,&
   ixmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2
double precision :: x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,ndim) 
!-----------------------------------------------------------------------------
ixMmin1=ixGmin1+dixB;ixMmin2=ixGmin2+dixB;ixMmax1=ixGmax1-dixB
ixMmax2=ixGmax2-dixB;
!ix^L=ixM^L^LADD1;
ixmin1=ixGmin1+1;ixmin2=ixGmin2+1;ixmax1=ixGmax1-1;ixmax2=ixGmax2-1;

select case (typeaxial)
case ("slabtest")

   pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2) = dx1*dx2

   if (need_only_volume) return

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)= dx2
   pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = dx2
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dx1
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dx1
   

   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)=dx1
   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)=dx2;

case ("spherical")

   do idims=1,min(ndim,2)
      select case(idims)
      case(1)
         do ix = ixGextmin1,ixGextmax1
            x(ix,ixGextmin2:ixGextmax2,1)=xmin1+(dble(ix-dixB)-half)*dx1
         end do
      case(2)
         do ix = ixGextmin2,ixGextmax2
            x(ixGextmin1:ixGextmax1,ix,2)=xmin2+(dble(ix-dixB)-half)*dx2
         end do
      end select
   end do


   if(typespherical==0) then
     pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2)&
        =(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)**2&
        +dx1**2/12.0d0)*dx1  &
              *two*dabs(dsin(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,&
                 2)))*dsin(half*dx2)
   else
     pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2)&
        =(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)**2)*dx1  &
              *dabs(dsin(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)))&
                 *dx2
   endif



   if (need_only_volume) return

   pgeogrid%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1:ndim)&
      =x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1:ndim)

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;

   if(typespherical==0) then
       pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(x(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,1)+half*dx1)**2  &
              *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))*dsin(half*dx2)
   else
       pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(x(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,1)+half*dx1)**2  &
              *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))*dx2
   endif

   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)*dx1 &
              *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+half*dx2)

   



   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   if(typespherical==0) then
       pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=x(ixCmin1:ixCmax1,&
          ixCmin2:ixCmax2,1)**2  *two*dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          2))*dsin(half*dx2)
   else
      pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=x(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)**2  *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))&
         *dx2
   endif


   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)*dx1 &
              *dsin(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2))

   

   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)=dx1


    pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,2)&
       =x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)*dx2
   

case ("cylindrical")


   do ix = ixGextmin1,ixGextmax1
      x(ix,ixGextmin2:ixGextmax2,1)=xmin1+(dble(ix-dixB)-half)*dx1
   end do

   pgeogrid%dvolume(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2)=dabs(half&
      *((x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)+half*dx1)**2&
      -(x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)-half*dx1)**2))*dx2



   if (need_only_volume) return

   pgeogrid%x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)&
      =x(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;

   pgeogrid%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dabs(x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1)+half*dx1)*dx2
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   if (-2==2) pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*dx1
   if (-1==2) pgeogrid%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dx1
   

   ixCmin1=ixmin1-kr(1,1);ixCmin2=ixmin2-kr(2,1); ixCmax1=ixmax1
   ixCmax2=ixmax2;
 !!pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*dx2
   pgeogrid%surface1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dabs(x(ixCmin1:ixCmax1,&
      ixCmin2:ixCmax2,1))*dx2
   
   ixCmin1=ixmin1-kr(1,2);ixCmin2=ixmin2-kr(2,2); ixCmax1=ixmax1
   ixCmax2=ixmax2;
   if (-2==2) pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)*dx1
   if (-1==2) pgeogrid%surface2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=dx1
   


   pgeogrid%dx(ixGextmin1:ixGextmax1,ixGextmin2:ixGextmax2,1)=dx1


   
   


case default

   call mpistop("Sorry, typeaxial unknown")
   
end select

end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

include 'amrvacdef.f'

integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,ixmax2, idir
double precision :: q(ixImin1:ixImax1,ixImin2:ixImax2), gradq(ixImin1:ixImax1,&
   ixImin2:ixImax2)

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2),invdx
integer :: jxmin1,jxmin2,jxmax1,jxmax2, hxmin1,hxmin2,hxmax1,hxmax2, ixCmin1,&
   ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2 

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then

   jxmin1=ixmin1+kr(idir,1);jxmin2=ixmin2+kr(idir,2);jxmax1=ixmax1+kr(idir,1)
   jxmax2=ixmax2+kr(idir,2);
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
   hxmax2=ixmax2-kr(idir,2);
   gradq(ixmin1:ixmax1,ixmin2:ixmax2) = half*(q(jxmin1:jxmax1,jxmin2:jxmax2)&
      -q(hxmin1:hxmax1,hxmin2:hxmax2))*invdx

else
   hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
   hxmax2=ixmax2-kr(idir,2);
   ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmax1=ixmax1;ixCmax2=ixmax2;
   jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
   jxCmax1=ixCmax1+kr(idir,1);jxCmax2=ixCmax2+kr(idir,2);
   select case(idir)
   case(1)
      qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*half*(q(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         +q(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
      gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qC(ixmin1:ixmax1,ixmin2:ixmax2)&
         -qC(hxmin1:hxmax1,hxmin2:hxmax2))/mygeo%dvolume(ixmin1:ixmax1,&
         ixmin2:ixmax2)
      ! Substract difference divergence and gradient
 !gradq(ixmin1:ixmax1,ixmin2:ixmax2)=gradq(ixmin1:ixmax1,ixmin2:ixmax2)-q(ixmin1:ixmax1,ixmin2:ixmax2) &
 !*(mygeo%surfaceC1(ixmin1:ixmax1,ixmin2:ixmax2)-mygeo%surfaceC1(hxmin1:hxmax1,hxmin2:hxmax2)) &
      !              /mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2) 
   case(2)
      qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC2(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*half*(q(ixCmin1:ixCmax1,ixCmin2:ixCmax2)&
         +q(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
      gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qC(ixmin1:ixmax1,ixmin2:ixmax2)&
         -qC(hxmin1:hxmax1,hxmin2:hxmax2))/mygeo%dvolume(ixmin1:ixmax1,&
         ixmin2:ixmax2)
      ! Substract difference divergence and gradient
 !gradq(ixmin1:ixmax1,ixmin2:ixmax2)=gradq(ixmin1:ixmax1,ixmin2:ixmax2)-q(ixmin1:ixmax1,ixmin2:ixmax2) &
 !*(mygeo%surfaceC2(ixmin1:ixmax1,ixmin2:ixmax2)-mygeo%surfaceC2(hxmin1:hxmax1,hxmin2:hxmax2)) &
      !              /mygeo%dvolume(ixmin1:ixmax1,ixmin2:ixmax2) 
   end select
end if

end subroutine gradient
!=============================================================================
subroutine gradientS(q,ixImin1,ixImin2,ixImax1,ixImax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

include 'amrvacdef.f'

integer :: ixImin1,ixImin2,ixImax1,ixImax2, ixmin1,ixmin2,ixmax1,ixmax2, idir
double precision :: q(ixImin1:ixImax1,ixImin2:ixImax2), gradq(ixImin1:ixImax1,&
   ixImin2:ixImax2)
double precision :: dxdim

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2)
double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2):: qL,qR,dqC,ldq,&
   invdx
integer                          :: hxmin1,hxmin2,hxmax1,hxmax2,ixCmin1,&
   ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,gxCmin1,gxCmin2,&
   gxCmax1,gxCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,idummy
character*79 :: savetypelimiter,savetypegradlimiter,save2typelimiter
!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
hxmin1=ixmin1-kr(idir,1);hxmin2=ixmin2-kr(idir,2);hxmax1=ixmax1-kr(idir,1)
hxmax2=ixmax2-kr(idir,2);
ixCmin1=hxmin1;ixCmin2=hxmin2;ixCmax1=ixmax1;ixCmax2=ixmax2;
jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
jxCmax1=ixCmax1+kr(idir,1);jxCmax2=ixCmax2+kr(idir,2);
gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2);gxCmax1=jxCmax1
gxCmax2=jxCmax2;
hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
hxCmax1=gxCmax1+kr(idir,1);hxCmax2=gxCmax2+kr(idir,2);
idummy=0

savetypelimiter=typelimiter
savetypegradlimiter=typegradlimiter
! set the gradient limiter here
typelimiter=typegradlimiter
qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(hxCmin1:hxCmax1,hxCmin2:hxCmax2)
qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = q(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
if (typelimiter/='ppm') then
   dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2)&
      -qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
   save2typelimiter=typelimiter
   if(save2typelimiter=='koren') typelimiter='korenL'
   if(save2typelimiter=='cada')  typelimiter='cadaL'
   if(save2typelimiter=='cada3') typelimiter='cada3L'
   dxdim=dxlevel(idir)
   call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
      gxCmax1,gxCmax2,idummy,idir,ldq,dxdim)
   qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) &
      + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
   if(save2typelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
        gxCmax1,gxCmax2,idummy,idir,ldq,dxdim)
   endif
   if(save2typelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
        gxCmax1,gxCmax2,idummy,idir,ldq,dxdim)
   endif
   if(save2typelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
        gxCmax1,gxCmax2,idummy,idir,ldq,dxdim)
   endif
   typelimiter=save2typelimiter
   qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) &
      - half*ldq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
else
   call PPMlimitervar(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,ixMhi2,&
      idir,q,q,qL,qR)
endif
! set the method limiter back
typelimiter=savetypelimiter
typegradlimiter=savetypegradlimiter

if (slab) then
   gradq(ixmin1:ixmax1,ixmin2:ixmax2)=half*(qR(ixmin1:ixmax1,ixmin2:ixmax2)&
      -qL(hxmin1:hxmax1,hxmin2:hxmax2))*invdx
else
   select case(idir)
   case(1)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qR(ixmin1:ixmax1,ixmin2:ixmax2)&
       -qL(hxmin1:hxmax1,hxmin2:hxmax2))/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       idir) 
   case(2)
    gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(qR(ixmin1:ixmax1,ixmin2:ixmax2)&
       -qL(hxmin1:hxmax1,hxmin2:hxmax2))/mygeo%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       idir) 
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,divq)

! Calculate divergence of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
    divq(ixImin1:ixImax1,ixImin2:ixImax2)

double precision :: qC(ixImin1:ixImax1,ixImin2:ixImax2), invdx(1:ndim)
integer :: jxOmin1,jxOmin2,jxOmax1,jxOmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
    ixCmin1,ixCmin2,ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, idims,&
    ixmin1,ixmin2,ixmax1,ixmax2 
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
if (slab) then
  do idims=1,ndim

     jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
     jxOmax1=ixOmax1+kr(idims,1);jxOmax2=ixOmax2+kr(idims,2);
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+half*(qvec(jxOmin1:jxOmax1,jxOmin2:jxOmax2,idims)&
        -qvec(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idims))*invdx(idims)

  end do
else
  do idims=1,ndim
     hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
     hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
     ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
     jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
     jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
     select case(idims)
     case(1)
        qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC1(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)&
           +qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,idims))
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
           -qC(hxOmin1:hxOmax1,hxOmin2:hxOmax2) 
     case(2)
        qC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC2(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*half*(qvec(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idims)&
           +qvec(jxCmin1:jxCmax1,jxCmin2:jxCmax2,idims))
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qC(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
           -qC(hxOmin1:hxOmax1,hxOmin2:hxOmax2) 
      end select
  end do
  divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
     /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
end if


end subroutine divvector 
!=============================================================================
subroutine curlvector(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,curlvec,idirmin,idirmin0,ndir0)

! Calculate curl of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   idirmin,ixmin1,ixmin2,ixmax1,ixmax2,idir,jdir,kdir,hxOmin1,hxOmin2,hxOmax1,&
   hxOmax2,jxOmin1,jxOmin2,jxOmax1,jxOmax2,ndir0,idirmin0
double precision :: qvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir0),&
   curlvec(ixImin1:ixImax1,ixImin2:ixImax2,idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
   ixImin2:ixImax2),surface(ixImin1:ixImax1,ixImin2:ixImax2),&
   mydx(ixImin1:ixImax1,ixImin2:ixImax2)
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop("Error in curl: Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx=1.d0/dxlevel
idirmin=4
curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixmin1:ixmax1,ixmin2:ixmax2)=qvec(ixmin1:ixmax1,ixmin2:ixmax2,kdir)
      hxOmin1=ixOmin1-kr(jdir,1);hxOmin2=ixOmin2-kr(jdir,2)
      hxOmax1=ixOmax1-kr(jdir,1);hxOmax2=ixOmax2-kr(jdir,2);
      jxOmin1=ixOmin1+kr(jdir,1);jxOmin2=ixOmin2+kr(jdir,2)
      jxOmax1=ixOmax1+kr(jdir,1);jxOmax2=ixOmax2+kr(jdir,2);
      select case(typeaxial)
        case('slab')


         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
            jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx(jdir)


        case('spherical')
         select case(jdir)
            case(1)
             tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                ixImin2:ixImax2)*mygeo%x(ixImin1:ixImax1,ixImin2:ixImax2,1)
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2))&
                /(mygeo%x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)&
                *mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
    case(2)
             mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mygeo%dx(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,2)
             if(idir==1) then
                tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                   ixImin2:ixImax2)*dsin(mygeo%x(ixImin1:ixImax1,&
                   ixImin2:ixImax2,2))
                mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsin(mygeo%x&
                   (ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*mydx(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2)
             endif
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp(jxOmin1:jxOmax1,&
                jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,hxOmin2:hxOmax2))&
                /mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
 
         end select
        case('cylindrical')
         if(-2==3) then
           select case(jdir)
              case(1)
               mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mygeo%dx(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,1)
               if(idir==3) then
                  tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*mygeo%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1)
                  mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mygeo%x&
                     (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*mydx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
               endif
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp&
                  (jxOmin1:jxOmax1,jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      case(2)
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*(tmp&
                  (jxOmin1:jxOmax1,jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  2)
   
           end select
         end if
         if(-1==3) then
           select case(jdir)
              case(1)
               mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mygeo%dx(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,1)
               if(idir==2) then
                  tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmp(ixImin1:ixImax1,&
                     ixImin2:ixImax2)*mygeo%x(ixImin1:ixImax1,ixImin2:ixImax2,&
                     1)
                  mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mygeo%x&
                     (ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*mydx(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2)
               endif
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-half*(tmp&
                  (jxOmin1:jxOmax1,jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mydx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      case(2)
               tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-half*(tmp&
                  (jxOmin1:jxOmax1,jxOmin2:jxOmax2)-tmp(hxOmin1:hxOmax1,&
                  hxOmin2:hxOmax2))/mygeo%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  2)
   
           end select
         end if
      end select
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=curlvec&
            (ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)+tmp2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)
      else
         curlvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=curlvec&
            (ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)-tmp2(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

include 'amrvacdef.f'

integer :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
double precision :: qvec(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:ndir),&
    divq(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

double precision,dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2):: qL,qR,dqC,ldq
double precision :: dxdim, invdx(1:ndim)

integer :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
   jxCmin1,jxCmin2,jxCmax1,jxCmax2,idims,ixmin1,ixmin2,ixmax1,ixmax2,gxCmin1,&
   gxCmin2,gxCmax1,gxCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,idummy
character*79, save :: savetypelimiter,savetypegradlimiter,save2typelimiter
!-----------------------------------------------------------------------------
ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;

if (ixImin1>ixmin1.or.ixImax1<ixmax1.or.ixImin2>ixmin2.or.ixImax2&
   <ixmax2) call mpistop("Error in divvectorS: Non-conforming input limits")

idummy=0
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
do idims=1,ndim
   hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
   hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);
   ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmax1=ixOmax1;ixCmax2=ixOmax2;
   jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
   jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
   gxCmin1=ixCmin1-kr(idims,1);gxCmin2=ixCmin2-kr(idims,2);gxCmax1=jxCmax1
   gxCmax2=jxCmax2;
   hxCmin1=gxCmin1+kr(idims,1);hxCmin2=gxCmin2+kr(idims,2)
   hxCmax1=gxCmax1+kr(idims,1);hxCmax2=gxCmax2+kr(idims,2);
   savetypelimiter=typelimiter
   savetypegradlimiter=typegradlimiter
   ! set the gradient limiter here
   typelimiter=savetypegradlimiter
  !! {if(gxCmin^D<ixGlo1)then
  !!   gxCmin^D=ixGlo1
  !!   hxCmin^D=hxCmin^D+1
  !! endif \}
  !! {if(hxCmax^D>ixGhi1)then
  !!   hxCmax^D=ixGhi1
  !!   gxCmax^D=gxCmax^D-1
  !! endif \}
   qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
      idims)
   qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2) = qvec(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
      idims)
   if(typelimiter/='ppm') then
      dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2)= qR(gxCmin1:gxCmax1,&
         gxCmin2:gxCmax2)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2)
      save2typelimiter=typelimiter
      if(save2typelimiter=='koren') typelimiter='korenL'
      if(save2typelimiter=='cada')  typelimiter='cadaL'
      if(save2typelimiter=='cada3') typelimiter='cada3L'
      dxdim=dxlevel(idims)
      call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
         gxCmax1,gxCmax2,idummy,idims,ldq,dxdim)
      qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) + half*ldq(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      if(save2typelimiter=='koren')then
         typelimiter='korenR'
         call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
            gxCmax1,gxCmax2,idummy,idims,ldq,dxdim)
       endif
      if(save2typelimiter=='cada')then
         typelimiter='cadaR'
         call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
            gxCmax1,gxCmax2,idummy,idims,ldq,dxdim)
       endif
      if(save2typelimiter=='cada3')then
         typelimiter='cada3R'
         call dwlimiter2(dqC,ixImin1,ixImin2,ixImax1,ixImax2,gxCmin1,gxCmin2,&
            gxCmax1,gxCmax2,idummy,idims,ldq,dxdim)
       endif
       typelimiter=save2typelimiter
      qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = qR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2) - half*ldq(jxCmin1:jxCmax1,jxCmin2:jxCmax2)
   else if (typelimiter .eq. 'ppm') then
      dqC(ixImin1:ixImax1,ixImin2:ixImax2)=qvec(ixImin1:ixImax1,&
         ixImin2:ixImax2,idims)
      call PPMlimitervar(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
         ixMhi2,idims,dqC,dqC,qL,qR)
   else
      call mpistop('typelimiter unknown in divvectorS')
   endif
   ! set the method limiter back
   typegradlimiter=typelimiter
   typelimiter=savetypelimiter

   if (slab) then
     divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+half*(qR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
        -qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2))*invdx(idims)
   else
     select case(idims)
     case(1)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC1(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC1(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
           -qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2) 
     case(2)
        qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC2(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=mygeo%surfaceC2(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+qR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)&
           -qL(hxOmin1:hxOmax1,hxOmin2:hxOmax2) 
      end select
   end if
end do
if(.not.slab) divq(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divq(ixOmin1:ixOmax1,&
   ixOmin2:ixOmax2)/mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

end subroutine divvectorS
!=============================================================================
subroutine extremaq(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,q,nshift,qMax,qMin)

include 'amrvacdef.f'

integer,intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: q(ixImin1:ixImax1,ixImin2:ixImax2)
integer,intent(in)           :: nshift

double precision, intent(out) :: qMax(ixImin1:ixImax1,ixImin2:ixImax2),&
   qMin(ixImin1:ixImax1,ixImin2:ixImax2)

integer           :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
   ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
   ishift,i,j 
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
 if (ishift==1) then
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(q(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(q(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 else
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
   qMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(qMax(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
   qMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(qMin(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2),q(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
      q(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
 end do

 
enddo

end subroutine  extremaq
!=============================================================================
subroutine extremaw(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,w,nshift,wMax,wMin)

include 'amrvacdef.f'

integer,intent(in)            :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
integer,intent(in)            :: nshift

double precision, intent(out) :: wMax(ixImin1:ixImax1,ixImin2:ixImax2,&
   1:nwflux),wMin(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux)

integer          :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
   ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
   ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
 ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
 ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
 if (ishift==1) then
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
 else
    wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= max(wMax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
    wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= min(wMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2,&
       1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,1:nwflux))
 end if
 
 idims=1
 jdims=idims+1
 do i=-1,1
   ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
   ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
   ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
   ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
   ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
   ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
   wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= &
     max(wMax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,&
        ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,&
        1:nwflux))
   wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)= &
     min(wMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux),w(ixsRmin1:ixsRmax1,&
        ixsRmin2:ixsRmax2,1:nwflux),w(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2,&
        1:nwflux))
 end do

 
enddo

end subroutine  extremaw
!=============================================================================
subroutine extremaa(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,a,nshift,aMin)

include 'amrvacdef.f'

integer,intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2
double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2)
integer,intent(in)           :: nshift

double precision, intent(out) :: aMin(ixImin1:ixImax1,ixImin2:ixImax2)

integer          :: ixsmin1,ixsmin2,ixsmax1,ixsmax2,ixsRmin1,ixsRmin2,&
   ixsRmax1,ixsRmax2,ixsLmin1,ixsLmin2,ixsLmax1,ixsLmax2,idims,jdims,kdims,&
   ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
  idims=1
  ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmin2=ixOmin2+ishift*kr(idims,2)
  ixsRmax1=ixOmax1+ishift*kr(idims,1);ixsRmax2=ixOmax2+ishift*kr(idims,2);
  ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmin2=ixOmin2-ishift*kr(idims,2)
  ixsLmax1=ixOmax1-ishift*kr(idims,1);ixsLmax2=ixOmax2-ishift*kr(idims,2);
  aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(a(ixsRmin1:ixsRmax1,&
     ixsRmin2:ixsRmax2),a(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
     a(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
  
  idims=1
  jdims=idims+1
  do i=-1,1
    ixsmin1=ixOmin1+i*ishift*kr(idims,1);ixsmin2=ixOmin2+i*ishift*kr(idims,2)
    ixsmax1=ixOmax1+i*ishift*kr(idims,1);ixsmax2=ixOmax2+i*ishift*kr(idims,2);
    ixsRmin1=ixsmin1+ishift*kr(jdims,1);ixsRmin2=ixsmin2+ishift*kr(jdims,2)
    ixsRmax1=ixsmax1+ishift*kr(jdims,1);ixsRmax2=ixsmax2+ishift*kr(jdims,2);
    ixsLmin1=ixsmin1-ishift*kr(jdims,1);ixsLmin2=ixsmin2-ishift*kr(jdims,2)
    ixsLmax1=ixsmax1-ishift*kr(jdims,1);ixsLmax2=ixsmax2-ishift*kr(jdims,2);
    aMin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(aMin(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),a(ixsRmin1:ixsRmax1,ixsRmin2:ixsRmax2),&
       a(ixsLmin1:ixsLmax1,ixsLmin2:ixsLmax2))
  end do
 
  
end do

end subroutine extremaa
!=============================================================================
subroutine locate_in_table(xpoint,table,nxtp,ixtp,res)
!  Fast search table to find xpoint's location in the table
!  record the distance between xpoint and the its cloest element table(ixtp)
! INPUT:
!  xpoint : the point you want to locate in the table
!  table : 1D table you want to find location in
!  nxtp : number of elements in the table
! OUTPUT:
!  ixtp : closest element's index in table
!  res : offset(distance) from the closest element in table
include 'amrvacdef.f'

integer, intent(in) :: nxtp
double precision,intent(in)   :: xpoint,table(nxtp)
double precision, intent(out) :: res
integer, intent(out) :: ixtp

integer :: jl,jc,jh
!-----------------------------------------------------------------------------
if(xpoint < table(1)) then
  ixtp=1
  res=(xpoint-table(1))/(table(2)-table(1))
else if (xpoint > table(nxtp)) then
  ixtp=nxtp
  res=(xpoint-table(nxtp))/(table(nxtp)-table(nxtp-1))
else
  jl=0
  jh=nxtp+1
  do
    if (jh-jl <= 1) exit
    jc=(jh+jl)/2
    if (xpoint >= table(jc)) then
        jl=jc
    else
        jh=jc
    end if
  end do
  res=(xpoint-table(jh-1))/(table(jh)-table(jh-1))
  if(res<=0.5d0) then
    ixtp=jh-1
  else
    ixtp=jh
    res=res-1.d0
  endif
end if
end subroutine locate_in_table
!=============================================================================
