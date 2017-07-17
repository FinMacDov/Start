INCLUDE:vacusr.gravity.t
INCLUDE:vacusr.viscosity.t

!=============================================================================
subroutine specialini(ix^L,w)

include 'vacdef.f'

integer:: ix^L
double precision:: w(ixG^T,1:nw)

double precision:: rhoin,xcent^D,radius
double precision:: inirho,iniene
double precision:: onemor,inix,ddx
double precision:: p_1,p_2

integer:: iii_,iix_1,info,i,j
double precision:: pi,comi,eneu,sum,mode,bmax,l
character*79 atmfilename

integer:: ix_1,ix_2

double precision:: ixc_1,ixc_2
double precision:: rfc,a1,a2


double precision:: bz,xc, bmin,bmax

!-----------------------------------------------------------------------------

{^IFMPI if (ipe .ne. 0) read(unitpar,'(a)') atmfilename}

{^IFMPI if (ipe .eq. 0) then}
write(*,*)'param load'
read(unitpar,'(a)') atmfilename
write(*,*) 'from file ',atmfilename
open(33,file=atmfilename,status='old')

iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)

{^IFMPI endif}

{^IFMPI call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(iniene,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}

do ix_1=ixGhi1,ixGlo1,-1
{^IFMPI if (ipe .eq. 0)} read(33,*) inix,inirho

{^IFMPI if (ipe .eq. 0)} print*,'inix,inirho=',inix,inirho


{^IFMPI call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(inix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}
{^IFMPI call MPI_BCAST(inirho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierrmpi)}

!reading in data files and seting values for conservatives variables. 
 do ix_2=ixGlo2,ixGhi2

   x(ix_1,ix_2,1)=inix*1000.d0
   w(ix_1,ix_2,rho_)=inirho
   w(ix_1,ix_2,e_)=iniene !This is energy at the highest part of the domain.
   w(ix_1,ix_2,m1_)=0.0
   w(ix_1,ix_2,m2_)=0.0

 enddo

!Define your energy and calc the pressure profile working down the domain. It would be better to do this AMRVAC as you call the function to convert from prim to conv and also call the numerical shceme which you have chosen.

enddo
{^IFMPI if (ipe .eq. 0)} close(33)

{^IFMPI if (ipe .eq. 0) } print*,'grav=',eqpar(grav0_),eqpar(grav1_),eqpar(grav2_)



print*, '###################',ix^L

call primitive(ix^L,w)

do ix_2=ixGlo2,ixGhi2
 do ix_1=ixGhi1-1,ixGlo1,-1 
!I dont understand the point in comi
comi=-abs(x(ix_1+1,ix_2,1)-x(ix_1,ix_2,1))

! I dont understand this line. Units wise it works you are add pressure to the system
w(ix_1,ix_2,p_)=w(ix_1+1,ix_2,p_)+w(ix_1,ix_2,rho_)*comi*1.d0*eqpar(grav1_)

if (ix_2 .eq. ixGlo2) print*,'eee=',w(ix_1,ix_2,rho_),w(ix_1,ix_2,p_)

 enddo
enddo


!goto 200

do ix_2=ixGlo2,ixGhi2
 do ix_1=ixGlo1+2,ixGhi1-2
       
       w(ix_1,ix_2,rho_)=-(1.D0/eqpar(grav1_))*(1.D0/(12.D0*(x(ix_1+1,ix_2,1)-x(ix_1,ix_2,1))))*(w(ix_1+2,ix_2,p_)-8.D0*w(ix_1+1,ix_2,p_)+8.D0*w(ix_1-1,ix_2,p_)-w(ix_1-2,ix_2,p_))
               
if (ix_2 .eq. ixGlo2) print*,'eeee=',w(ix_1,ix_2,rho_),w(ix_1,ix_2,p_)

     enddo
   enddo

!lower boundary
do ix_1=ixmin1+4,ixmin1+2,-1
  do ix_2=ixmin2,ixmax2
!    do ix_3=ixmin3,ixmax3
!        p_1=w(ix_1+2,ix_2,ix_3,p_)-8.d0*w(ix_1+1,ix_2,ix_3,p_)+8.d0*w(ix_1-1,ix_2,ix_3,p_)
!        p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
!        w(ix_1-2,ix_2,ix_3,p_) = 12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1
         p_1=w(ix_1+2,ix_2,p_)-8.d0*w(ix_1+1,ix_2,p_)+8.d0*w(ix_1-1,ix_2,p_)
         p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
         w(ix_1-2,ix_2,p_) = 12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1
!       enddo
    enddo
 enddo


!upper boundary
do ix_1=ixmax1-4,ixmax1-2
   do ix_2=ixmin2,ixmax2
!      do ix_3=ixmin3,ixmax3
         
!          p_1=w(ix_1-2,ix_2,ix_3,p_)-8.d0*w(ix_1-1,ix_2,ix_3,p_)+8.d0*w(ix_1+1,ix_2,ix_3,p_)
!          p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
!          w(ix_1+2,ix_2,ix_3,p_) = -12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1-1,ix_2,ix_3,1))*p_2+p_1
           p_1=w(ix_1-2,ix_2,p_)-8.d0*w(ix_1-1,ix_2,p_)+8.d0*w(ix_1+1,ix_2,p_)
           p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
           w(ix_1+2,ix_2,p_) = -12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1
      enddo
   enddo
!enddo


200 continue



call conserve(ix^L,w)




w(ix^S,eb_)=w(ix^S,e_)
w(ix^S,e_)=0.d0

w(ix^S,rhob_)=w(ix^S,rho_)
w(ix^S,rho_)=0.d0
