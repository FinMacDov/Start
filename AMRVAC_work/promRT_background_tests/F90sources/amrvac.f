!> AMRVAC solves a set of hyperbolic equations
!> \f$\vec{u}_t + \nabla_x \cdot \vec{f}(\vec{u}) = \vec{s}\f$
!> using adaptive mesh refinement.
program amrvac


! following line may avoid mystery problems in the combination
!  where we use ifort compiler and MPT (SGIs MPI implementation)
!DEC$ ATTRIBUTES NOINLINE :: read_snapshot


include 'amrvacdef.f'

integer :: itin
double precision :: time0, time_in, tin
!-----------------------------------------------------------------------------
call comm_start

time_advance=.false.

time0=MPI_WTIME()
time_bc=zero


! get input parameters
call readcommandline
call readparameters
call initialize_vars
call init_comm_types


! Begin of the code
! -----------------
if (snapshotini/=-1) then
   ! restart from previous file or dat file conversion
   ! get input data from previous VAC/AMRVAC run
   itin=it
   tin=t

   ! order: first physics dependent part then user for defaults.
   call initglobal



   if (typeparIO==1)then
     call read_snapshot
   else
     call read_snapshotnopar
   end if



   ! modify globals
   if (changeglobals) call initglobal



   if (convert) then
       if (npe/=1.and.(.not.(index(convert_type,'mpi')>=1)) &
          .and. convert_type .ne. 'user')  call mpistop&
          ("non-mpi conversion only uses 1 cpu")
       call generate_plotfile

       call comm_finalize
       stop
   end if

   if (itreset) it=itin
   if (treset) t=tin
   ! modify initial condition
   if (firstprocess) call modify_IC
   ! reset AMR grid
   if (resetgrid) call settree
   ! Select active grids
   call selectgrids

else
   ! order: first physics dependent part then user for defaults.
   call initglobal



   ! form and initialize all grids at level one
   call initlevelone
   ! set up and initialize finer level grids, if needed
   call settree



end if

! set up boundary flux conservation arrays
if (levmax>levmin) call allocateBflux

if (mype==0) then
   print*,'-------------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Startup phase took : ',MPI_WTIME()-time0,' sec'
   print*,'-------------------------------------------------------------------------------'
end if



time_advance=.true.

! do time integration of all grids on all levels
call timeintegration

if (mype==0) then
   print*,'-------------------------------------------------------------------------------'
   write(*,'(a,f17.3,a)')' Finished AMRVAC in : ',MPI_WTIME()-time0,' sec'
   print*,'-------------------------------------------------------------------------------'
end if


call comm_finalize

end program amrvac
!=============================================================================
subroutine timeintegration

use mod_timing
include 'amrvacdef.f'

integer :: level, ifile, fixcount

logical, external :: timetosave, fixgrid

!-----------------------------------------------------------------------------
time_in=MPI_WTIME()
fixcount=1

n_saves(filelog_:fileout_) = snapshotini

do ifile=nfile,1,-1
   if (time_accurate) tsavelast(ifile)=t
   itsavelast(ifile)=it
end do

itmin=it
! the next two are used to keep track of the performance during runtime:
itTimeLast=it
timeLast=MPI_WTIME()

call getbc(t,0.d0,ixGlo1,ixGlo2,ixGhi1,ixGhi2,pw,pwCoarse,pgeo,pgeoCoarse,&
   .false.,0,nwflux+nwaux)

!  ------ start of integration loop. ------------------
timeloop0=MPI_WTIME()
timefirstbc=timeloop0-time_in
time_bc=0.d0
if (mype==0) then
   write(*,'(a,f12.3,a)')' BCs before Advance took : ',timefirstbc,' sec'
end if

time_evol : do

   call setdt
   if(fixprocess) call process(it,t)

   timeio0=MPI_WTIME()
   do ifile=nfile,1,-1
      if(timetosave(ifile)) call saveamrfile(ifile)
   end do

   timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

   ! exit time loop criteria
   if (it>=itmax) exit time_evol
   if (time_accurate .and. t>=tmax) exit time_evol

   call advance(it)

   if((.not.time_accurate).or.(residmin>smalldouble)) then
      call getresidual(it)
   endif 
   if (residual<residmin .or. residual>residmax) exit time_evol

   ! resetting of tree BEFORE IO and setdt
   timegr0=MPI_WTIME()
   if(ditregrid>1) then
     if(fixcount<ditregrid) then
       fixcount=fixcount+1
     else
       if (mxnest>1 .and. .not.(fixgrid(0))) call resettree
       fixcount=1
     endif
   else
     if (mxnest>1 .and. .not.(fixgrid(0))) call resettree
   endif
   timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

   it = it + 1
   if (time_accurate) t = t + dt
   if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
   
   if(it>9000000)then
     it = slowsteps+10
     itmin=0
     itsavelast(:)=0
   end if
   
end do time_evol

timeloop=MPI_WTIME()-timeloop0

if (mype==0) then
   write(*,'(a,f12.3,a)')' Total timeloop took        : ',timeloop,' sec'
   write(*,'(a,f12.3,a)')' Time spent on Regrid+Update: ',timegr_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timegr_tot&
      /timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on IO in loop   : ',timeio_tot,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timeio_tot&
      /timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on BC           : ',time_bc,' sec'
   write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*time_bc&
      /timeloop,' %'
   write(*,'(a,f12.3,a)')' Time spent on run          : ',timeloop&
      -timeio_tot,' sec'
end if

timeio0=MPI_WTIME()
do ifile=nfile,1,-1
   if(itsavelast(ifile)<it)call saveamrfile(ifile)
enddo
if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)




if (mype==0) then
   write(*,'(a,f12.3,a)')' Total time spent on IO     : ',timeio_tot,' sec'
   write(*,'(a,f12.3,a)')' Total timeintegration took : ',MPI_WTIME()&
      -time_in,' sec'
end if

end subroutine timeintegration
!=============================================================================
logical function timetosave(ifile)

! Save times are defined by either tsave(isavet(ifile),ifile) or
! itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
! Other conditions may be included.

include 'amrvacdef.f'

integer:: ifile
logical:: oksave
!-----------------------------------------------------------------------------
oksave=.false.
if (it==itsave(isaveit(ifile),ifile)) then
   oksave=.true.
   isaveit(ifile)=isaveit(ifile)+1
end if
if (it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.
if (time_accurate) then
   if (t>=tsave(isavet(ifile),ifile)) then
      oksave=.true.
      isavet(ifile)=isavet(ifile)+1
   end if
   if (t>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
     oksave=.true.
     n_saves(ifile) = n_saves(ifile) + 1
   endif
end if
if (oksave) then
   if (time_accurate) tsavelast(ifile) =t
   itsavelast(ifile)=it
end if
timetosave=oksave

return
end function timetosave

!> Return true if the AMR grid should not be adapted any more. This is
!> controlled by tfixgrid or itfixgrid. Other conditions may be included.
!> @todo Fix dummy argument?
logical function fixgrid(dummy)

  include 'amrvacdef.f'
  integer :: dummy              !< Unused dummy argument

  fixgrid= (t>=tfixgrid .or. it>=itfixgrid)

end function fixgrid

!> Perform global initialization, in three steps:
!> * Call initglobaldata from the selected physics module
!> * Call initglobaldata_usr which is supplied by the user
!> * Call checkglobaldata from the selected physics module
subroutine initglobal

  include 'amrvacdef.f'

  call initglobaldata
  call initglobaldata_usr
  call checkglobaldata

end subroutine initglobal
