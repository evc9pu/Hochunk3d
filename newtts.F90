program tts

  ! radiative transfer in a disk+envelope around a star.
  ! height of disk, z, goes as z1*r**b.
  ! density goes as c*r**-a.
  ! opacity is dust with either rayleigh or H-G phase function.
  ! output is stokes parameters as a function of angle
  !
  ! history:
  ! 00/03/19 (mjw):  update version stamp for changes on 03/19.  Also,
  ! change version name to VGER (from mctherm)
! 2012/03/08:  add WRIMSG_MOD and initialize (serially) here
! 2012/03/10:  small tweaks, CPU output, cleanup some screen output
! 2012/03/11:  tweak random number treatment
!              include timer info on entire run 
! 2012/03/22:  add new random number generator ("remove" calls to set_random_seed and
!                 replace logic for ran() ).  added USER_SET_GENERATOR, ECUYER_COTE_MOD, and
!                 RANDOM_STANDARD_UNIFORM_MOD...all modules with a few small changes.
!                 **CURRENTLY limited to 31 processors in this implementation**
! http://mathforum.org/kb/thread.jspa?threadID=86518&messageID=421870
! https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware.aspx?Software_Id=27
!

  use tts_mod
  use grid_mod
  use output_mod
  use constants
  use wrimsg_mod
  use user_set_generator
#ifdef MPI
  use ttsre_mpi_mod
#endif

  implicit none

  integer :: debug_level,iwav, ipeelorig,nporig,iter,ncells,i1
  real*8 :: isrf_int, cpusec, cpusec_old, cpusec_start,delta_cpusec
  character(len=80) :: atname,dustname(ndg),cfllog,cmsgnm
  character(len=3) :: converge

#ifdef MPI

! ... parallel version stamp and processor info
  character :: CPAR_STAMP*40,CPROC*40
  integer :: i,len_proc
  data CPAR_STAMP/'PARALLEL VERSION 1.1a [2012.03.22]'/

! ... misc scalars
  integer :: np_reset

! define a few useful things
  call ttsre_mpi_init()
  call initialize_wrimsg_parallel()
!      call allocate_mars_mpi()

! ... Initialize MPI for this process, getting MPI_COMM_WORLD process
!     ID and total number of processes.
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )


!
! ... create file tag from MYID (for child processes) and open
!     log file for output of IDVOUT and IDVERR streams
!
  if (myid .eq. 0) then
! ****************************************************************************
! ****************************** Parent PROCESS  *****************************
! ****************************************************************************
! ... send output to screen; user can redirect if desired
     write(cfllog,'(a)') 'p_disk-parent.log'
! ... create new communicator for child processes, so pass
!     undefined color for parent
     call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,&
          myid-1,COMM_REDUCE,ierr)
  else
! ****************************************************************************
! ******************************* Child PROCESSES  ***************************
! ****************************************************************************

! ... Create file tag for each process, using cascade logic (i.e., don't
!     need to test ".ge." case as it should have triggered previous "if")
     if (myid .lt. 10000) then
        write(ext,'(i4.4)') myid
     else
        stop 'HA HA HA...increase logic block of CFLLOG'
     endif
     cfllog = 'p_disk-child'//ext//'.log'

! ... Create new communicator for child processes, defining rank
!     to be in the same order as MPI_COMM_WOLRD (i.e., myid=1
!     will become myid_reduce=0)
     call MPI_COMM_SPLIT(MPI_COMM_WORLD,0,myid-1,COMM_REDUCE,&
          ierr)
  endif

  open(unit=IDVOUT,file=cfllog,status='unknown')
  write(IDVOUT,'(a,a)') version,' - ttsre version'
  write(IDVOUT,'(a,i4,a,i4,a)') &
       '(MPI_COMM_WORLD): Process number ', myid, ' of ',&
       numprocs,' total (proc. no. start with 0)'

! ... trap too many processors for current RNG
  if (numprocs > 31) then
     write(0,*) 'CANNOT USE MORE THAN 31 PROCESSES AT THIS TIME.'
     goto 0990
  endif

! COMM_REDUCE is artifact...will need to remove eventually
  if (myid > 0) then
     call MPI_COMM_RANK( COMM_REDUCE, myid_reduce, ierr )
     call MPI_COMM_SIZE( COMM_REDUCE, numprocs_reduce, ierr )
     write(IDVOUT,'(a,i4,a,i4,a)') &
          '(COMM_REDUCE): Process number ', myid_reduce, ' of ',&
          numprocs_reduce,' total (proc. no. start with 0)'
  endif

! ... get processor info and writeout version stamp
  call MPI_GET_PROCESSOR_NAME(CPROC, len_proc, ierr)
  CPROC =  'Processor = '//CPROC(1:28)
  write(IDVOUT,'(a)') CPROC
  write(IDVOUT,'(a)') CPAR_STAMP

  
#else

  call initialize_wrimsg_serial()
  open(unit=12,file='disk.dat',status='unknown')
  write(12,'(a,a)') version,' - ttsre version'

#endif

  debug_level = 0






  ! input and setup
  call setup(atname,dustname,i1)

#ifdef MPI
  
  write(cmsgnm,'(a,i4)') 'Setting Random Number Generator for stream ',myid+1
  call user_set_all(12345678,87654321,myid+1)

#else

  write(cmsgnm,*) 'Setting Random Number Generator for stream 1'
  call user_set_all(12345678,87654321,1)

#endif
  call wrimsg('TTS',cmsgnm)

  ipeelorig=ipeel
  nporig=np
  if (ilucy.eq.1) then

     ipeel=0

     ! for all but last iteration only need as many photons to get
     ! good temperature.  determined in chakrabarti & whitney 07
     ! that np_temp ~ number of grid cells.
     ! user inputs np for final iteration to get desired SED fidelity

     ncells = (nrg-1)*(ntg-1)*(npg-1)
     if (ntg.lt.10.and.npg.lt.10) then
        np=100*ncells
     else
        np=3.*ncells
     end if

   !  if (np.lt.nporig/10) np=nporig/10
     if (np.gt.nporig) np=nporig
     if (np.lt.npmin) np = npmin

     WRITE(0,*) 'newtts_lucy, np per iteration, np_final',np,nporig
     write(12,*) 'newtts_lucy, np per iteration, np_final',np,nporig

  end if

#ifdef MPI
  write(cmsgnm,'(a,i9)') 'Reducing NP so that NP_new * num_processors = NP_old =  ', np
  call wrimsg('TTS',cmsgnm)
  np = np / numprocs

#endif


  npsav=np
  write(IDVOUT,*) 'np,nporig',np,nporig
  iter=0

  write(IDVOUT,*) 'calling dustinit'
  call dustinit(dustname)
  write(IDVOUT,*) 'done with dustinit, calling draine_cdf'

  call draine_cdf() ! small grains
  call isrf_cdf(isrf_int,isrf_Av)   !need this for energy density normalization
  write(IDVOUT,*) 'tts, done with draine_cdf'
  write(IDVOUT,*) 'tts, done with isrf_cdf'
  if (iout.eq.1) then
     ! interstellar radiation field (external illumation)
   !  call isrf_cdf(isrf_int,isrf_Av)
     l_isrf=isrf_int*4.d0*pi*(rmax*rsol*rstar)**2/4.d0*isrf_scl
     print*,'isrtf_int,rmax*rsol*rstar,isrf_scl'
     print*,isrf_int,rmax*rsol*rstar,isrf_scl

     ! I think you divide by 2 because at each point on the sphere, only
     ! half the radiation is headed into the sphere.  I think....
     ! if (isrf_frac.gt.1.d0) then
     ! print*,'you should not be here unless roundoff error'
     ! print*,'isrf_frac ',isrf_frac
     ! isrf_frac=1.d0
     ! print*,'all luminosity is in external radiation, now = ',
     ! $           l_isrf
     ! end if

     ! testing
     ! if doing external illumination only
     ! isrf_frac=1.0
     ! lstars=0.d0
     print*,'l_isrf (lsun)',l_isrf/lsun
  else
     l_isrf=0.d0
     ! isrf_frac=0.d0
  end if

  call atmosinit(atname)
  write(IDVOUT,*) 'done with atmosinit'

  iwav=1

#ifdef MPI
! ... Although "race" conditions shouldn't be occurring, let's wait
!     for everyone here.
  write(IDVOUT,*) 'Random number seed barrier...'
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if (myid .eq. 0) then
! ****************************************************************************
! ****************************** PARENT PROCESS  *****************************
! ****************************************************************************
     write(IDVOUT,'(a,i8)') 'PARENT>> base random number seed = ',&
          i1

     do i=1,numprocs-1
        call MPI_RECV(iready,1,MPI_INTEGER,MPI_ANY_SOURCE,&
             tag_ready,MPI_COMM_WORLD, status, ierr)
        sender = status(MPI_SOURCE)
        ibuf(1) = i1*((sender+1)**3)
        write(cmsgnm,'(a,i11,a,i4)') &
             'PARENT>> Sending random number seed (i1) = ',&
             ibuf(1),' to process ',sender
        call wrimsg('TTS',cmsgnm)
        call MPI_SEND(ibuf,4,MPI_INTEGER,sender,&
             tag_start,MPI_COMM_WORLD,ierr)
     end do

  else
! ****************************************************************************
! ******************************* CHILD PROCESSES  ***************************
! ****************************************************************************
     call MPI_SEND(iready,1,MPI_INTEGER,0,tag_ready,&
          MPI_COMM_WORLD, ierr)
     write(IDVOUT,'(a)') 'CHILD>> Sending IREADY signal to Parent'
     call MPI_RECV(ibuf,4,MPI_INTEGER,0,&
          tag_start,MPI_COMM_WORLD, status, ierr)
     i1 = -abs(ibuf(1))
     write(cmsgnm,'(a,i11)') &
          'CHILD>> Receiving random number seed (i1) from parent: ',i1
     call wrimsg('TTS',cmsgnm)
  endif
#endif


  call cpu_time(cpusec_start)

  ! start iteration here, since setup_wave uses inner radius which
  ! will change as dust destruction radius gets updated

  call initarr(.true.)
  call setup_wave()


  if (ilucy.eq.0) then
     converge='yes'
  else
     converge='no'
  end if

  do while (converge.eq.'no')

     iter=iter+1

     call cpu_time(cpusec)
     write(cmsgnm,'(a,i2,a,f11.2)') 'Starting Lucy iteration ',iter,' at CPU time ',cpusec
     call wrimsg('TTS',cmsgnm)

! reset np using npsav, this seems to only be necessary for MPI usage, but since npsav was in the serial 
! code, we'll do it in both versions
     np = npsav
     write(cmsgnm,'(a,i9)') 'Reseting NP in Lucy Iteration Loop to ',npsav
     call wrimsg('TTS',cmsgnm)

     write(12,*) ''
     write(12,*) ' iteration # ',iter

     ! radiative transfer
     call disk(converge,iter,i1)

     cpusec_old = cpusec
     call cpu_time(cpusec)
     write(cmsgnm,'(a,i2,a,f11.2,a)') 'Elapsed time for Lucy iteration ',iter,': ',&
          cpusec-cpusec_old,' (sec CPU time)'
     call wrimsg('TTS',cmsgnm)

     ! if number of iterations is not sufficient, force converge to 'no'
     if(iter.lt.n_iter_min.and.ilucy.eq.1) converge='no'
     if(iter.ge.n_iter_max.and.ilucy.eq.1) then
        write(cmsgnm,'(a)') 'at user-specified max iteration--CAUTION, T may not have converged!'
        call wrimsg('TTS',cmsgnm)
        converge='yes'
      end if

     call initarr(.false.)

  end do

  ! final iteration (or first if ilucy=0)
  ipeel=ipeelorig
  np=nporig
#ifdef MPI
  write(cmsgnm,'(a,i9)') 'Reducing NP so that NP_new * num_processors = NP_old =  ', np
  call wrimsg('TTS',cmsgnm)
  np = np / numprocs

#endif
  iter=iter+1

  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2)') 'Initializing images and calling DISK at CPU time ',cpusec
  call wrimsg('TTS',cmsgnm)

  call initimages()

  call disk(converge,iter,i1)

  cpusec_old = cpusec
  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2,a)') 'Elapsed time for final RT call to DISK',cpusec-cpusec_old,' (sec CPU time)'
  call wrimsg('TTS',cmsgnm)


#ifdef MPI
! only parent process needs to continue, all others go to barrier and wait
  if (myid > 0) goto 0990
#endif

  ! convert stokes arrays to flux arrays
  call diskflux()

! output images/spectra
  write(cmsgnm,'(a)') 'Done with diskflux, calling output.'
  call wrimsg('TTS',cmsgnm)
  call output()

  write(cmsgnm,'(a,i2)') 'total number of iterations: ',iter
  call wrimsg('TTS',cmsgnm)

! final timing call
  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2,a)') 'Elapsed CPU time for execution of main program: ',&
       cpusec-cpusec_start,' (sec)'
  call wrimsg('TTS',cmsgnm)



#ifdef MPI
! ... Clean up MPI baggage for this process
0990 call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(irc)
  if (myid.eq.0) then
     WRITE(IDVOUT,'(a)') 'Parent process terminating....'
  else
     WRITE(IDVOUT,'(a)') 'Child process terminating....'
  endif
  WRITE(IDVOUT,'(a)') 'Have a nice day.'
  close(IDVOUT)
#else

  close(12)

#endif


end program tts
