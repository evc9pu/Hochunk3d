module grid_mod

  implicit none
  save

  ! *****  set ntg to an odd number!!!!! *******
  ! will work on a fix later...

  ! Grid dimensions
  integer :: nrg, ntg, npg

  ! number of opacity species (dust/gas/sgs)
  integer, parameter :: ndg = 8
  logical,parameter :: is_any(ndg) = [1,1,1,1,1,1,1,1] == 1
  logical,parameter :: is_disk(ndg) = [1,1,0,0,1,1,0,0] == 1
  logical,parameter :: is_envelope(ndg) = [0,0,1,0,0,0,1,0] == 1
  logical,parameter :: is_cavity(ndg) = [0,0,0,1,0,0,0,1] == 1
  logical,parameter :: is_sg(ndg) = [0,0,0,0,1,1,1,1] == 1
  logical,parameter :: is_thermal(ndg) = .not. is_sg
  logical,parameter :: is_thermal_disk(ndg) = [1,1,0,0,0,0,0,0] == 1
  
  integer :: fractmask(ndg)

  real*8,parameter :: texp = 1.5d0
  real*8,parameter :: rexp = 3.0d0
  ! real*8,parameter :: rexp = 2.0   ! constant density 1-D geometry

  integer,allocatable :: killed(:,:,:)

  real*8,allocatable :: densarr(:,:,:,:)
  real*8,allocatable :: nabs(:,:,:,:)
  real*8,allocatable :: densvararr(:,:,:)
  real*8,allocatable :: massarr(:,:,:,:)
  real*8,allocatable :: tdust(:,:,:,:)
  real*8,allocatable :: dtauabs(:,:,:,:)
  real*8,allocatable :: tdust2(:,:,:,:)
  real*8,allocatable :: tdust2dw(:,:,:)
  real*8,allocatable :: jmean(:,:,:),imfract(:,:,:)
  real*8,allocatable :: vcell(:,:,:)
  real*8,allocatable :: rarr(:),r2arr(:)
  real*8,allocatable :: sigarr(:,:,:)
  real*8,allocatable :: thetarr(:),tmptharr(:)
  real*8,allocatable :: sintarr(:),costarr(:)
  real*8,allocatable :: tan2arr(:),phiarr(:)
  real*8,allocatable :: aarr(:),barr(:),tauarr(:)
  real*8,allocatable :: ravearr(:),thetavearr(:),phiavearr(:)

  integer,allocatable :: diffdir(:,:,:)
  integer,allocatable :: imean(:,:,:),dbig(:,:,:),dsg(:,:,:)
  integer :: iwarn
  logical,allocatable :: diffus(:,:,:),nodust(:,:,:,:)

contains

  subroutine grid_allocate()

    implicit none

    write(*,'("Setting grid dimensions to (",I0,",",I0,",",I0,")")') nrg, ntg, npg

    if(nrg*ntg*npg > 1e8) stop "Number of grid cells requested exceeds 100 million"

    if(nrg==0.or.ntg==0.or.npg==0) then
       stop "ERROR: grid dimensions not set correctly"
    end if

    allocate(killed(nrg, ntg, npg))

    ! Physical quantities
    allocate(densarr(nrg,ntg,npg,ndg))
    allocate(densvararr(nrg,ntg,npg))
    allocate(massarr(nrg,ntg,npg,ndg))
    allocate(tdust(nrg,ntg,npg,ndg))
    allocate(dtauabs(nrg,ntg,npg,ndg))
    allocate(tdust2(nrg,ntg,npg,ndg))
    allocate(tdust2dw(nrg,ntg,npg))
    allocate(jmean(nrg,ntg,npg),imfract(nrg,ntg,npg))
    allocate(nabs(nrg,ntg,npg,ndg),diffdir(nrg,ntg,npg))
    allocate(imean(nrg,ntg,npg),dbig(nrg,ntg,npg),dsg(nrg,ntg,npg))
    allocate(diffus(nrg,ntg,npg),nodust(nrg,ntg,npg,ndg))

    ! Geometrical quantities
    allocate(vcell(nrg,ntg,npg))
    allocate(rarr(nrg),r2arr(nrg))
    allocate(sigarr(nrg,npg,ndg))
    allocate(thetarr(ntg),tmptharr(ntg))
    allocate(sintarr(ntg),costarr(ntg))
    allocate(tan2arr(ntg),phiarr(npg))
    allocate(aarr(npg),barr(npg),tauarr(nrg))
    allocate(ravearr(nrg-1),thetavearr(ntg-1),phiavearr(npg-1))

  end subroutine grid_allocate

end module grid_mod
