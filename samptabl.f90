subroutine table

  ! makes a table of integrated probability distribution function.
  ! for limb darkening

  use tabl_mod
  implicit none

  integer :: i

  delmu=1.0d0/dble(ntab)
  xmu(1)=0.0d0

  do i=2,ntab
     xmu(i)=delmu*dble(i)
     prob(i)=0.5d0*(xmu(i)**3.d0 + xmu(i)**2.d0)
  end do

  write(6,*) 'in table, prob(n), should equal 1, ',prob(ntab)
  print*,'ntab',ntab

  return
end subroutine table

! *********************************************************

! *********************************************************

subroutine darkening(random,mu)

  ! interpolates from table (subroutine table) to sample
  ! mu from probability distribution function.
  ! mu ranges from 0 to 1

  use tabl_mod
  implicit none

  integer :: n
  real*8 :: random,mu

  ! do i=1,ntab
  ! check=random-prob(i)
  ! if (check.lt.0.0) then
  ! n=i
  ! go to 10
  ! end if
  ! end do

  ! locate index in prob array
  call locate(prob,ntab,random,n)

  ! write(6,*) 'table for mu must be wrong...'

  do
     ! mu=xmu(n)-(prob(n)-random)/(prob(n)-prob(n-1))*delmu
     mu=xmu(n+1)-(prob(n+1)-random)/(prob(n+1)-prob(n))*delmu
     if (mu.gt.1.0d0.or.mu.lt.0.d0) then
        print*,'dammit'
     else
        exit
     end if
  end do

  return
end subroutine darkening

