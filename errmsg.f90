subroutine errmsg(cstats,csubrt,cmsgnm)
  implicit none
  character(len=*),intent(in) :: cmsgnm,cstats,csubrt
  if(cstats=='FATAL') then
     write(*,'(" >>>>> FATAL ERROR IN PROCEDURE: ",a)') csubrt
     write(*,'(1X,A)') cmsgnm
     write(*,'(" >>>>> EXECUTION ABORTED ")')
     stop
  else if(cstats=='WARNING') then
     write(*,'(" >>>>> WARNING IN PROCEDURE: ",a)') csubrt
     write(*,'(1X,A)') cmsgnm
  else
     write(*,'("ERROR: unknown cstats in errmsg: ",A)') cstats
     stop
  end if
end subroutine errmsg
