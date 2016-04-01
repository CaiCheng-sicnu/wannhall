MODULE input
  !
  USE constants,   only : dp
  !
  IMPLICIT NONE
  !
  character(len=80) seed    ! 
  real(dp)    sigma, ef
  integer     idir, jdir
  !
 CONTAINS
  !
 SUBROUTINE read_input
  !
  USE constants, only : dp, fin
  USE para
  !
  IMPLICIT NONE
  !
  if (inode.eq.0) then
    open(unit=fin, file="wannhall.inp")
    read(fin, *) seed
    read(fin, *) ef, sigma
    read(fin, *) idir, jdir
    close(fin)
  endif
  CALL para_sync(ef)
  CALL para_sync(sigma)
  CALL para_sync(idir)
  CALL para_sync(jdir)
  !
 END SUBROUTINE
  !
END MODULE
