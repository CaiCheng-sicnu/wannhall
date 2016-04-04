MODULE input
  !
  USE constants,   only : dp
  !
  IMPLICIT NONE
  !
  character(len=80) seed    ! 
  real(dp)    sigma
  integer     idir, jdir
  !
 CONTAINS
  !
 SUBROUTINE read_input
  !
  USE constants, only : dp, fin, Ry_to_eV
  USE para
  !
  IMPLICIT NONE
  !
  if (inode.eq.0) then
    open(unit=fin, file="wannhall.inp")
    read(fin, *) seed
    read(fin, *) sigma
    sigma=sigma/Ry_to_eV
    read(fin, *) idir, jdir
    close(fin)
  endif
  CALL para_sync(sigma)
  CALL para_sync(idir)
  CALL para_sync(jdir)
  !
 END SUBROUTINE
  !
END MODULE
