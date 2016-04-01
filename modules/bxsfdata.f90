MODULE bxsfdata
  !
  USE constants,   only : dp
  !
  IMPLICIT NONE
  !
  integer nkx, nky, nkz, nbnd
  real(dp), allocatable :: kvec(:, :)
  real(dp), allocatable :: ek(:, :)
  real(dp) ef
  !
 CONTAINS
  !
 SUBROUTINE read_bxsf(fn)
  !
  USE constants,  only : dp, fin, stdout
  USE para
  !
  IMPLICIT NONE
  !
  character(len=80) fn
  integer ikx, iky, ikz, ibnd, ii, nkpt
  !
  if (inode.eq.0) then
    !
    write(stdout, *) " # Reading calculated mesh data from "//trim(fn)
    !
    open(unit=fin, file=trim(fn))
    !
    do ii=1, 10
      read(fin, *)
    enddo
    !
    read(fin, *) 
    !
    do ii=1, 3
      read(fin, *)
    enddo
    !
    read(fin, *) nbnd
    !
    read(fin, *) nkx, nky, nkz
    !
    nkx=nkx-1
    nky=nky-1
    nkz=nkz-1
    !
    nkpt=nkx*nky*nkz
    !
  endif
  !
  CALL para_sync(nkpt)
  CALL para_sync(nbnd)
  !
  allocate(kvec(1:3, 1:nkpt))
  allocate(ek(1:nbnd, 1:nkpt))
  !
  if (inode.eq.0) then
    do ii=1, 4
      read(fin, *)
    enddo
    !
    do ibnd=1, nbnd
      read(fin, *)
      do ikx=0, nkx
      do iky=0, nky
      do ikz=0, nkz
        if ((ikx<nkx).and.(iky<nky).and.(ikz<nkz)) then
          ii=ikx*nky*nkz+iky*nkz+ikz+1
          kvec(1, ii)=ikx*1.d0/nkx
          kvec(2, ii)=iky*1.d0/nky
          kvec(3, ii)=ikz*1.d0/nkz
          read(fin, *) ek(ibnd, ii)
        else
          read(fin, *)
        endif
      enddo !ikz
      enddo !iky
      enddo !ikx
    enddo  ! ibnd
    !
    close(unit=fin)
    !
    write(stdout, *) " # done."
    !
  endif
  !
  CALL para_sync(kvec, 3, nkpt)
  CALL para_sync(ek, nbnd, nkpt)
  !
 END SUBROUTINE
  !
 SUBROUTINE finalize_bxsf()
  !
  IMPLICIT NONE
  !
  if (allocated(kvec)) deallocate(kvec)
  if (allocated(ek)) deallocate(ek)
  !
 END SUBROUTINE
  !
END MODULE
