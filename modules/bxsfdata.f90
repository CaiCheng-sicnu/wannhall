MODULE bxsfdata
  !
  USE constants,   only : dp
  !
  IMPLICIT NONE
  !
  integer nkx, nky, nkz, nbnd
  real(dp), allocatable :: kvec(:, :)
  real(dp), allocatable :: ek(:, :)
  real(dp), dimension(1:3, 1:3) :: acell
  real(dp) ef, omega
  !
 CONTAINS
  !
 SUBROUTINE read_bxsf(fn)
  !
  USE constants,  only : dp, fin, stdout, bohr_to_angstrom, Ry_to_eV, twopi
  USE para
  !
  IMPLICIT NONE
  !
  character(len=80) fn
  integer ikx, iky, ikz, ibnd, ii, nkpt
  real(dp), dimension(1:3, 1:3) :: bvec
  !
  if (inode.eq.0) then
    !
    write(stdout, *) " # Reading calculated mesh data from "//trim(fn)
    !
    open(unit=fin, file=trim(fn))
    !
    do ii=1, 8
      read(fin, *)
    enddo
    !
    read(fin, '(20X,1F)') ef
    !
    write(stdout, '(1A,1F16.9,1A)', advance='no') " # Fermi level: ", ef, "eV, "
    !
    ef=ef/Ry_to_eV
    !
    write(stdout, '(1F16.9,1A)') ef, "Ry"
    !
    do ii=1, 5
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
    read(fin, *)
    do ii=1, 3
      read(fin, *) bvec(:, ii)
    enddo
    !
    bvec(:, :)=bvec(:, :)*bohr_to_angstrom/twopi
    !
    CALL cross_product(acell(:, 1), bvec(:, 2), bvec(:, 3))
    CALL cross_product(acell(:, 2), bvec(:, 3), bvec(:, 1))
    CALL cross_product(acell(:, 3), bvec(:, 1), bvec(:, 2))
    omega=1.d0/sum(acell(:, 1)*bvec(:, 1))
    !
    acell(:, :)=acell(:, :)*omega
    !
    write(stdout, *) " # Original unit cell: (in bohr)"
    do ii=1, 3
      write(stdout, '(3F22.16)') acell(:, ii)
    enddo
    !
    write(stdout, '(1A,1F22.16)') " # unit cell volume: ", omega
    !
    do ibnd=1, nbnd
      !
      read(fin, *) fn
      !
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
      !
    enddo  ! ibnd
    !
    ek(:, :)=ek(:, :)/Ry_to_eV
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
 SUBROUTINE cross_product(v1, v2, v3)
  !
  USE constants, only : dp
  !
  IMPLICIT NONE
  !
  real(dp), dimension(1:3) :: v1, v2, v3
  !
  v1(1)=v2(2)*v3(3)-v2(3)*v3(2)
  v1(2)=v2(3)*v3(1)-v2(1)*v3(3)
  v1(3)=v2(1)*v3(2)-v2(2)*v3(1)
  !
 END SUBROUTINE
  !
END MODULE
