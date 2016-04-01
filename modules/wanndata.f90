include 'lapack.f90'

!
!   wanndata.f90
!   
!
!   Created by Chao Cao on 01/03/14.
!   Copyright 2013 __MyCompanyName__. All rights reserved.
!

MODULE wanndata
  !
  use constants
  !
  IMPLICIT NONE
  !
  INTEGER norb
  INTEGER nrpt
  !
  COMPLEX(DP), ALLOCATABLE :: ham(:,:,:)
  !  ham(norb, norb, nrpt)
  REAL(DP), ALLOCATABLE :: wgt(:)
  !  wgt(nrpt)
  REAL(DP), ALLOCATABLE :: rvec(:,:)
  !  rvec(1:3, nrpt)
CONTAINS

SUBROUTINE read_ham(fn)
!
  USE para
  USE constants
  !
  IMPLICIT NONE
  !
  CHARACTER(len=80) fn
  INTEGER irpt, ii, t1, t2, t3, t4, t5
  INTEGER, ALLOCATABLE :: wt(:)
  REAL(DP) a, b
  !
  if (inode.eq.0) then
    !
    write(stdout, *) " # Reading spin unpolarized Hamiltonian from "//trim(fn)
    !
    open(unit=fin, file=trim(fn))
    !
    read(fin, *)
    read(fin, *) t1
    read(fin, *) t2
    !
    norb = t1
    nrpt = t2
    write(stdout, *) " #  Dimensions:"
    write(stdout, *) "    # of orbitals:", norb
    write(stdout, *) "    # of real-space grid:", nrpt
    !
  endif  ! inode.eq.0
  !
  CALL para_sync(norb)
  CALL para_sync(nrpt)
  !
  allocate(ham(1:norb, 1:norb, 1:nrpt))
  allocate(wgt(1:nrpt))
  allocate(rvec(1:3, 1:nrpt))
  !
  if (inode.eq.0) then
    !
    allocate(wt(1:nrpt))
    read(fin, '(15I5)') (wt(irpt),irpt=1,nrpt)
    wgt(:)=wt(:)
    deallocate(wt)
    !
    do irpt=1, nrpt
      do ii=1, norb*norb
        read(fin, *) t1, t2, t3, t4, t5, a, b
        if (ii.eq.1) then
          rvec(1, irpt)=t1
          rvec(2, irpt)=t2
          rvec(3, irpt)=t3
        endif
      ham(t4, t5, irpt)=CMPLX(a,b)
      enddo
    enddo
    !
    close(unit=fin)
    write(stdout, *) " # Done."
    !
  endif ! inode.eq.0
  !
  CALL para_sync(ham, norb, norb, nrpt)
  CALL para_sync(wgt, nrpt)
  CALL para_sync(rvec, 3, nrpt)
  !
END SUBROUTINE

SUBROUTINE finalize_wann()
  !
  IMPLICIT NONE
  !
  if (allocated(ham)) deallocate(ham)
  if (allocated(wgt)) deallocate(wgt)
  if (allocated(rvec)) deallocate(rvec)
  !
END SUBROUTINE

SUBROUTINE calc_ek(eig, kvec)
  !
  use lapack95,  only : heev
  use constants, only : dp, twopi, cmplx_0, cmplx_i
  !
  implicit none
  !
  real(dp), dimension(1:norb) :: eig
  real(dp), dimension(1:3) :: kvec
  !
  real(dp)   rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:, :)
  !
  integer ir, info
  !
  allocate(work(1:norb, 1:norb))
  !
  work(:, :)=cmplx_0
  do ir=1, nrpt
    rdotk=SUM(kvec(:)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/wgt(ir)
    work(:, :)=work(:, :) + fact*ham(:, :, ir)
  enddo
  !
  call heev(work(:, :), eig, 'N', 'U', info)
  !
  deallocate(work)
  !
END SUBROUTINE

SUBROUTINE calc_dedk(v, kvec, idir)
  !
  use constants, only : dp, eps6
  !
  implicit none
  !
  real(dp), dimension(1:norb) :: v
  real(dp), dimension(1:3) :: kvec
  integer idir
  !
  real(dp), allocatable :: t(:, :)
  real(dp), dimension(1:3) :: kv1, kv2
  !
  allocate(t(1:norb, 1:2))
  !
  kv1=kvec
  kv2=kvec
  !
  kv1(idir)=kvec(idir)+eps6
  kv2(idir)=kvec(idir)-eps6
  !
  CALL calc_ek(t(:,1), kv1)
  CALL calc_ek(t(:,2), kv2)
  !
  v(:)=(t(:,1)-t(:,2))/(2*eps6)
  !
  deallocate(t)
  !
END SUBROUTINE

SUBROUTINE calc_d2edk2(mu, kvec, idir, jdir)
  !
  use constants, only : dp, eps6
  !
  implicit none
  !
  real(dp), dimension(1:norb) :: mu
  real(dp), dimension(1:3) :: kvec
  integer idir, jdir
  !
  real(dp), allocatable :: t(:, :)
  real(dp), dimension(1:3) :: kv
  !
  allocate(t(1:norb, 1:4))
  !
  if (idir.eq.jdir) then
    !
    kv=kvec
    CALL calc_ek(t(:,2), kv)
    !
    kv(idir)=kvec(idir)+eps6
    CALL calc_ek(t(:,3), kv)
    !
    kv(idir)=kvec(idir)-eps6
    CALL calc_ek(t(:,1), kv)
    !
    mu(:)=(t(:,1)+t(:,3)-2*t(:,2))*2.5E11
    !
  else
    !
    kv=kvec
    kv(idir)=kvec(idir)+eps6
    kv(jdir)=kvec(jdir)+eps6
    CALL calc_ek(t(:,1), kv)
    !
    kv(idir)=kvec(idir)-eps6
    CALL calc_ek(t(:,2), kv)
    !
    kv(jdir)=kvec(jdir)-eps6
    CALL calc_ek(t(:,4), kv)
    !
    kv(idir)=kvec(idir)+eps6
    CALL calc_ek(t(:,3), kv)
    !
    mu(:)=(t(:,1)+t(:,4)-t(:,2)-t(:,3))*2.5E11
    !
  endif
  !
  deallocate(t)
  !
END SUBROUTINE
  !
END MODULE

