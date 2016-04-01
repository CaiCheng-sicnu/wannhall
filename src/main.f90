PROGRAM wannhall
  !
  USE constants,  only : dp, stdin, stdout
  USE para,       only : init_para, first_k, last_k, finalize_para, distribute_k, para_merge, inode
  USE wanndata,   only : read_ham, calc_ek, calc_dedk, calc_d2edk2, finalize_wann
  USE bxsfdata,   only : read_bxsf, nkx, nky, nkz, kvec, ek, nbnd, finalize_bxsf
  USE input,      only : ef, sigma, idir, jdir, seed, read_input
  !
  IMPLICIT NONE
  !
  integer ik, nkpt, ii
  !
  real(dp) sigma_xx, sigma_xyz
  real(dp), allocatable :: vx(:), vy(:), mu_yy(:), mu_xy(:)
  real(dp) delta
  logical isfirst
  character(len=80) fn
  !
  CALL init_para
  !
  CALL read_input
  !
  fn=trim(seed)//"_hr.dat"
  CALL read_ham(fn)
  !
  fn=trim(seed)//".bxsf"
  CALL read_bxsf(fn)
  !
  nkpt=nkx*nky*nkz
  !
  CALL distribute_k(nkpt)
  !
  allocate(vx(1:nbnd), vy(1:nbnd), mu_yy(1:nbnd), mu_xy(1:nbnd))
  !
  sigma_xx=0.d0
  sigma_xyz=0.d0
  !
  if (inode.eq.0) write(stdout, *) " # Start actual calculation ..."
  !
  do ik=first_k, last_k
    !
    if (inode.eq.0) write(stdout, '(1A)', ADVANCE='no') "."
    isfirst=.true.
    do ii=1, nbnd
      if ( abs(ef-ek(ii, ik))<5.0*sigma ) then
        if (isfirst) then
          CALL calc_dedk(vx, kvec(:, ik), idir)
          CALL calc_dedk(vy, kvec(:, ik), jdir)
          CALL calc_d2edk2(mu_yy, kvec(:, ik), jdir, jdir)
          CALL calc_d2edk2(mu_xy, kvec(:, ik), idir, jdir)
          isfirst=.false.
        endif
        !
        sigma_xx=sigma_xx+vx(ii)*vx(ii)*delta(ek(ii, ik)-ef, sigma)
        sigma_xyz=sigma_xyz+vx(ii)*(vx(ii)*mu_yy(ii)-vy(ii)*mu_xy(ii))*delta(ek(ii, ik)-ef, sigma)
        !
      endif
    enddo
    !
  enddo
  !
  if (inode.eq.0) write(stdout, *) " Merging results..."
  !
  deallocate(vx, vy, mu_xy, mu_yy)
  !
  CALL para_merge(sigma_xx)
  CALL para_merge(sigma_xyz)
  !
  sigma_xx=sigma_xx/nkpt
  sigma_xyz=sigma_xyz/nkpt
  !
  if (inode.eq.0) then
    write(stdout, '(1A,1F22.16)') "Hall conductance: ", sigma_xx*sigma_xx/sigma_xyz
  endif
  !
  CALL finalize_wann
  CALL finalize_bxsf
  CALL finalize_para
  !
END PROGRAM
