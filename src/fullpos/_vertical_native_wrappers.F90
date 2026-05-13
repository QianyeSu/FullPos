subroutine fullpos_pressure_ppq_c(ncol, nlev, nout, values, ak, bk, ps, levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:)
  real(kind=jprb), allocatable :: field(:,:), out(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppq(kproma, kst, kprof, kflev, klevp, klolev, kppm, klevb, prpres, ldbelo, ldblow, &
                   prxp, prxpd, pqf, pqpp, ldbob, pbob0hack)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kst, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      real(kind=jprb), intent(in) :: prpres(kproma,klevp)
      logical, intent(in) :: ldbelo(kproma,klevp), ldblow(klevp)
      real(kind=jprb), intent(in) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: pqf(kproma,kflev)
      real(kind=jprb), intent(out) :: pqpp(kproma,klevp)
      logical, optional, intent(in) :: ldbob
      real(kind=jprb), optional, intent(in) :: pbob0hack
    end subroutine ppq
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), field(ncol,nlev), out(ncol,nout))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      field(j,k) = values(j,k)
    enddo
    do k = 1, nout
      prpres(j,k) = levels(k)
    enddo
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppq(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, prpres, ld_belo, ld_blow, prxp, prxpd, field, out)

  output = out
end subroutine fullpos_pressure_ppq_c


subroutine fullpos_pressure_ppuv_c(ncol, nlev, nout, u_values, v_values, ak, bk, ps, levels, u_output, v_output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: u_values(ncol, nlev), v_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), levels(nout)
  real(c_double), intent(out) :: u_output(ncol, nout), v_output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:), lnpres(:,:)
  real(kind=jprb), allocatable :: uf(:,:), vf(:,:), uout(:,:), vout(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppuv(kproma, kstart, kprof, kflev, klevp, klolev, kppm, klevb, ldbelo, ldblow, &
                    prpres, prxp, prxpd, puf, pvf, pupp, pvpp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      logical, intent(in) :: ldbelo(kproma,klevp), ldblow(klevp)
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: puf(kproma,kflev), pvf(kproma,kflev)
      real(kind=jprb), intent(out) :: pupp(kproma,klevp), pvpp(kproma,klevp)
    end subroutine ppuv
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), lnpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), uf(ncol,nlev), vf(ncol,nlev), uout(ncol,nout), vout(ncol,nout))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      uf(j,k) = u_values(j,k)
      vf(j,k) = v_values(j,k)
    enddo
    do k = 1, nout
      prpres(j,k) = levels(k)
      lnpres(j,k) = log(levels(k))
    enddo
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppuv(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, ld_belo, ld_blow, lnpres, prxp, prxpd, uf, vf, uout, vout)

  u_output = uout
  v_output = vout
end subroutine fullpos_pressure_ppuv_c


subroutine fullpos_pressure_ppt_c(ncol, nlev, nout, values, ak, bk, ps, levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:), lnpres(:,:)
  real(kind=jprb), allocatable :: field(:,:), out(:,:), porog(:), ptstar(:), pt0(:), r0(:,:)
  real(kind=jprb), allocatable :: sttf(:,:), stzf(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppt(kproma, kstart, kprof, kflev, klevp, klolev, kppm, klevb, pres, plnpres, &
                   ldbelo, ldbels, ldblow, ldbles, ldextr, porog, prxp, prxpd, ptstar, pt0, pr0, ptf, ptpp, &
                   psttf, pr2, prtpp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      real(kind=jprb), intent(in) :: pres(kproma,klevp), plnpres(kproma,klevp)
      logical, intent(in) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp), ldextr
      real(kind=jprb), intent(in) :: porog(kproma), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: ptstar(kproma), pt0(kproma), pr0(kproma,kflev), ptf(kproma,kflev)
      real(kind=jprb), intent(inout) :: ptpp(kproma,klevp)
      real(kind=jprb), optional, intent(in) :: psttf(kproma,0:kflev)
      real(kind=jprb), optional, intent(in) :: pr2(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: prtpp(kproma,klevp)
    end subroutine ppt
    subroutine ppsta(cdatm, kproma, kstart, kprof, klevp, klolev, prpres, plnpres, pstt, pstfi)
      use parkind1, only: jpim, jprb
      character(len=*), intent(in) :: cdatm
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, klevp, klolev
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), plnpres(kproma,klevp)
      real(kind=jprb), intent(out) :: pstt(kproma,klevp), pstfi(kproma,klevp)
    end subroutine ppsta
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), lnpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), field(ncol,nlev), out(ncol,nout), porog(ncol), ptstar(ncol), pt0(ncol), r0(ncol,nlev))
  allocate(sttf(ncol,0:nlev), stzf(ncol,0:nlev))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      field(j,k) = values(j,k)
      r0(j,k) = 287.0596736665907_jprb
    enddo
    do k = 1, nout
      prpres(j,k) = levels(k)
      lnpres(j,k) = log(levels(k))
      out(j,k) = 0.0_jprb
    enddo
    porog(j) = 0.0_jprb
    ptstar(j) = values(j,nlev)
    pt0(j) = ptstar(j)
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppsta('PPREF', ncol, 1, ncol, nlev + 1, 1, prxp(1,0,2), prxp(1,0,4), sttf, stzf)
  call ppt(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, prpres, lnpres, ld_belo, ld_bels, ld_blow, ld_bles, &
           .true., porog, prxp, prxpd, ptstar, pt0, r0, field, out, sttf)

  output = out
end subroutine fullpos_pressure_ppt_c


subroutine fullpos_column_pressure_ppq_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:)
  real(kind=jprb), allocatable :: field(:,:), out(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppq(kproma, kst, kprof, kflev, klevp, klolev, kppm, klevb, prpres, ldbelo, ldblow, &
                   prxp, prxpd, pqf, pqpp, ldbob, pbob0hack)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kst, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      real(kind=jprb), intent(in) :: prpres(kproma,klevp)
      logical, intent(in) :: ldbelo(kproma,klevp), ldblow(klevp)
      real(kind=jprb), intent(in) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: pqf(kproma,kflev)
      real(kind=jprb), intent(out) :: pqpp(kproma,klevp)
      logical, optional, intent(in) :: ldbob
      real(kind=jprb), optional, intent(in) :: pbob0hack
    end subroutine ppq
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), field(ncol,nlev), out(ncol,nout))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      field(j,k) = values(j,k)
    enddo
    do k = 1, nout
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
      prpres(j,k) = target_pressures(j,k)
    enddo
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppq(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, prpres, ld_belo, ld_blow, prxp, prxpd, field, out)

  output = out
end subroutine fullpos_column_pressure_ppq_c


subroutine fullpos_column_pressure_ppuv_c(ncol, nlev, nout, u_values, v_values, ak, bk, ps, target_pressures, u_output, v_output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: u_values(ncol, nlev), v_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
  real(c_double), intent(out) :: u_output(ncol, nout), v_output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:), lnpres(:,:)
  real(kind=jprb), allocatable :: uf(:,:), vf(:,:), uout(:,:), vout(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppuv(kproma, kstart, kprof, kflev, klevp, klolev, kppm, klevb, ldbelo, ldblow, &
                    prpres, prxp, prxpd, puf, pvf, pupp, pvpp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      logical, intent(in) :: ldbelo(kproma,klevp), ldblow(klevp)
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: puf(kproma,kflev), pvf(kproma,kflev)
      real(kind=jprb), intent(out) :: pupp(kproma,klevp), pvpp(kproma,klevp)
    end subroutine ppuv
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), lnpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), uf(ncol,nlev), vf(ncol,nlev), uout(ncol,nout), vout(ncol,nout))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      uf(j,k) = u_values(j,k)
      vf(j,k) = v_values(j,k)
    enddo
    do k = 1, nout
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
      prpres(j,k) = target_pressures(j,k)
      lnpres(j,k) = log(target_pressures(j,k))
    enddo
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppuv(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, ld_belo, ld_blow, lnpres, prxp, prxpd, uf, vf, uout, vout)

  u_output = uout
  v_output = vout
end subroutine fullpos_column_pressure_ppuv_c


subroutine fullpos_column_pressure_ppt_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:), lnpres(:,:)
  real(kind=jprb), allocatable :: field(:,:), out(:,:), porog(:), ptstar(:), pt0(:), r0(:,:)
  real(kind=jprb), allocatable :: sttf(:,:), stzf(:,:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                      klevb, ldbelo, ldbels, ldblow, ldbles)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
      logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
    end subroutine ppflev
    subroutine ppt(kproma, kstart, kprof, kflev, klevp, klolev, kppm, klevb, pres, plnpres, &
                   ldbelo, ldbels, ldblow, ldbles, ldextr, porog, prxp, prxpd, ptstar, pt0, pr0, ptf, ptpp, &
                   psttf, pr2, prtpp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, klolev, kppm
      integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
      real(kind=jprb), intent(in) :: pres(kproma,klevp), plnpres(kproma,klevp)
      logical, intent(in) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp), ldextr
      real(kind=jprb), intent(in) :: porog(kproma), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
      real(kind=jprb), intent(in) :: ptstar(kproma), pt0(kproma), pr0(kproma,kflev), ptf(kproma,kflev)
      real(kind=jprb), intent(inout) :: ptpp(kproma,klevp)
      real(kind=jprb), optional, intent(in) :: psttf(kproma,0:kflev)
      real(kind=jprb), optional, intent(in) :: pr2(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: prtpp(kproma,klevp)
    end subroutine ppt
    subroutine ppsta(cdatm, kproma, kstart, kprof, klevp, klolev, prpres, plnpres, pstt, pstfi)
      use parkind1, only: jpim, jprb
      character(len=*), intent(in) :: cdatm
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, klevp, klolev
      real(kind=jprb), intent(in) :: prpres(kproma,klevp), plnpres(kproma,klevp)
      real(kind=jprb), intent(out) :: pstt(kproma,klevp), pstfi(kproma,klevp)
    end subroutine ppsta
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), lnpres(ncol,nout), klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout))
  allocate(ld_blow(nout), ld_bles(nout), field(ncol,nlev), out(ncol,nout), porog(ncol), ptstar(ncol), pt0(ncol), r0(ncol,nlev))
  allocate(sttf(ncol,0:nlev), stzf(ncol,0:nlev))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      field(j,k) = values(j,k)
      r0(j,k) = 287.0596736665907_jprb
    enddo
    do k = 1, nout
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
      prpres(j,k) = target_pressures(j,k)
      lnpres(j,k) = log(target_pressures(j,k))
      out(j,k) = 0.0_jprb
    enddo
    porog(j) = 0.0_jprb
    ptstar(j) = values(j,nlev)
    pt0(j) = ptstar(j)
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)
  call ppflev(ncol, 1, ncol, nlev, nout, kppm, prpres, presh, presf, klevb, ld_belo, ld_bels, ld_blow, ld_bles)
  call ppsta('PPREF', ncol, 1, ncol, nlev + 1, 1, prxp(1,0,2), prxp(1,0,4), sttf, stzf)
  call ppt(ncol, 1, ncol, nlev, nout, 1, kppm, klevb, prpres, lnpres, ld_belo, ld_bels, ld_blow, ld_bles, &
           .true., porog, prxp, prxpd, ptstar, pt0, r0, field, out, sttf)

  output = out
end subroutine fullpos_column_pressure_ppt_c


subroutine fullpos_apache_column_pressure_tuvq_c(ncol, nlev, nout, temperature, u_values, v_values, q_values, &
                                                  ak, bk, ps, target_pressures, target_surface_pressure, enable_lescale, &
                                                  t_output, u_output, v_output, q_output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use intdyn_mod, only: yytxyb
  use yomcver, only: tcver
  use yomct0, only: lecmwf
  use yomppvi, only: lescale, lescale_t, lescale_q, lescale_u, lescale_pd, lescale_gfl, lrppuv_cstext, lrppuv_callitpq
  use yomphy, only: yrphy
  use yomvert, only: tvab, tveta, tvfe, tvertical_geom, vp00, toppres
  use type_gflflds, only: type_igflfld
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  integer(c_int), intent(in) :: enable_lescale
  real(c_double), intent(in) :: temperature(ncol, nlev), u_values(ncol, nlev), v_values(ncol, nlev), q_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout), target_surface_pressure(ncol)
  real(c_double), intent(out) :: t_output(ncol, nout), u_output(ncol, nout), v_output(ncol, nout), q_output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  integer(kind=jpim), parameter :: kppm = 4, kppgfl = 7
  type(type_igflfld) :: ydkc
  type(tvab) :: vab
  type(tveta) :: veta
  type(tvfe) :: vfe
  type(tcver) :: cver
  type(tvertical_geom) :: vertical_geom
  integer(kind=jpim), allocatable :: ilevb(:,:,:)
  logical, allocatable :: ldbelo(:,:), ldbels(:,:), ldblow(:), ldbles(:), ldgfl(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), xyb(:,:,:), prxp(:,:,:), prxpd(:,:,:), prpres(:,:), lnpres(:,:)
  real(kind=jprb), allocatable :: t1(:,:), u1(:,:), v1(:,:), q1(:,:), pspd1(:,:), pww1(:,:), pgfl1(:,:,:)
  real(kind=jprb), allocatable :: t2(:,:), u2(:,:), v2(:,:), pspd2(:,:), pwh2(:,:,:), pwf2(:,:), pgfl2(:,:,:), pfi2(:,:), prh2(:,:), peorg(:,:)
  real(kind=jprb), allocatable :: pts1(:), pts2(:), peps2(:), poro1(:), poro2(:), prhbnds(:), pr1(:,:), pkapa1(:,:), prhf1(:,:)
  real(kind=jprb), allocatable :: pepr2f(:,:), pepr2h(:,:,:), plnpr1(:,:), palph1(:,:)

  interface
    subroutine gprcp(kproma, kstart, kprof, kflev, pq, pqi, pql, pqr, pqs, pqg, pcp, pr, pkap, pgfl, kgfltyp, ldthermact)
      use parkind1, only: jpim, jprb
      use yom_ygfl, only: ygfl
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), optional, intent(in) :: pq(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqi(kproma,kflev), pql(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqr(kproma,kflev), pqs(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqg(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pgfl(kproma,kflev,ygfl%ndim)
      real(kind=jprb), optional, intent(out) :: pcp(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: pkap(kproma,kflev)
      integer(kind=jpim), optional, intent(in) :: kgfltyp
      logical, optional, intent(in) :: ldthermact
    end subroutine gprcp
    subroutine gphpre(kproma, kflev, kstart, kprof, ydvab, ydcver, presh, pxyb, presf)
      use parkind1, only: jpim, jprb
      use yomcver, only: tcver
      use yomvert, only: tvab
      use intdyn_mod, only: yytxyb
      integer(kind=jpim), intent(in) :: kproma, kflev, kstart, kprof
      type(tvab), intent(in) :: ydvab
      type(tcver), intent(in) :: ydcver
      real(kind=jprb), intent(inout) :: presh(kproma,0:kflev)
      real(kind=jprb), optional, intent(out) :: pxyb(kproma,kflev,yytxyb%ndim)
      real(kind=jprb), optional, intent(out) :: presf(kproma,kflev)
    end subroutine gphpre
    subroutine apache(ldneige,kproma,kstart,kprof,klev,koplev,kppgfl,kppm,ydkc,ldpdep,ldwhalf,ldwfull,ldgfl,ldobs,ldspch,prhbnds,&
                      poro1,poro2,pdpblc,pts1,pts2,peps2,pepr2f,pepr2h,ydvert_geom,pt1,pu1,pv1,pspd1,pww1,pgfl1,&
                      pres1h,pres1f,plnpr1,palph1,pr1,pkapa1,prhf1,pt2,pu2,pv2,pspd2,pwh2,pwf2,pgfl2,pfi2,prh2,peorg)
      use parkind1, only: jpim, jprb
      use type_gflflds, only: type_igflfld
      use yomvert, only: tvertical_geom
      logical, intent(in) :: ldneige, ldpdep, ldwhalf, ldwfull, ldobs, ldspch
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, klev, koplev, kppgfl, kppm
      type(type_igflfld), intent(in) :: ydkc
      logical, intent(in) :: ldgfl(kppgfl)
      real(kind=jprb), intent(in) :: prhbnds(2), poro1(kproma), poro2(kproma), pdpblc
      real(kind=jprb), intent(in) :: pts1(kproma), pts2(kproma), peps2(kproma)
      real(kind=jprb), intent(in) :: pepr2f(kproma,koplev), pepr2h(kproma,koplev,0:1)
      type(tvertical_geom), intent(in) :: ydvert_geom
      real(kind=jprb), intent(in) :: pt1(kproma,klev), pu1(kproma,klev), pv1(kproma,klev)
      real(kind=jprb), intent(inout) :: pspd1(kproma,klev), pww1(kproma,0:klev), pgfl1(kproma,klev,kppgfl)
      real(kind=jprb), intent(in) :: pres1h(kproma,0:klev), pres1f(kproma,klev), plnpr1(kproma,klev), palph1(kproma,klev)
      real(kind=jprb), intent(in) :: pr1(kproma,klev), pkapa1(kproma,klev), prhf1(kproma,klev)
      real(kind=jprb), intent(inout) :: pt2(kproma,koplev), pu2(kproma,koplev), pv2(kproma,koplev), pspd2(kproma,koplev)
      real(kind=jprb), intent(inout) :: pwh2(kproma,koplev,0:1), pwf2(kproma,koplev), pgfl2(kproma,koplev,kppgfl)
      real(kind=jprb), intent(inout) :: pfi2(kproma,koplev), prh2(kproma,koplev)
      real(kind=jprb), intent(out) :: peorg(kproma,koplev)
    end subroutine apache
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  lescale = enable_lescale /= 0_c_int
  lescale_t = lescale
  lescale_q = lescale
  lescale_u = lescale
  lescale_pd = .false.
  lescale_gfl = .false.
  lrppuv_cstext = .false.
  lrppuv_callitpq = .true.
  yrphy%lneige = .false.
  lecmwf = .true.

  allocate(vab%vah(0:nlev), vab%vbh(0:nlev), vab%valh(0:nlev), vab%vc(nlev), vab%vaf(0:nlev))
  allocate(vab%vbf(0:nlev), vab%vdela(nlev), vab%vdelb(nlev))
  allocate(veta%vfe_rdetah(0:nlev))
  allocate(presh(ncol,0:nlev), presf(ncol,nlev), xyb(ncol,nlev,yytxyb%ndim), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(prpres(ncol,nout), lnpres(ncol,nout), ilevb(ncol,nout,kppm), ldbelo(ncol,nout), ldbels(ncol,nout))
  allocate(ldblow(nout), ldbles(nout), ldgfl(kppgfl))
  allocate(t1(ncol,nlev), u1(ncol,nlev), v1(ncol,nlev), q1(ncol,nlev), pspd1(ncol,nlev), pww1(ncol,0:nlev), pgfl1(ncol,nlev,kppgfl))
  allocate(t2(ncol,nout), u2(ncol,nout), v2(ncol,nout), pspd2(ncol,nout), pwh2(ncol,nout,0:1), pwf2(ncol,nout), pgfl2(ncol,nout,kppgfl), pfi2(ncol,nout), prh2(ncol,nout), peorg(ncol,nout))
  allocate(pts1(ncol), pts2(ncol), peps2(ncol), poro1(ncol), poro2(ncol), prhbnds(2), pr1(ncol,nlev), pkapa1(ncol,nlev), prhf1(ncol,nlev))
  allocate(pepr2f(ncol,nout), pepr2h(ncol,nout,0:1), plnpr1(ncol,nlev), palph1(ncol,nlev))

  vp00 = 100000.0_jprb
  toppres = 0.1_jprb
  cver%laprxpk = .true.
  cver%ndlnpr = 0
  cver%rhydr0 = log(2.0_jprb)
  cver%lvertfe = .false.
  vertical_geom%lnonhyd_geom = .false.
  vertical_geom%yrvab = vab
  vertical_geom%yrveta = veta
  vertical_geom%yrvfe = vfe
  vertical_geom%yrcver = cver

  ydkc%iq = 1
  ydkc%il = 2
  ydkc%ii = 3
  ydkc%irr = 4
  ydkc%is = 5
  ydkc%ig = 6
  ydkc%ih = 7
  ldgfl = .false.
  ldgfl(1) = .true.
  prhbnds(1) = 0.0_jprb
  prhbnds(2) = 1.0_jprb
  pspd1 = 0.0_jprb
  pww1 = 0.0_jprb
  pgfl1 = 0.0_jprb
  t2 = 0.0_jprb
  u2 = 0.0_jprb
  v2 = 0.0_jprb
  pspd2 = 0.0_jprb
  pwh2 = 0.0_jprb
  pwf2 = 0.0_jprb
  pgfl2 = 0.0_jprb
  pfi2 = 0.0_jprb
  prh2 = 0.0_jprb
  peorg = 0.0_jprb

  do k = 0, nlev
    vab%vah(k) = ak(k + 1)
    vab%vbh(k) = bk(k + 1)
    vab%valh(k) = vab%vah(k) / vp00
    if (k > 0) then
      vab%vc(k) = vab%vah(k) * vab%vbh(k - 1) - vab%vah(k - 1) * vab%vbh(k)
      vab%vaf(k) = 0.5_jprb * (vab%vah(k) + vab%vah(k - 1))
      vab%vbf(k) = 0.5_jprb * (vab%vbh(k) + vab%vbh(k - 1))
      vab%vdela(k) = vab%vah(k) - vab%vah(k - 1)
      vab%vdelb(k) = vab%vbh(k) - vab%vbh(k - 1)
    endif
  enddo
  vab%vaf(0) = vab%vah(0)
  vab%vbf(0) = vab%vbh(0)
  vertical_geom%yrvab = vab

  do j = 1, ncol
    if (ps(j) <= 0.0_c_double) then
      ierr = 2_c_int
      return
    endif
    if (target_surface_pressure(j) <= 0.0_c_double) then
      ierr = 4_c_int
      return
    endif
    poro1(j) = 0.0_jprb
    poro2(j) = 0.0_jprb
    pts1(j) = temperature(j,nlev)
    pts2(j) = temperature(j,nlev)
    presh(j,nlev) = ps(j)
    peps2(j) = target_surface_pressure(j)
    do k = 1, nlev
      t1(j,k) = temperature(j,k)
      u1(j,k) = u_values(j,k)
      v1(j,k) = v_values(j,k)
      q1(j,k) = q_values(j,k)
      pspd1(j,k) = 0.0_jprb
      prhf1(j,k) = 0.0_jprb
      pgfl1(j,k,1) = q_values(j,k)
    enddo
  enddo

  call gphpre(ncol, nlev, 1, ncol, vab, cver, presh, pxyb=xyb, presf=presf)
  plnpr1 = xyb(:,:,yytxyb%m_lnpr)
  palph1 = xyb(:,:,yytxyb%m_alph)

  call gprcp(ncol, 1, ncol, nlev, pq=q1, pr=pr1, pkap=pkapa1, ldthermact=.true.)

  do j = 1, ncol
    do k = 1, nout
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 3_c_int
        return
      endif
      pepr2f(j,k) = target_pressures(j,k)
      pepr2h(j,k,0) = target_pressures(j,k)
      pepr2h(j,k,1) = target_pressures(j,k)
      prpres(j,k) = target_pressures(j,k)
      lnpres(j,k) = log(target_pressures(j,k))
    enddo
  enddo

  call apache(.false., ncol, 1, ncol, nlev, nout, kppgfl, kppm, ydkc, &
              .false., .false., .false., ldgfl, .false., .true., prhbnds, &
              poro1, poro2, 0.0_jprb, pts1, pts2, peps2, pepr2f, pepr2h, &
              vertical_geom, &
              t1, u1, v1, pspd1, pww1, pgfl1, &
              presh, presf, plnpr1, palph1, pr1, pkapa1, prhf1, &
              t2, u2, v2, pspd2, pwh2, pwf2, pgfl2, pfi2, prh2, peorg)

  t_output = t2
  u_output = u2
  v_output = v2
  q_output = pgfl2(:,:,1)
end subroutine fullpos_apache_column_pressure_tuvq_c


subroutine fullpos_hybrid_pressure_ppq_c(ncol, nlev, nout, values, ak, bk, ps, target_ak, target_bk, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol)
  real(c_double), intent(in) :: target_ak(nout + 1), target_bk(nout + 1)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: j, k
  real(c_double), allocatable :: target_pressures(:,:)

  interface
    subroutine fullpos_column_pressure_ppq_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr) bind(C)
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in) :: ncol, nlev, nout
      real(c_double), intent(in) :: values(ncol, nlev)
      real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
      real(c_double), intent(out) :: output(ncol, nout)
      integer(c_int), intent(out) :: ierr
    end subroutine fullpos_column_pressure_ppq_c
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(target_pressures(ncol,nout))
  do j = 1, ncol
    do k = 1, nout
      target_pressures(j,k) = 0.5_c_double * &
        (target_ak(k) + target_bk(k) * ps(j) + target_ak(k + 1) + target_bk(k + 1) * ps(j))
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
    enddo
  enddo

  call fullpos_column_pressure_ppq_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr)
end subroutine fullpos_hybrid_pressure_ppq_c


subroutine fullpos_hybrid_pressure_ppuv_c(ncol, nlev, nout, u_values, v_values, ak, bk, ps, target_ak, target_bk, &
                                          u_output, v_output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: u_values(ncol, nlev), v_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol)
  real(c_double), intent(in) :: target_ak(nout + 1), target_bk(nout + 1)
  real(c_double), intent(out) :: u_output(ncol, nout), v_output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: j, k
  real(c_double), allocatable :: target_pressures(:,:)

  interface
    subroutine fullpos_column_pressure_ppuv_c(ncol, nlev, nout, u_values, v_values, ak, bk, ps, target_pressures, &
                                              u_output, v_output, ierr) bind(C)
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in) :: ncol, nlev, nout
      real(c_double), intent(in) :: u_values(ncol, nlev), v_values(ncol, nlev)
      real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
      real(c_double), intent(out) :: u_output(ncol, nout), v_output(ncol, nout)
      integer(c_int), intent(out) :: ierr
    end subroutine fullpos_column_pressure_ppuv_c
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(target_pressures(ncol,nout))
  do j = 1, ncol
    do k = 1, nout
      target_pressures(j,k) = 0.5_c_double * &
        (target_ak(k) + target_bk(k) * ps(j) + target_ak(k + 1) + target_bk(k + 1) * ps(j))
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
    enddo
  enddo

  call fullpos_column_pressure_ppuv_c(ncol, nlev, nout, u_values, v_values, ak, bk, ps, target_pressures, &
                                      u_output, v_output, ierr)
end subroutine fullpos_hybrid_pressure_ppuv_c


subroutine fullpos_hybrid_pressure_ppt_c(ncol, nlev, nout, values, ak, bk, ps, target_ak, target_bk, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol)
  real(c_double), intent(in) :: target_ak(nout + 1), target_bk(nout + 1)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: j, k
  real(c_double), allocatable :: target_pressures(:,:)

  interface
    subroutine fullpos_column_pressure_ppt_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr) bind(C)
      use iso_c_binding, only: c_int, c_double
      integer(c_int), intent(in) :: ncol, nlev, nout
      real(c_double), intent(in) :: values(ncol, nlev)
      real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), target_pressures(ncol, nout)
      real(c_double), intent(out) :: output(ncol, nout)
      integer(c_int), intent(out) :: ierr
    end subroutine fullpos_column_pressure_ppt_c
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(target_pressures(ncol,nout))
  do j = 1, ncol
    do k = 1, nout
      target_pressures(j,k) = 0.5_c_double * &
        (target_ak(k) + target_bk(k) * ps(j) + target_ak(k + 1) + target_bk(k + 1) * ps(j))
      if (target_pressures(j,k) <= 0.0_c_double) then
        ierr = 2_c_int
        return
      endif
    enddo
  enddo

  call fullpos_column_pressure_ppt_c(ncol, nlev, nout, values, ak, bk, ps, target_pressures, output, ierr)
end subroutine fullpos_hybrid_pressure_ppt_c


subroutine fullpos_eta_pressures_c(ncol, nlev, nout, ak, bk, ps, eta_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use yomcver, only: tcver
  use yomvert, only: tvab, vp00, toppres
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), eta_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, ilev
  type(tvab) :: vab
  type(tcver) :: cver
  real(kind=jprb), allocatable :: levels(:), surface(:)
  real(kind=jprb), allocatable :: presf(:,:), preshb(:,:), presht(:,:)

  interface
    subroutine ppleta(ydvab, ydcver, kproma, kst, knd, koplev, pxlev, ps, presf, preshb, presht)
      use parkind1, only: jpim, jprb
      use yomcver, only: tcver
      use yomvert, only: tvab
      type(tvab), intent(in) :: ydvab
      type(tcver), intent(in) :: ydcver
      integer(kind=jpim), intent(in) :: kproma, kst, knd, koplev
      real(kind=jprb), intent(in) :: pxlev(koplev), ps(kproma)
      real(kind=jprb), intent(out) :: presf(kproma,koplev), preshb(kproma,koplev), presht(kproma,koplev)
    end subroutine ppleta
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(vab%vah(0:nlev), vab%vbh(0:nlev), vab%valh(0:nlev), vab%vc(nlev), vab%vaf(0:nlev))
  allocate(vab%vbf(0:nlev), vab%vdela(nlev), vab%vdelb(nlev))
  allocate(levels(nout), surface(ncol), presf(ncol,nout), preshb(ncol,nout), presht(ncol,nout))

  vp00 = 100000.0_jprb
  toppres = 0.1_jprb
  cver%laprxpk = .true.
  cver%ndlnpr = 0
  cver%rhydr0 = log(2.0_jprb)
  cver%lvertfe = .false.

  do k = 0, nlev
    vab%vah(k) = ak(k + 1)
    vab%vbh(k) = bk(k + 1)
    vab%valh(k) = vab%vah(k) / vp00
  enddo
  vab%vaf(0) = vab%vah(0)
  vab%vbf(0) = vab%vbh(0)
  do k = 1, nlev
    vab%vc(k) = vab%vah(k) * vab%vbh(k - 1) - vab%vah(k - 1) * vab%vbh(k)
    vab%vaf(k) = 0.5_jprb * (vab%vah(k) + vab%vah(k - 1))
    vab%vbf(k) = 0.5_jprb * (vab%vbh(k) + vab%vbh(k - 1))
    vab%vdela(k) = vab%vah(k) - vab%vah(k - 1)
    vab%vdelb(k) = vab%vbh(k) - vab%vbh(k - 1)
  enddo

  do j = 1, ncol
    if (ps(j) <= 0.0_c_double) then
      ierr = 2_c_int
      return
    endif
    surface(j) = ps(j)
  enddo
  do k = 1, nout
    ilev = int(eta_levels(k), kind=jpim)
    if (ilev < 1 .or. ilev > nlev) then
      ierr = 3_c_int
      return
    endif
    if (abs(eta_levels(k) - real(ilev, kind=c_double)) > 1.0e-9_c_double) then
      ierr = 4_c_int
      return
    endif
    levels(k) = real(ilev, kind=jprb)
  enddo

  call ppleta(vab, cver, ncol, 1, ncol, nout, levels, surface, presf, preshb, presht)

  do k = 1, nout
    do j = 1, ncol
      if (presf(j,k) <= 0.0_jprb) then
        ierr = 5_c_int
        return
      endif
      output(j,k) = presf(j,k)
    enddo
  enddo
end subroutine fullpos_eta_pressures_c


subroutine fullpos_theta_pressures_c(ncol, nlev, nout, temperature, ak, bk, ps, theta_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: temperature(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), theta_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, lev
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), temp(:,:), kappa(:,:), theta(:,:)
  real(kind=jprb), allocatable :: target_pressure(:), target_lnpressure(:)

  interface
    subroutine gptet(kproma, kstart, kprof, kflev, presf, pt, pkap, pteta)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), intent(in) :: presf(kproma,kflev), pt(kproma,kflev), pkap(kproma,kflev)
      real(kind=jprb), intent(out) :: pteta(kproma,kflev)
    end subroutine gptet
    subroutine pplteta(kproma, kstart, kprof, kflev, prpresf, pteta, pxteta, prpres, prlnpres)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), intent(in) :: prpresf(kproma,kflev), pteta(kproma,kflev), pxteta
      real(kind=jprb), intent(inout) :: prpres(kproma)
      real(kind=jprb), intent(out) :: prlnpres(kproma)
    end subroutine pplteta
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), temp(ncol,nlev), kappa(ncol,nlev), theta(ncol,nlev))
  allocate(target_pressure(ncol), target_lnpressure(ncol))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) then
        ierr = 2_c_int
        return
      endif
      temp(j,k) = temperature(j,k)
      kappa(j,k) = 287.0596736665907_jprb / 1004.7095955714683_jprb
    enddo
  enddo

  call gptet(ncol, 1, ncol, nlev, presf, temp, kappa, theta)

  do lev = 1, nout
    if (theta_levels(lev) <= 0.0_c_double) then
      ierr = 3_c_int
      return
    endif
    target_pressure(:) = 0.0_jprb
    call pplteta(ncol, 1, ncol, nlev, presf, theta, theta_levels(lev), target_pressure, target_lnpressure)
    do j = 1, ncol
      if (target_pressure(j) <= 0.0_jprb) then
        ierr = 4_c_int
        return
      endif
      output(j,lev) = target_pressure(j)
    enddo
  enddo
end subroutine fullpos_theta_pressures_c


subroutine fullpos_temperature_pressures_c(ncol, nlev, nout, temperature, ak, bk, ps, temperature_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use yomsta, only: tsta
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: temperature(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), temperature_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, lev
  type(tsta) :: sta
  real(kind=jprb), allocatable :: presh(:,:), temp(:,:), rgas(:,:), geoph(:,:), geopf(:,:), out_geop(:,:), pst(:,:)

  interface
    subroutine ppltw(ydsta, kproma, kst, knd, kflev, pgeo, pt, pxtemp, pxgeo)
      use parkind1, only: jpim, jprb
      use yomsta, only: tsta
      type(tsta), intent(in) :: ydsta
      integer(kind=jpim), intent(in) :: kproma, kst, knd, kflev
      real(kind=jprb), intent(in) :: pgeo(kproma,kflev), pt(kproma,kflev), pxtemp
      real(kind=jprb), intent(out) :: pxgeo(kproma)
    end subroutine ppltw
    subroutine fpps(kproma, kst, knd, kflev, koplev, pin_geoph, pout_geop, pt, pr, presh, pst, psppp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kst, knd, kflev, koplev
      real(kind=jprb), intent(in) :: pin_geoph(kproma,0:kflev), pout_geop(kproma,koplev)
      real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev), presh(kproma,0:kflev), pst(kproma)
      real(kind=jprb), intent(out) :: psppp(kproma,koplev)
    end subroutine fpps
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), temp(ncol,nlev), rgas(ncol,nlev), geoph(ncol,0:nlev), geopf(ncol,nlev))
  allocate(out_geop(ncol,nout), pst(ncol,nout), sta%stz(nlev))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    geoph(j,nlev) = 0.0_jprb
    do k = 1, nlev
      temp(j,k) = temperature(j,k)
      rgas(j,k) = 287.0596736665907_jprb
    enddo
    do k = nlev, 1, -1
      geoph(j,k - 1) = geoph(j,k) + rgas(j,k) * temp(j,k) * log(presh(j,k) / presh(j,k - 1))
      geopf(j,k) = 0.5_jprb * (geoph(j,k - 1) + geoph(j,k))
    enddo
  enddo
  do k = 1, nlev
    sta%stz(k) = 0.5_jprb * (geopf(1,k) / 9.80665_jprb)
  enddo

  do lev = 1, nout
    call ppltw(sta, ncol, 1, ncol, nlev, geopf, temp, temperature_levels(lev), out_geop(1,lev))
  enddo

  call fpps(ncol, 1, ncol, nlev, nout, geoph, out_geop, temp, rgas, presh, temp(1,nlev), pst)

  do lev = 1, nout
    do j = 1, ncol
      if (pst(j,lev) <= 0.0_jprb) then
        ierr = 2_c_int
        return
      endif
      output(j,lev) = pst(j,lev)
    enddo
  enddo
end subroutine fullpos_temperature_pressures_c


subroutine fullpos_height_above_orography_pressures_c(ncol, nlev, nout, temperature, specific_humidity, &
                                                       has_specific_humidity, ak, bk, ps, orography_geopotential, &
                                                       height_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use intdyn_mod, only: yytxyb
  use yomcst, only: rd, rg
  use yomcver, only: tcver
  use yomvert, only: tvab, tvertical_geom, vp00, toppres
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout, has_specific_humidity
  real(c_double), intent(in) :: temperature(ncol, nlev), specific_humidity(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), orography_geopotential(ncol)
  real(c_double), intent(in) :: height_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, lev
  type(tvab) :: vab
  type(tcver) :: cver
  type(tvertical_geom) :: vertical_geom
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), xyb(:,:,:)
  real(kind=jprb), allocatable :: temp(:,:), q(:,:), rgas(:,:), kappa(:,:)
  real(kind=jprb), allocatable :: geoph(:,:), geopf(:,:), target_geop(:,:), surface_temp(:), target_pressure(:,:)

  interface
    subroutine gphpre(kproma, kflev, kstart, kprof, ydvab, ydcver, presh, pxyb, presf)
      use parkind1, only: jpim, jprb
      use yomcver, only: tcver
      use yomvert, only: tvab
      use intdyn_mod, only: yytxyb
      integer(kind=jpim), intent(in) :: kproma, kflev, kstart, kprof
      type(tvab), intent(in) :: ydvab
      type(tcver), intent(in) :: ydcver
      real(kind=jprb), intent(inout) :: presh(kproma,0:kflev)
      real(kind=jprb), optional, intent(out) :: pxyb(kproma,kflev,yytxyb%ndim)
      real(kind=jprb), optional, intent(out) :: presf(kproma,kflev)
    end subroutine gphpre
    subroutine gprcp(kproma, kstart, kprof, kflev, pq, pqi, pql, pqr, pqs, pqg, &
                     pcp, pr, pkap, pgfl, kgfltyp, ldthermact)
      use parkind1, only: jpim, jprb
      use yom_ygfl, only: ygfl
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), optional, intent(in) :: pq(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqi(kproma,kflev), pql(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqr(kproma,kflev), pqs(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqg(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pgfl(kproma,kflev,ygfl%ndim)
      real(kind=jprb), optional, intent(out) :: pcp(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: pkap(kproma,kflev)
      integer(kind=jpim), optional, intent(in) :: kgfltyp
      logical, optional, intent(in) :: ldthermact
    end subroutine gprcp
    subroutine gpgeo(kproma, kstart, kprof, kflev, phi, phif, pt, pr, plnpr, palph, pvgeom)
      use parkind1, only: jpim, jprb
      use yomvert, only: tvertical_geom
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), intent(inout) :: phi(kproma,0:kflev)
      real(kind=jprb), intent(out) :: phif(kproma,kflev)
      real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), intent(in) :: plnpr(kproma,kflev), palph(kproma,kflev)
      type(tvertical_geom), intent(in) :: pvgeom
    end subroutine gpgeo
    subroutine fpps(kproma, kst, knd, kflev, koplev, pin_geoph, pout_geop, pt, pr, presh, pst, psppp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kst, knd, kflev, koplev
      real(kind=jprb), intent(in) :: pin_geoph(kproma,0:kflev), pout_geop(kproma,koplev)
      real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev), presh(kproma,0:kflev), pst(kproma)
      real(kind=jprb), intent(out) :: psppp(kproma,koplev)
    end subroutine fpps
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(vab%vah(0:nlev), vab%vbh(0:nlev), vab%valh(0:nlev), vab%vc(nlev), vab%vaf(0:nlev))
  allocate(vab%vbf(0:nlev), vab%vdela(nlev), vab%vdelb(nlev))
  allocate(presh(ncol,0:nlev), presf(ncol,nlev), xyb(ncol,nlev,yytxyb%ndim))
  allocate(temp(ncol,nlev), q(ncol,nlev), rgas(ncol,nlev), kappa(ncol,nlev))
  allocate(geoph(ncol,0:nlev), geopf(ncol,nlev), target_geop(ncol,nout), surface_temp(ncol), target_pressure(ncol,nout))

  vp00 = 100000.0_jprb
  toppres = 0.1_jprb
  cver%laprxpk = .true.
  cver%ndlnpr = 0
  cver%rhydr0 = log(2.0_jprb)
  cver%lvertfe = .false.
  vertical_geom%lnonhyd_geom = .false.
  vertical_geom%yrcver = cver

  do k = 0, nlev
    vab%vah(k) = ak(k + 1)
    vab%vbh(k) = bk(k + 1)
    vab%valh(k) = vab%vah(k) / vp00
  enddo
  vab%vaf(0) = vab%vah(0)
  vab%vbf(0) = vab%vbh(0)
  do k = 1, nlev
    vab%vc(k) = vab%vah(k) * vab%vbh(k - 1) - vab%vah(k - 1) * vab%vbh(k)
    vab%vaf(k) = 0.5_jprb * (vab%vah(k) + vab%vah(k - 1))
    vab%vbf(k) = 0.5_jprb * (vab%vbh(k) + vab%vbh(k - 1))
    vab%vdela(k) = vab%vah(k) - vab%vah(k - 1)
    vab%vdelb(k) = vab%vbh(k) - vab%vbh(k - 1)
  enddo
  vertical_geom%yrvab = vab

  do j = 1, ncol
    if (ps(j) <= 0.0_c_double) then
      ierr = 2_c_int
      return
    endif
    presh(j,nlev) = ps(j)
    geoph(j,nlev) = orography_geopotential(j)
    surface_temp(j) = temperature(j,nlev)
    do k = 1, nlev
      temp(j,k) = temperature(j,k)
      q(j,k) = specific_humidity(j,k)
    enddo
  enddo

  do lev = 1, nout
    if (height_levels(lev) < 0.0_c_double) then
      ierr = 3_c_int
      return
    endif
    do j = 1, ncol
      target_geop(j,lev) = orography_geopotential(j) + rg * height_levels(lev)
    enddo
  enddo

  call gphpre(ncol, nlev, 1, ncol, vab, cver, presh, pxyb=xyb, presf=presf)

  if (has_specific_humidity /= 0_c_int) then
    call gprcp(ncol, 1, ncol, nlev, pq=q, pr=rgas, pkap=kappa, ldthermact=.true.)
  else
    rgas(:,:) = rd
  endif

  call gpgeo(ncol, 1, ncol, nlev, geoph, geopf, temp, rgas, &
             xyb(:,:,yytxyb%m_lnpr), xyb(:,:,yytxyb%m_alph), vertical_geom)
  call fpps(ncol, 1, ncol, nlev, nout, geoph, target_geop, temp, rgas, presh, surface_temp, target_pressure)

  do lev = 1, nout
    do j = 1, ncol
      if (target_pressure(j,lev) <= 0.0_jprb) then
        ierr = 4_c_int
        return
      endif
      output(j,lev) = target_pressure(j,lev)
    enddo
  enddo
end subroutine fullpos_height_above_orography_pressures_c


subroutine fullpos_height_above_sea_pressures_c(ncol, nlev, nout, temperature, specific_humidity, &
                                                has_specific_humidity, ak, bk, ps, orography_geopotential, &
                                                height_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use intdyn_mod, only: yytxyb
  use yomcst, only: rd, rg
  use yomcver, only: tcver
  use yomvert, only: tvab, tvertical_geom, vp00, toppres
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout, has_specific_humidity
  real(c_double), intent(in) :: temperature(ncol, nlev), specific_humidity(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), orography_geopotential(ncol)
  real(c_double), intent(in) :: height_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, lev
  type(tvab) :: vab
  type(tcver) :: cver
  type(tvertical_geom) :: vertical_geom
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), xyb(:,:,:)
  real(kind=jprb), allocatable :: temp(:,:), q(:,:), rgas(:,:), kappa(:,:)
  real(kind=jprb), allocatable :: geoph(:,:), geopf(:,:), target_geop(:,:), surface_temp(:), target_pressure(:,:)

  interface
    subroutine gphpre(kproma, kflev, kstart, kprof, ydvab, ydcver, presh, pxyb, presf)
      use parkind1, only: jpim, jprb
      use yomcver, only: tcver
      use yomvert, only: tvab
      use intdyn_mod, only: yytxyb
      integer(kind=jpim), intent(in) :: kproma, kflev, kstart, kprof
      type(tvab), intent(in) :: ydvab
      type(tcver), intent(in) :: ydcver
      real(kind=jprb), intent(inout) :: presh(kproma,0:kflev)
      real(kind=jprb), optional, intent(out) :: pxyb(kproma,kflev,yytxyb%ndim)
      real(kind=jprb), optional, intent(out) :: presf(kproma,kflev)
    end subroutine gphpre
    subroutine gprcp(kproma, kstart, kprof, kflev, pq, pqi, pql, pqr, pqs, pqg, &
                     pcp, pr, pkap, pgfl, kgfltyp, ldthermact)
      use parkind1, only: jpim, jprb
      use yom_ygfl, only: ygfl
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), optional, intent(in) :: pq(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqi(kproma,kflev), pql(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqr(kproma,kflev), pqs(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqg(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pgfl(kproma,kflev,ygfl%ndim)
      real(kind=jprb), optional, intent(out) :: pcp(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: pkap(kproma,kflev)
      integer(kind=jpim), optional, intent(in) :: kgfltyp
      logical, optional, intent(in) :: ldthermact
    end subroutine gprcp
    subroutine gpgeo(kproma, kstart, kprof, kflev, phi, phif, pt, pr, plnpr, palph, pvgeom)
      use parkind1, only: jpim, jprb
      use yomvert, only: tvertical_geom
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), intent(inout) :: phi(kproma,0:kflev)
      real(kind=jprb), intent(out) :: phif(kproma,kflev)
      real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), intent(in) :: plnpr(kproma,kflev), palph(kproma,kflev)
      type(tvertical_geom), intent(in) :: pvgeom
    end subroutine gpgeo
    subroutine fpps(kproma, kst, knd, kflev, koplev, pin_geoph, pout_geop, pt, pr, presh, pst, psppp)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kst, knd, kflev, koplev
      real(kind=jprb), intent(in) :: pin_geoph(kproma,0:kflev), pout_geop(kproma,koplev)
      real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev), presh(kproma,0:kflev), pst(kproma)
      real(kind=jprb), intent(out) :: psppp(kproma,koplev)
    end subroutine fpps
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(vab%vah(0:nlev), vab%vbh(0:nlev), vab%valh(0:nlev), vab%vc(nlev), vab%vaf(0:nlev))
  allocate(vab%vbf(0:nlev), vab%vdela(nlev), vab%vdelb(nlev))
  allocate(presh(ncol,0:nlev), presf(ncol,nlev), xyb(ncol,nlev,yytxyb%ndim))
  allocate(temp(ncol,nlev), q(ncol,nlev), rgas(ncol,nlev), kappa(ncol,nlev))
  allocate(geoph(ncol,0:nlev), geopf(ncol,nlev), target_geop(ncol,nout), surface_temp(ncol), target_pressure(ncol,nout))

  vp00 = 100000.0_jprb
  toppres = 0.1_jprb
  cver%laprxpk = .true.
  cver%ndlnpr = 0
  cver%rhydr0 = log(2.0_jprb)
  cver%lvertfe = .false.
  vertical_geom%lnonhyd_geom = .false.
  vertical_geom%yrcver = cver

  do k = 0, nlev
    vab%vah(k) = ak(k + 1)
    vab%vbh(k) = bk(k + 1)
    vab%valh(k) = vab%vah(k) / vp00
  enddo
  vab%vaf(0) = vab%vah(0)
  vab%vbf(0) = vab%vbh(0)
  do k = 1, nlev
    vab%vc(k) = vab%vah(k) * vab%vbh(k - 1) - vab%vah(k - 1) * vab%vbh(k)
    vab%vaf(k) = 0.5_jprb * (vab%vah(k) + vab%vah(k - 1))
    vab%vbf(k) = 0.5_jprb * (vab%vbh(k) + vab%vbh(k - 1))
    vab%vdela(k) = vab%vah(k) - vab%vah(k - 1)
    vab%vdelb(k) = vab%vbh(k) - vab%vbh(k - 1)
  enddo
  vertical_geom%yrvab = vab

  do j = 1, ncol
    if (ps(j) <= 0.0_c_double) then
      ierr = 2_c_int
      return
    endif
    presh(j,nlev) = ps(j)
    geoph(j,nlev) = orography_geopotential(j)
    surface_temp(j) = temperature(j,nlev)
    do k = 1, nlev
      temp(j,k) = temperature(j,k)
      q(j,k) = specific_humidity(j,k)
    enddo
  enddo

  do lev = 1, nout
    if (height_levels(lev) < -10000.0_c_double) then
      ierr = 3_c_int
      return
    endif
    do j = 1, ncol
      target_geop(j,lev) = rg * height_levels(lev)
    enddo
  enddo

  call gphpre(ncol, nlev, 1, ncol, vab, cver, presh, pxyb=xyb, presf=presf)

  if (has_specific_humidity /= 0_c_int) then
    call gprcp(ncol, 1, ncol, nlev, pq=q, pr=rgas, pkap=kappa, ldthermact=.true.)
  else
    rgas(:,:) = rd
  endif

  call gpgeo(ncol, 1, ncol, nlev, geoph, geopf, temp, rgas, &
             xyb(:,:,yytxyb%m_lnpr), xyb(:,:,yytxyb%m_alph), vertical_geom)
  call fpps(ncol, 1, ncol, nlev, nout, geoph, target_geop, temp, rgas, presh, surface_temp, target_pressure)

  do lev = 1, nout
    do j = 1, ncol
      if (target_pressure(j,lev) <= 0.0_jprb) then
        ierr = 4_c_int
        return
      endif
      output(j,lev) = target_pressure(j,lev)
    enddo
  enddo
end subroutine fullpos_height_above_sea_pressures_c


subroutine fullpos_potential_vorticity_pressures_c(ncol, nlev, nout, pv_values, ak, bk, ps, coriolis, pv_levels, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev, nout
  real(c_double), intent(in) :: pv_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), coriolis(ncol), pv_levels(nout)
  real(c_double), intent(out) :: output(ncol, nout)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k, lev
  integer(kind=jpim), parameter :: kppm = 4
  integer(kind=jpim), allocatable :: klevb(:,:,:)
  logical, allocatable :: ld_belo(:,:), ld_bels(:,:), ld_blow(:), ld_bles(:)
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), prxp(:,:,:), prxpd(:,:,:)
  real(kind=jprb), allocatable :: pvf(:,:), prpres(:), prlnpres(:), pcp2d(:)

  interface
    subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
      real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
      real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    end subroutine ppinit
    subroutine ppltp(pvcap, kiterpv, ldisopv, kproma, kstart, kprof, kflev, prpresf, ptoup, pxtoup, ldcp, pcorio, prpres, prlnpres, pcp2d)
      use parkind1, only: jpim, jprb
      real(kind=jprb), intent(in) :: pvcap
      integer(kind=jpim), intent(in) :: kiterpv, kproma, kstart, kprof, kflev
      logical, intent(in) :: ldisopv, ldcp
      real(kind=jprb), intent(in) :: prpresf(kproma,kflev), ptoup(kproma,kflev), pxtoup, pcorio(kproma)
      real(kind=jprb), intent(out) :: prpres(kproma), prlnpres(kproma), pcp2d(kproma)
    end subroutine ppltp
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1 .or. nout <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(presh(ncol,0:nlev), presf(ncol,nlev), prxp(ncol,0:nlev,kppm), prxpd(ncol,0:nlev,kppm))
  allocate(klevb(ncol,nout,kppm), ld_belo(ncol,nout), ld_bels(ncol,nout), ld_blow(nout), ld_bles(nout))
  allocate(pvf(ncol,nlev), prpres(ncol), prlnpres(ncol), pcp2d(ncol))

  do j = 1, ncol
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) <= 0.0_jprb) presh(j,k) = 0.1_jprb
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) presf(j,k) = 0.1_jprb
      pvf(j,k) = pv_values(j,k)
    enddo
  enddo

  call ppinit(ncol, 1, ncol, nlev, kppm, presh, presf, prxp, prxpd)

  do lev = 1, nout
    if (pv_levels(lev) <= 0.0_c_double) then
      ierr = 2_c_int
      return
    endif
    call ppltp(1.0e-10_jprb, 1, .false., ncol, 1, ncol, nlev, presf, pvf, pv_levels(lev), .false., coriolis, &
               prpres, prlnpres, pcp2d)
    output(:, lev) = prpres(:)
  enddo
end subroutine fullpos_potential_vorticity_pressures_c


subroutine fullpos_gprcp_kappa_c(ncol, nlev, q_values, kappa_values, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: ncol, nlev
  real(c_double), intent(in) :: q_values(ncol, nlev)
  real(c_double), intent(out) :: kappa_values(ncol, nlev)
  integer(c_int), intent(out) :: ierr

  real(kind=jprb), allocatable :: q(:,:), kappa(:,:)

  interface
    subroutine gprcp(kproma, kstart, kprof, kflev, pq, pqi, pql, pqr, pqs, pqg, &
                     pcp, pr, pkap, pgfl, kgfltyp, ldthermact)
      use parkind1, only: jpim, jprb
      use yom_ygfl, only: ygfl
      integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
      real(kind=jprb), optional, intent(in) :: pq(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqi(kproma,kflev), pql(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqr(kproma,kflev), pqs(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pqg(kproma,kflev)
      real(kind=jprb), optional, intent(in) :: pgfl(kproma,kflev,ygfl%ndim)
      real(kind=jprb), optional, intent(out) :: pcp(kproma,kflev), pr(kproma,kflev)
      real(kind=jprb), optional, intent(out) :: pkap(kproma,kflev)
      integer(kind=jpim), optional, intent(in) :: kgfltyp
      logical, optional, intent(in) :: ldthermact
    end subroutine gprcp
  end interface

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 0) then
    ierr = 1_c_int
    return
  endif

  allocate(q(ncol,nlev), kappa(ncol,nlev))
  q = q_values
  call gprcp(ncol, 1, ncol, nlev, pq=q, pkap=kappa, ldthermact=.true.)
  kappa_values = kappa
end subroutine fullpos_gprcp_kappa_c


subroutine fullpos_diagnose_potential_vorticity_c(ncol, nlev, u_values, v_values, temperature, &
                                                  relative_vorticity, temperature_meridional_gradient, &
                                                  temperature_zonal_gradient, surface_pressure_meridional_gradient, &
                                                  surface_pressure_zonal_gradient, kappa_values, ak, bk, ps, coriolis, &
                                                  potential_vorticity, potential_temperature, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  use yomvert, only: tvab, vp00
  implicit none

  integer(c_int), intent(in) :: ncol, nlev
  real(c_double), intent(in) :: u_values(ncol, nlev), v_values(ncol, nlev), temperature(ncol, nlev)
  real(c_double), intent(in) :: relative_vorticity(ncol, nlev)
  real(c_double), intent(in) :: temperature_meridional_gradient(ncol, nlev), temperature_zonal_gradient(ncol, nlev)
  real(c_double), intent(in) :: surface_pressure_meridional_gradient(ncol), surface_pressure_zonal_gradient(ncol)
  real(c_double), intent(in) :: kappa_values(ncol, nlev)
  real(c_double), intent(in) :: ak(nlev + 1), bk(nlev + 1), ps(ncol), coriolis(ncol)
  real(c_double), intent(out) :: potential_vorticity(ncol, nlev), potential_temperature(ncol, nlev)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: j, k
  external :: gppvo
  type(tvab) :: vab
  real(kind=jprb), allocatable :: presh(:,:), presf(:,:), rdelp(:,:), kappa(:,:)
  real(kind=jprb), allocatable :: u(:,:), v(:,:), t(:,:), vort(:,:), tm(:,:), tl(:,:), spm(:), spl(:), cor(:)
  real(kind=jprb), allocatable :: pvo(:,:), theta(:,:)

  ierr = 0_c_int
  if (ncol <= 0 .or. nlev <= 1) then
    ierr = 1_c_int
    return
  endif

  allocate(vab%vah(0:nlev), vab%vbh(0:nlev), vab%valh(0:nlev), vab%vc(nlev), vab%vaf(0:nlev))
  allocate(vab%vbf(0:nlev), vab%vdela(nlev), vab%vdelb(nlev))
  allocate(presh(ncol,0:nlev), presf(ncol,nlev), rdelp(ncol,nlev), kappa(ncol,nlev))
  allocate(u(ncol,nlev), v(ncol,nlev), t(ncol,nlev), vort(ncol,nlev), tm(ncol,nlev), tl(ncol,nlev))
  allocate(spm(ncol), spl(ncol), cor(ncol), pvo(ncol,nlev), theta(ncol,nlev))

  vp00 = 100000.0_jprb
  do k = 0, nlev
    vab%vah(k) = ak(k + 1)
    vab%vbh(k) = bk(k + 1)
    vab%valh(k) = vab%vah(k) / vp00
  enddo
  do k = 1, nlev
    vab%vc(k) = vab%vah(k) * vab%vbh(k - 1) - vab%vah(k - 1) * vab%vbh(k)
    vab%vaf(k) = 0.5_jprb * (vab%vah(k) + vab%vah(k - 1))
    vab%vbf(k) = 0.5_jprb * (vab%vbh(k) + vab%vbh(k - 1))
    vab%vdela(k) = vab%vah(k) - vab%vah(k - 1)
    vab%vdelb(k) = vab%vbh(k) - vab%vbh(k - 1)
  enddo
  vab%vaf(0) = vab%vah(0)
  vab%vbf(0) = vab%vbh(0)

  do j = 1, ncol
    spm(j) = surface_pressure_meridional_gradient(j)
    spl(j) = surface_pressure_zonal_gradient(j)
    cor(j) = coriolis(j)
    do k = 0, nlev
      presh(j,k) = ak(k + 1) + bk(k + 1) * ps(j)
      if (presh(j,k) < 0.0_jprb) then
        ierr = 2_c_int
        return
      endif
      if (k > 0 .and. presh(j,k) <= 0.0_jprb) then
        ierr = 2_c_int
        return
      endif
    enddo
    do k = 1, nlev
      presf(j,k) = 0.5_jprb * (presh(j,k - 1) + presh(j,k))
      if (presf(j,k) <= 0.0_jprb) then
        ierr = 3_c_int
        return
      endif
      rdelp(j,k) = 1.0_jprb / abs(presh(j,k) - presh(j,k - 1))
      kappa(j,k) = kappa_values(j,k)
      u(j,k) = u_values(j,k)
      v(j,k) = v_values(j,k)
      t(j,k) = temperature(j,k)
      vort(j,k) = relative_vorticity(j,k)
      tm(j,k) = temperature_meridional_gradient(j,k)
      tl(j,k) = temperature_zonal_gradient(j,k)
    enddo
  enddo

  call gppvo(vab, ncol, 1, ncol, nlev, presf, rdelp, kappa, cor, vort, u, v, t, tm, tl, spm, spl, pvo, theta)

  potential_vorticity = pvo
  potential_temperature = theta
end subroutine fullpos_diagnose_potential_vorticity_c
