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
