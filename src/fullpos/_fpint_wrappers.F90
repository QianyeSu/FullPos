subroutine fullpos_fpint4_c(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, &
                            kder, kdmp, ldml_i, ldnrst_i, kbox, kl0, pwxx, pwxy, &
                            ldmask_i, pbuf, prow, pundef) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf
  integer(c_int), intent(in) :: kder(kfields), kdmp(kfldbuf), ldnrst_i(kfields)
  integer(c_int), intent(in) :: kbox(kgpend), kl0(kfproma, kfprow), ldmask_i(kfields)
  integer(c_int), intent(in) :: ldml_i
  real(c_double), intent(in) :: pwxx(kfproma, 4), pwxy(kfproma, 12)
  real(c_double), intent(in) :: pbuf(kaslb1 * kfldbuf), pundef
  real(c_double), intent(out) :: prow(kfproma, kfields)

  logical :: ldml
  logical :: ldnrst(kfields), ldmask(kfields)

  interface
    subroutine fpint4(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, kder, &
                      kdmp, ldml, ldnrst, kbox, kl0, pwxx, pwxy, ldmask, pbuf, prow, pundef)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kaslb1, kfprow, kfields, kfproma, kfldbuf
      integer(kind=jpim), intent(in) :: kgpst, kgpend
      integer(kind=jpim), intent(in) :: kder(kfields), kdmp(kfldbuf)
      logical, intent(in) :: ldml
      logical, intent(in) :: ldnrst(kfields)
      integer(kind=jpim), intent(in) :: kbox(kgpend), kl0(kfproma, kfprow)
      real(kind=jprb), intent(in) :: pwxx(kfproma, 4), pwxy(kfproma, 5:16)
      logical, intent(in) :: ldmask(kfields)
      real(kind=jprb), intent(in) :: pbuf(kaslb1 * kfldbuf)
      real(kind=jprb), intent(out) :: prow(kfproma, kfields)
      real(kind=jprb), intent(in), optional :: pundef
    end subroutine fpint4
  end interface

  ldml = ldml_i /= 0
  ldnrst = ldnrst_i /= 0
  ldmask = ldmask_i /= 0
  call fpint4(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, kder, &
              kdmp, ldml, ldnrst, kbox, kl0, pwxx, pwxy, ldmask, pbuf, prow, pundef)
end subroutine fullpos_fpint4_c


subroutine fullpos_fpint12_c(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, &
                             kder, kdmp, ldml_i, ldnrst_i, ldmono_i, kbox, kl0, pwxx, &
                             pwxy, ldmask_i, pbuf, prow, pundef) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf
  integer(c_int), intent(in) :: kder(kfields), kdmp(kfldbuf), ldnrst_i(kfields), ldmono_i(kfields)
  integer(c_int), intent(in) :: kbox(kgpend), kl0(kfproma, kfprow), ldmask_i(kfields)
  integer(c_int), intent(in) :: ldml_i
  real(c_double), intent(in) :: pwxx(kfproma, 12), pwxy(kfproma, 12)
  real(c_double), intent(in) :: pbuf(kaslb1 * kfldbuf), pundef
  real(c_double), intent(out) :: prow(kfproma, kfields)

  logical :: ldml
  logical :: ldnrst(kfields), ldmono(kfields), ldmask(kfields)

  interface
    subroutine fpint12(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, kder, &
                       kdmp, ldml, ldnrst, ldmono, kbox, kl0, pwxx, pwxy, ldmask, &
                       pbuf, prow, pundef)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kaslb1, kfprow, kfields, kfproma, kfldbuf
      integer(kind=jpim), intent(in) :: kgpst, kgpend
      integer(kind=jpim), intent(in) :: kder(kfields), kdmp(kfldbuf)
      logical, intent(in) :: ldml
      logical, intent(in) :: ldnrst(kfields), ldmono(kfields)
      integer(kind=jpim), intent(in) :: kbox(kgpend), kl0(kfproma, kfprow)
      real(kind=jprb), intent(in) :: pwxx(kfproma, 12), pwxy(kfproma, 5:16)
      logical, intent(in) :: ldmask(kfields)
      real(kind=jprb), intent(in) :: pbuf(kaslb1 * kfldbuf)
      real(kind=jprb), intent(out) :: prow(kfproma, kfields)
      real(kind=jprb), intent(in), optional :: pundef
    end subroutine fpint12
  end interface

  ldml = ldml_i /= 0
  ldnrst = ldnrst_i /= 0
  ldmono = ldmono_i /= 0
  ldmask = ldmask_i /= 0
  call fpint12(kaslb1, kfprow, kfields, kgpst, kgpend, kfproma, kfldbuf, kder, &
               kdmp, ldml, ldnrst, ldmono, kbox, kl0, pwxx, pwxy, ldmask, &
               pbuf, prow, pundef)
end subroutine fullpos_fpint12_c


subroutine fullpos_fpavg_c(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, &
                           kfldbuf, kmask, ks0, ldmask_i, pbuf, pmask, pundef, prow) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf
  integer(c_int), intent(in) :: kmask(kfields), ks0(kfproma, kslwide * 2), ldmask_i(kfields)
  real(c_double), intent(in) :: pbuf(kaslb1 * kfldbuf), pmask(kfproma, kmasks), pundef
  real(c_double), intent(out) :: prow(kfproma, kfields)

  logical :: ldmask(kfields)

  interface
    subroutine fpavg(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf, &
                     kmask, ks0, ldmask, pbuf, pmask, pundef, prow)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kaslb1, kslwide, kfields, kmasks
      integer(kind=jpim), intent(in) :: kfproma, kfldbuf, kgpst, kgpend
      integer(kind=jpim), intent(in) :: kmask(kfields), ks0(kfproma, kslwide * 2)
      logical, intent(in) :: ldmask(kfields)
      real(kind=jprb), intent(in) :: pbuf(kaslb1 * kfldbuf), pmask(kfproma, kmasks), pundef
      real(kind=jprb), intent(out) :: prow(kfproma, kfields)
    end subroutine fpavg
  end interface

  ldmask = ldmask_i /= 0
  call fpavg(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf, &
             kmask, ks0, ldmask, pbuf, pmask, pundef, prow)
end subroutine fullpos_fpavg_c


subroutine fullpos_fpnear_c(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, &
                            kfldbuf, kmask, ks0, ldmask_i, pbuf, pmask, pundef, prow) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf
  integer(c_int), intent(in) :: kmask(kfields), ks0(kfproma, kslwide * 2), ldmask_i(kfields)
  real(c_double), intent(in) :: pbuf(kaslb1 * kfldbuf), pmask(kfproma, kmasks), pundef
  real(c_double), intent(out) :: prow(kfproma, kfields)

  logical :: ldmask(kfields)

  interface
    subroutine fpnear(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf, &
                      kmask, ks0, ldmask, pbuf, pmask, pundef, prow)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kaslb1, kslwide, kfields, kmasks
      integer(kind=jpim), intent(in) :: kfproma, kfldbuf, kgpst, kgpend
      integer(kind=jpim), intent(in) :: kmask(kfields), ks0(kfproma, kslwide * 2)
      logical, intent(in) :: ldmask(kfields)
      real(kind=jprb), intent(in) :: pbuf(kaslb1 * kfldbuf), pmask(kfproma, kmasks), pundef
      real(kind=jprb), intent(out) :: prow(kfproma, kfields)
    end subroutine fpnear
  end interface

  ldmask = ldmask_i /= 0
  call fpnear(kaslb1, kslwide, kfields, kmasks, kgpst, kgpend, kfproma, kfldbuf, &
              kmask, ks0, ldmask, pbuf, pmask, pundef, prow)
end subroutine fullpos_fpnear_c


subroutine fullpos_horizontal_halo_c(nlat, nsrc_points, ntarget_points, kslwide, use_near, &
                                    values, nloen, source_lats, target_lats, target_lons, &
                                    output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: nlat, nsrc_points, ntarget_points, kslwide, use_near
  integer(c_int), intent(in) :: nloen(nlat)
  real(c_double), intent(in) :: values(nsrc_points), source_lats(nlat)
  real(c_double), intent(in) :: target_lats(ntarget_points), target_lons(ntarget_points)
  real(c_double), intent(out) :: output(ntarget_points)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: kdgsa, kdgen, kfrstloff, kproma, kend, kfields, kmasks
  integer(kind=jpim) :: kfldbuf, kgpst, kgpend, kaslb1, kslw
  integer(kind=jpim) :: j, slot, col, n, total_ext, offset_src, offset_ext, row, lon0
  integer(kind=jpim), allocatable :: kla(:), klo(:,:), ks0(:,:), kmask(:)
  integer(kind=jpim), allocatable :: ldmask_i(:)
  logical :: ldmask(1)
  real(kind=jprb) :: pi, p4jp, pundef
  real(kind=jprb), allocatable :: pbuf_ext(:), pmask(:,:), prow(:,:)

  interface
    subroutine suhox1(kslwide, kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, &
                      pi, kloen, platin, plat, plon, kla, klo)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kslwide, kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: p4jp, pi
      integer(kind=jpim), intent(in) :: kloen(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: plat(kend), plon(kend)
      integer(kind=jpim), intent(out) :: kla(kproma)
      integer(kind=jpim), intent(out) :: klo(kproma,2*kslwide)
    end subroutine suhox1
  end interface

  ierr = 0_c_int
  if (kslwide <= 0_c_int) then
    ierr = 1_c_int
    return
  endif
  if (sum(nloen) /= nsrc_points) then
    ierr = 2_c_int
    return
  endif

  kdgsa = 1
  kdgen = nlat
  kfrstloff = 0
  kproma = ntarget_points
  kend = ntarget_points
  kslw = kslwide
  pi = acos(-1.0_jprb)
  p4jp = real(nlat, jprb) / pi
  pundef = 1.0e20_jprb

  allocate(kla(kproma), klo(kproma,2*kslw))
  call suhox1(kslw, kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, nloen, &
              source_lats, target_lats, target_lons, kla, klo)

  total_ext = 0
  do j = 1, nlat
    total_ext = total_ext + nloen(j) + 2 * kslw
  enddo
  allocate(pbuf_ext(total_ext))
  offset_src = 0
  offset_ext = 0
  do j = 1, nlat
    n = nloen(j)
    do col = 1, kslw
      pbuf_ext(offset_ext + col) = values(offset_src + mod(n - kslw + col - 1, n) + 1)
    enddo
    do col = 1, n
      pbuf_ext(offset_ext + kslw + col) = values(offset_src + col)
    enddo
    do col = 1, kslw
      pbuf_ext(offset_ext + kslw + n + col) = values(offset_src + mod(col - 1, n) + 1)
    enddo
    offset_src = offset_src + n
    offset_ext = offset_ext + n + 2 * kslw
  enddo

  allocate(ks0(kproma,2*kslw))
  do j = 1, kproma
    do slot = 1, 2 * kslw
      row = kla(j) + slot - kslw
      if (row < 1) row = 1
      if (row > nlat) row = nlat
      lon0 = mod(klo(j,slot), nloen(row))
      ks0(j,slot) = row_offset_halo(nloen, nlat, row, kslw) + kslw + lon0 + 1
    enddo
  enddo

  kfields = 1
  kmasks = 1
  kfldbuf = 1
  kgpst = 1
  kgpend = kproma
  kaslb1 = total_ext
  allocate(kmask(1), ldmask_i(1), pmask(kproma,1), prow(kproma,1))
  kmask = 1
  ldmask_i = 0
  ldmask = .false.
  pmask = 1.0_jprb

  if (use_near /= 0) then
    call fpnear(kaslb1, kslw, kfields, kmasks, kgpst, kgpend, kproma, kfldbuf, &
                kmask, ks0, ldmask, pbuf_ext, pmask, pundef, prow)
  else
    call fpavg(kaslb1, kslw, kfields, kmasks, kgpst, kgpend, kproma, kfldbuf, &
               kmask, ks0, ldmask, pbuf_ext, pmask, pundef, prow)
  endif
  output = prow(:,1)

contains
  integer(kind=jpim) function row_offset_halo(kloen_local, nlat_local, row, halo)
    integer(c_int), intent(in) :: kloen_local(nlat_local)
    integer(c_int), intent(in) :: nlat_local
    integer(kind=jpim), intent(in) :: row, halo
    integer(kind=jpim) :: jj
    row_offset_halo = 0
    do jj = 1, row - 1
      row_offset_halo = row_offset_halo + kloen_local(jj) + 2 * halo
    enddo
  end function row_offset_halo
end subroutine fullpos_horizontal_halo_c


subroutine fullpos_horizontal_halo_batch_c(nlat, nsrc_points, nfields, ntarget_points, kslwide, use_near, &
                                          values, nloen, source_lats, target_lats, target_lons, &
                                          output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: nlat, nsrc_points, nfields, ntarget_points, kslwide, use_near
  integer(c_int), intent(in) :: nloen(nlat)
  real(c_double), intent(in) :: values(nsrc_points,nfields), source_lats(nlat)
  real(c_double), intent(in) :: target_lats(ntarget_points), target_lons(ntarget_points)
  real(c_double), intent(out) :: output(ntarget_points,nfields)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: kdgsa, kdgen, kfrstloff, kproma, kend, kfields, kmasks
  integer(kind=jpim) :: kfldbuf, kgpst, kgpend, kaslb1, kslw
  integer(kind=jpim) :: fld, j, slot, col, n, total_ext, offset_src, offset_ext, row, lon0
  integer(kind=jpim), allocatable :: kla(:), klo(:,:), ks0(:,:), kmask(:)
  logical, allocatable :: ldmask(:)
  real(kind=jprb) :: pi, p4jp, pundef
  real(kind=jprb), allocatable :: pbuf_ext(:), pmask(:,:), prow(:,:)

  interface
    subroutine suhox1(kslwide, kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, &
                      pi, kloen, platin, plat, plon, kla, klo)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kslwide, kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: p4jp, pi
      integer(kind=jpim), intent(in) :: kloen(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: plat(kend), plon(kend)
      integer(kind=jpim), intent(out) :: kla(kproma)
      integer(kind=jpim), intent(out) :: klo(kproma,2*kslwide)
    end subroutine suhox1
  end interface

  ierr = 0_c_int
  if (kslwide <= 0_c_int .or. nfields <= 0_c_int) then
    ierr = 1_c_int
    return
  endif
  if (sum(nloen) /= nsrc_points) then
    ierr = 2_c_int
    return
  endif

  kdgsa = 1
  kdgen = nlat
  kfrstloff = 0
  kproma = ntarget_points
  kend = ntarget_points
  kslw = kslwide
  pi = acos(-1.0_jprb)
  p4jp = real(nlat, jprb) / pi
  pundef = 1.0e20_jprb

  allocate(kla(kproma), klo(kproma,2*kslw))
  call suhox1(kslw, kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, nloen, &
              source_lats, target_lats, target_lons, kla, klo)

  total_ext = 0
  do j = 1, nlat
    total_ext = total_ext + nloen(j) + 2 * kslw
  enddo
  kfields = nfields
  kfldbuf = kfields
  allocate(pbuf_ext(total_ext * kfields))
  do fld = 1, kfields
    offset_src = 0
    offset_ext = (fld - 1) * total_ext
    do j = 1, nlat
      n = nloen(j)
      do col = 1, kslw
        pbuf_ext(offset_ext + col) = values(offset_src + mod(n - kslw + col - 1, n) + 1, fld)
      enddo
      do col = 1, n
        pbuf_ext(offset_ext + kslw + col) = values(offset_src + col, fld)
      enddo
      do col = 1, kslw
        pbuf_ext(offset_ext + kslw + n + col) = values(offset_src + mod(col - 1, n) + 1, fld)
      enddo
      offset_src = offset_src + n
      offset_ext = offset_ext + n + 2 * kslw
    enddo
  enddo

  allocate(ks0(kproma,2*kslw))
  do j = 1, kproma
    do slot = 1, 2 * kslw
      row = kla(j) + slot - kslw
      if (row < 1) row = 1
      if (row > nlat) row = nlat
      lon0 = mod(klo(j,slot), nloen(row))
      ks0(j,slot) = row_offset_halo(nloen, nlat, row, kslw) + kslw + lon0 + 1
    enddo
  enddo

  kmasks = 1
  kgpst = 1
  kgpend = kproma
  kaslb1 = total_ext
  allocate(kmask(kfields), ldmask(kfields), pmask(kproma,1), prow(kproma,kfields))
  kmask = 1
  ldmask = .false.
  pmask = 1.0_jprb

  if (use_near /= 0) then
    call fpnear(kaslb1, kslw, kfields, kmasks, kgpst, kgpend, kproma, kfldbuf, &
                kmask, ks0, ldmask, pbuf_ext, pmask, pundef, prow)
  else
    call fpavg(kaslb1, kslw, kfields, kmasks, kgpst, kgpend, kproma, kfldbuf, &
               kmask, ks0, ldmask, pbuf_ext, pmask, pundef, prow)
  endif
  output = prow

contains
  integer(kind=jpim) function row_offset_halo(kloen_local, nlat_local, row, halo)
    integer(c_int), intent(in) :: kloen_local(nlat_local)
    integer(c_int), intent(in) :: nlat_local
    integer(kind=jpim), intent(in) :: row, halo
    integer(kind=jpim) :: jj
    row_offset_halo = 0
    do jj = 1, row - 1
      row_offset_halo = row_offset_halo + kloen_local(jj) + 2 * halo
    enddo
  end function row_offset_halo
end subroutine fullpos_horizontal_halo_batch_c


subroutine fullpos_horizontal_regular_c(nlat, nsrc_points, ntarget_points, kbinl, ldmono_i, values, nloen, &
                                        source_lats, target_lats, target_lons, output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: nlat, nsrc_points, ntarget_points, kbinl, ldmono_i
  integer(c_int), intent(in) :: nloen(nlat)
  real(c_double), intent(in) :: values(nsrc_points), source_lats(nlat)
  real(c_double), intent(in) :: target_lats(ntarget_points), target_lons(ntarget_points)
  real(c_double), intent(out) :: output(ntarget_points)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: kdgsa, kdgen, kfrstloff, kproma, kend, kfprow, kbinl_suhow1
  integer(kind=jpim) :: j, col, n, total_ext, offset_src, offset_ext
  integer(kind=jpim), allocatable :: kla(:), klon(:), klonn(:), klos(:), kloss(:)
  integer(kind=jpim), allocatable :: kl0(:,:), kbox(:), kder(:), kdmp(:)
  integer(kind=jpim) :: kgpst, kgpend, kfields, kfldbuf, kaslb1
  logical :: ldml
  logical :: ldnrst(1), ldmask(1), ldmono(1)
  real(kind=jprb) :: pi, p4jp, pundef
  real(kind=jprb) :: pcr(2)
  real(kind=jprb), allocatable :: pdlat(:,:), pdlo(:,:,:), pwxy(:,:), pwxx4(:,:), pwxx12(:,:)
  real(kind=jprb), allocatable :: prgmsd(:), pgmsf(:), pbuf_ext(:), prow(:,:)

  interface
    subroutine suhow1(kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, kloen, platin, &
                      kbinl, ldml, plat, plon, kla, klon, klonn, klos, kloss, pdlat, pdlo)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: p4jp, pi
      integer(kind=jpim), intent(in) :: kloen(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      integer(kind=jpim), intent(in) :: kbinl
      logical, intent(in) :: ldml
      real(kind=jprb), intent(in) :: plat(kend), plon(kend)
      integer(kind=jpim), intent(out) :: kla(kproma), klon(kproma), klonn(kproma)
      integer(kind=jpim), intent(out) :: klos(kproma), kloss(kproma)
      real(kind=jprb), intent(out) :: pdlat(4,kproma), pdlo(0:3,4,kproma)
    end subroutine suhow1

    subroutine suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, platin, prgmsd, pgmsf, &
                      kbinl, ldml, kla, pdlat, pdlo, pwxx, pwxy)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: pcr(2)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: prgmsd(kend), pgmsf(kend)
      integer(kind=jpim), intent(in) :: kbinl
      logical, intent(in) :: ldml
      integer(kind=jpim), intent(in) :: kla(kend)
      real(kind=jprb), intent(in) :: pdlat(4,kproma), pdlo(0:3,4,kproma)
      real(kind=jprb), intent(out) :: pwxx(kproma,kbinl), pwxy(kproma,5:16)
    end subroutine suhow2
  end interface

  ierr = 0_c_int
  if (kbinl /= 4_c_int .and. kbinl /= 12_c_int) then
    ierr = 1_c_int
    return
  endif
  if (sum(nloen) /= nsrc_points) then
    ierr = 2_c_int
    return
  endif

  kdgsa = 1
  kdgen = nlat
  kfrstloff = 0
  kproma = ntarget_points
  kend = ntarget_points
  kfprow = 4
  pi = acos(-1.0_jprb)
  p4jp = real(nlat, jprb) / pi
  ldml = .false.
  pcr = 0.0_jprb
  pundef = 1.0e20_jprb

  allocate(kla(kproma), klon(kproma), klonn(kproma), klos(kproma), kloss(kproma))
  allocate(pdlat(4,kproma), pdlo(0:3,4,kproma))
  allocate(prgmsd(kproma), pgmsf(kproma))
  prgmsd = 1.0_jprb
  pgmsf = 1.0_jprb

  ! SUHOW1 uses the historical selector 0 for 12-point interpolation,
  ! while SUHOW2 and FPINT12 use the explicit stencil size 12.
  if (kbinl == 12) then
    kbinl_suhow1 = 0
  else
    kbinl_suhow1 = kbinl
  endif

  call suhow1(kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, nloen, source_lats, &
              kbinl_suhow1, ldml, target_lats, target_lons, kla, klon, klonn, klos, kloss, &
              pdlat, pdlo)

  if (minval(kla) < 2 .or. maxval(kla) > nlat - 2) then
    ierr = 3_c_int
    return
  endif

  allocate(kl0(kproma,kfprow))
  total_ext = 0
  do j = 1, nlat
    total_ext = total_ext + nloen(j) + 3
  enddo
  allocate(pbuf_ext(total_ext))
  offset_src = 0
  offset_ext = 0
  do j = 1, nlat
    n = nloen(j)
    ! One western halo point, the native row, and two eastern halo points.
    ! FPINT12 reads KL0, KL0+1, KL0+2, and KL0+3 on each row.
    pbuf_ext(offset_ext + 1) = values(offset_src + n)
    do col = 1, n
      pbuf_ext(offset_ext + col + 1) = values(offset_src + col)
    enddo
    pbuf_ext(offset_ext + n + 2) = values(offset_src + 1)
    pbuf_ext(offset_ext + n + 3) = values(offset_src + min(2, n))
    offset_src = offset_src + n
    offset_ext = offset_ext + n + 3
  enddo

  do j = 1, kproma
    kl0(j,1) = row_offset_ext(nloen, nlat, kla(j)-1) + mod(klonn(j), nloen(kla(j)-1)) + 1
    kl0(j,2) = row_offset_ext(nloen, nlat, kla(j)  ) + mod(klon(j),  nloen(kla(j)))   + 1
    kl0(j,3) = row_offset_ext(nloen, nlat, kla(j)+1) + mod(klos(j),  nloen(kla(j)+1)) + 1
    kl0(j,4) = row_offset_ext(nloen, nlat, kla(j)+2) + mod(kloss(j), nloen(kla(j)+2)) + 1
  enddo

  allocate(pwxy(kproma,12))
  if (kbinl == 4) then
    allocate(pwxx4(kproma,4))
    call suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, source_lats, prgmsd, pgmsf, &
                kbinl, ldml, kla, pdlat, pdlo, pwxx4, pwxy)
  else
    allocate(pwxx12(kproma,12))
    call suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, source_lats, prgmsd, pgmsf, &
                kbinl, ldml, kla, pdlat, pdlo, pwxx12, pwxy)
  endif

  kfields = 1
  kfldbuf = 1
  kgpst = 1
  kgpend = kproma
  kaslb1 = total_ext
  allocate(kbox(kproma), kder(1), kdmp(1), prow(kproma,1))
  kbox = 1
  kder = 1
  kdmp = 0
  ldnrst = .false.
  ldmask = .false.
  ldmono = ldmono_i /= 0

  if (kbinl == 4) then
    call fpint4(kaslb1, kfprow, kfields, kgpst, kgpend, kproma, kfldbuf, kder, kdmp, &
                ldml, ldnrst, kbox, kl0, pwxx4, pwxy, ldmask, pbuf_ext, prow, pundef)
  else
    call fpint12(kaslb1, kfprow, kfields, kgpst, kgpend, kproma, kfldbuf, kder, kdmp, &
                 ldml, ldnrst, ldmono, kbox, kl0, pwxx12, pwxy, ldmask, pbuf_ext, prow, pundef)
  endif
  output = prow(:,1)

contains
  integer(kind=jpim) function row_offset_ext(kloen_local, nlat_local, row)
    integer(c_int), intent(in) :: kloen_local(nlat_local)
    integer(c_int), intent(in) :: nlat_local
    integer(kind=jpim), intent(in) :: row
    integer(kind=jpim) :: jj
    row_offset_ext = 0
    do jj = 1, row - 1
      row_offset_ext = row_offset_ext + kloen_local(jj) + 3
    enddo
  end function row_offset_ext
end subroutine fullpos_horizontal_regular_c


subroutine fullpos_horizontal_regular_batch_c(nlat, nsrc_points, nfields, ntarget_points, kbinl, ldmono_i, &
                                             values, nloen, source_lats, target_lats, target_lons, &
                                             output, ierr) bind(C)
  use iso_c_binding, only: c_int, c_double
  use parkind1, only: jpim, jprb
  implicit none

  integer(c_int), intent(in) :: nlat, nsrc_points, nfields, ntarget_points, kbinl, ldmono_i
  integer(c_int), intent(in) :: nloen(nlat)
  real(c_double), intent(in) :: values(nsrc_points,nfields), source_lats(nlat)
  real(c_double), intent(in) :: target_lats(ntarget_points), target_lons(ntarget_points)
  real(c_double), intent(out) :: output(ntarget_points,nfields)
  integer(c_int), intent(out) :: ierr

  integer(kind=jpim) :: kdgsa, kdgen, kfrstloff, kproma, kend, kfprow, kbinl_suhow1
  integer(kind=jpim) :: fld, j, col, n, total_ext, offset_src, offset_ext
  integer(kind=jpim), allocatable :: kla(:), klon(:), klonn(:), klos(:), kloss(:)
  integer(kind=jpim), allocatable :: kl0(:,:), kbox(:), kder(:), kdmp(:)
  integer(kind=jpim) :: kgpst, kgpend, kfields, kfldbuf, kaslb1
  logical :: ldml
  logical, allocatable :: ldnrst(:), ldmask(:), ldmono(:)
  real(kind=jprb) :: pi, p4jp, pundef
  real(kind=jprb) :: pcr(2)
  real(kind=jprb), allocatable :: pdlat(:,:), pdlo(:,:,:), pwxy(:,:), pwxx4(:,:), pwxx12(:,:)
  real(kind=jprb), allocatable :: prgmsd(:), pgmsf(:), pbuf_ext(:), prow(:,:)

  interface
    subroutine suhow1(kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, kloen, platin, &
                      kbinl, ldml, plat, plon, kla, klon, klonn, klos, kloss, pdlat, pdlo)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: p4jp, pi
      integer(kind=jpim), intent(in) :: kloen(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      integer(kind=jpim), intent(in) :: kbinl
      logical, intent(in) :: ldml
      real(kind=jprb), intent(in) :: plat(kend), plon(kend)
      integer(kind=jpim), intent(out) :: kla(kproma), klon(kproma), klonn(kproma)
      integer(kind=jpim), intent(out) :: klos(kproma), kloss(kproma)
      real(kind=jprb), intent(out) :: pdlat(4,kproma), pdlo(0:3,4,kproma)
    end subroutine suhow1

    subroutine suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, platin, prgmsd, pgmsf, &
                      kbinl, ldml, kla, pdlat, pdlo, pwxx, pwxy)
      use parkind1, only: jpim, jprb
      integer(kind=jpim), intent(in) :: kdgsa, kdgen, kfrstloff, kproma, kend
      real(kind=jprb), intent(in) :: pcr(2)
      real(kind=jprb), intent(in) :: platin(kdgsa+kfrstloff:kdgen+kfrstloff)
      real(kind=jprb), intent(in) :: prgmsd(kend), pgmsf(kend)
      integer(kind=jpim), intent(in) :: kbinl
      logical, intent(in) :: ldml
      integer(kind=jpim), intent(in) :: kla(kend)
      real(kind=jprb), intent(in) :: pdlat(4,kproma), pdlo(0:3,4,kproma)
      real(kind=jprb), intent(out) :: pwxx(kproma,kbinl), pwxy(kproma,5:16)
    end subroutine suhow2
  end interface

  ierr = 0_c_int
  if ((kbinl /= 4_c_int .and. kbinl /= 12_c_int) .or. nfields <= 0_c_int) then
    ierr = 1_c_int
    return
  endif
  if (sum(nloen) /= nsrc_points) then
    ierr = 2_c_int
    return
  endif

  kdgsa = 1
  kdgen = nlat
  kfrstloff = 0
  kproma = ntarget_points
  kend = ntarget_points
  kfprow = 4
  pi = acos(-1.0_jprb)
  p4jp = real(nlat, jprb) / pi
  ldml = .false.
  pcr = 0.0_jprb
  pundef = 1.0e20_jprb

  allocate(kla(kproma), klon(kproma), klonn(kproma), klos(kproma), kloss(kproma))
  allocate(pdlat(4,kproma), pdlo(0:3,4,kproma))
  allocate(prgmsd(kproma), pgmsf(kproma))
  prgmsd = 1.0_jprb
  pgmsf = 1.0_jprb

  if (kbinl == 12) then
    kbinl_suhow1 = 0
  else
    kbinl_suhow1 = kbinl
  endif

  call suhow1(kdgsa, kdgen, kfrstloff, kproma, kend, p4jp, pi, nloen, source_lats, &
              kbinl_suhow1, ldml, target_lats, target_lons, kla, klon, klonn, klos, kloss, &
              pdlat, pdlo)

  if (minval(kla) < 2 .or. maxval(kla) > nlat - 2) then
    ierr = 3_c_int
    return
  endif

  allocate(kl0(kproma,kfprow))
  total_ext = 0
  do j = 1, nlat
    total_ext = total_ext + nloen(j) + 3
  enddo
  kfields = nfields
  kfldbuf = kfields
  allocate(pbuf_ext(total_ext * kfields))
  do fld = 1, kfields
    offset_src = 0
    offset_ext = (fld - 1) * total_ext
    do j = 1, nlat
      n = nloen(j)
      pbuf_ext(offset_ext + 1) = values(offset_src + n, fld)
      do col = 1, n
        pbuf_ext(offset_ext + col + 1) = values(offset_src + col, fld)
      enddo
      pbuf_ext(offset_ext + n + 2) = values(offset_src + 1, fld)
      pbuf_ext(offset_ext + n + 3) = values(offset_src + min(2, n), fld)
      offset_src = offset_src + n
      offset_ext = offset_ext + n + 3
    enddo
  enddo

  do j = 1, kproma
    kl0(j,1) = row_offset_ext(nloen, nlat, kla(j)-1) + mod(klonn(j), nloen(kla(j)-1)) + 1
    kl0(j,2) = row_offset_ext(nloen, nlat, kla(j)  ) + mod(klon(j),  nloen(kla(j)))   + 1
    kl0(j,3) = row_offset_ext(nloen, nlat, kla(j)+1) + mod(klos(j),  nloen(kla(j)+1)) + 1
    kl0(j,4) = row_offset_ext(nloen, nlat, kla(j)+2) + mod(kloss(j), nloen(kla(j)+2)) + 1
  enddo

  allocate(pwxy(kproma,12))
  if (kbinl == 4) then
    allocate(pwxx4(kproma,4))
    call suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, source_lats, prgmsd, pgmsf, &
                kbinl, ldml, kla, pdlat, pdlo, pwxx4, pwxy)
  else
    allocate(pwxx12(kproma,12))
    call suhow2(kdgsa, kdgen, kfrstloff, kproma, kend, pcr, source_lats, prgmsd, pgmsf, &
                kbinl, ldml, kla, pdlat, pdlo, pwxx12, pwxy)
  endif

  kgpst = 1
  kgpend = kproma
  kaslb1 = total_ext
  allocate(kbox(kproma), kder(kfields), kdmp(kfields), prow(kproma,kfields))
  allocate(ldnrst(kfields), ldmask(kfields), ldmono(kfields))
  kbox = 1
  kder = 1
  kdmp = 0
  ldnrst = .false.
  ldmask = .false.
  ldmono = ldmono_i /= 0

  if (kbinl == 4) then
    call fpint4(kaslb1, kfprow, kfields, kgpst, kgpend, kproma, kfldbuf, kder, kdmp, &
                ldml, ldnrst, kbox, kl0, pwxx4, pwxy, ldmask, pbuf_ext, prow, pundef)
  else
    call fpint12(kaslb1, kfprow, kfields, kgpst, kgpend, kproma, kfldbuf, kder, kdmp, &
                 ldml, ldnrst, ldmono, kbox, kl0, pwxx12, pwxy, ldmask, pbuf_ext, prow, pundef)
  endif
  output = prow

contains
  integer(kind=jpim) function row_offset_ext(kloen_local, nlat_local, row)
    integer(c_int), intent(in) :: kloen_local(nlat_local)
    integer(c_int), intent(in) :: nlat_local
    integer(kind=jpim), intent(in) :: row
    integer(kind=jpim) :: jj
    row_offset_ext = 0
    do jj = 1, row - 1
      row_offset_ext = row_offset_ext + kloen_local(jj) + 3
    enddo
  end function row_offset_ext
end subroutine fullpos_horizontal_regular_batch_c
