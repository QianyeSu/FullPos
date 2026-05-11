module yomhook
  use parkind1, only: jpim, jprb
  implicit none
  logical :: lhook = .false.
  integer(kind=jpim), parameter :: jphook = kind(1.0_jprb)
contains
  subroutine dr_hook(name, flag, handle)
    character(len=*), intent(in) :: name
    integer, intent(in) :: flag
    real(kind=jprb), intent(inout) :: handle
    handle = 0.0_jprb
  end subroutine dr_hook
end module yomhook

module yomcst
  use parkind1, only: jprb
  implicit none
  real(kind=jprb), parameter :: ratm = 100000.0_jprb
  real(kind=jprb), parameter :: rg = 9.80665_jprb
  real(kind=jprb), parameter :: rd = 287.0596736665907_jprb
  real(kind=jprb), parameter :: rv = 461.5249933083879_jprb
  real(kind=jprb), parameter :: rcpd = 3.5_jprb * 287.0596736665907_jprb
  real(kind=jprb), parameter :: rcvd = 2.5_jprb * 287.0596736665907_jprb
  real(kind=jprb), parameter :: rcpv = 4.0_jprb * 461.5249933083879_jprb
  real(kind=jprb), parameter :: rcw = 4218.0_jprb
  real(kind=jprb), parameter :: rcs = 2106.0_jprb
  real(kind=jprb), parameter :: retv = rv / rd - 1.0_jprb
  real(kind=jprb), parameter :: rtt = 273.16_jprb
  real(kind=jprb), parameter :: rlvtt = 2.5008e6_jprb
  real(kind=jprb), parameter :: rlstt = 2.8345e6_jprb
  real(kind=jprb), parameter :: rgamw = (rcw - rcpv) / rv
  real(kind=jprb), parameter :: rbetw = rlvtt / rv + rgamw * rtt
  real(kind=jprb), parameter :: ralpw = log(611.14_jprb) + rbetw / rtt + rgamw * log(rtt)
  real(kind=jprb), parameter :: rgams = (rcs - rcpv) / rv
  real(kind=jprb), parameter :: rbets = rlstt / rv + rgams * rtt
  real(kind=jprb), parameter :: ralps = log(611.14_jprb) + rbets / rtt + rgams * log(rtt)
  real(kind=jprb), parameter :: ralpd = ralps - ralpw
  real(kind=jprb), parameter :: rbetd = rbets - rbetw
  real(kind=jprb), parameter :: rgamd = rgams - rgamw
end module yomcst

module yomlun
  implicit none
  integer :: nulout = 6
  integer :: nulerr = 0
  integer :: nulnam = 0
end module yomlun

module yomphy
  use parkind1, only: jpim
  implicit none
  type tphy
    logical :: lneige = .false.
  end type tphy
  type(tphy) :: yrphy
end module yomphy

module yom_ygfl
  use parkind1, only: jpim, jprb
  implicit none
  type gfl_component
    logical :: lthermact = .false.
    integer(kind=jpim) :: mp = 1
    integer(kind=jpim) :: mp1 = 1
    integer(kind=jpim) :: mp5 = 1
    integer(kind=jpim) :: mp9_ph = 1
    integer(kind=jpim) :: mp_sl1 = 1
    real(kind=jprb) :: r = 0.0_jprb
    real(kind=jprb) :: rcp = 0.0_jprb
  end type gfl_component

  type type_ygfl
    integer(kind=jpim) :: ndim = 0
    integer(kind=jpim) :: numflds = 0
    type(gfl_component), allocatable :: ycomp(:)
    type(gfl_component) :: yg
    type(gfl_component) :: yi
    type(gfl_component) :: yl
    type(gfl_component) :: yq
    type(gfl_component) :: yr
    type(gfl_component) :: ys
  end type type_ygfl

  type(type_ygfl) :: ygfl
end module yom_ygfl

module yomsta
  use parkind1, only: jpim, jprb
  implicit none
  type tsta
    real(kind=jprb), allocatable :: stpreh(:)
    real(kind=jprb), allocatable :: stpre(:)
    real(kind=jprb), allocatable :: stphi(:)
    real(kind=jprb), allocatable :: sttem(:)
    real(kind=jprb), allocatable :: stden(:)
    real(kind=jprb), allocatable :: stz(:)
  end type tsta
  real(kind=jprb) :: rdtdz1 = -0.0065_jprb
  real(kind=jprb) :: rdtdz2 = 0.0_jprb
  real(kind=jprb) :: rdtdz3 = 0.001_jprb
  real(kind=jprb) :: rdtdz4 = 0.0028_jprb
  real(kind=jprb) :: rdtdz5 = 0.0_jprb
  real(kind=jprb) :: rdtdz6 = -0.0028_jprb
  real(kind=jprb) :: rdtdz7 = -0.002_jprb
  real(kind=jprb) :: rdtdz8 = 0.0_jprb
  real(kind=jprb) :: rdtdz9 = 0.0_jprb
  real(kind=jprb) :: rztrop = 11000.0_jprb
  real(kind=jprb) :: rzstra = 20000.0_jprb
  real(kind=jprb) :: rzstr2 = 32000.0_jprb
  real(kind=jprb) :: rzstpo = 47000.0_jprb
  real(kind=jprb) :: rzmeso = 51000.0_jprb
  real(kind=jprb) :: rzmes2 = 71000.0_jprb
  real(kind=jprb) :: rzmepo = 84852.0_jprb
  real(kind=jprb) :: rzabov = 90000.0_jprb
  real(kind=jprb) :: rtsur = 288.15_jprb
  real(kind=jprb) :: rttrop = 216.65_jprb
  real(kind=jprb) :: rtstra = 216.65_jprb
  real(kind=jprb) :: rtstr2 = 228.65_jprb
  real(kind=jprb) :: rtstpo = 270.65_jprb
  real(kind=jprb) :: rtmeso = 270.65_jprb
  real(kind=jprb) :: rtmes2 = 214.65_jprb
  real(kind=jprb) :: rtmepo = 186.946_jprb
  real(kind=jprb) :: rtabov = 186.946_jprb
  real(kind=jprb) :: rptrop = 22632.06_jprb
  real(kind=jprb) :: rpstra = 5474.889_jprb
  real(kind=jprb) :: rpstr2 = 868.0187_jprb
  real(kind=jprb) :: rpstpo = 110.9063_jprb
  real(kind=jprb) :: rpmeso = 66.93887_jprb
  real(kind=jprb) :: rpmes2 = 3.956420_jprb
  real(kind=jprb) :: rpmepo = 0.3734_jprb
  real(kind=jprb) :: rpabov = 0.1_jprb
  real(kind=jprb) :: vdtdz1 = -0.0065_jprb
  real(kind=jprb) :: vdtdz2 = 0.0_jprb
  real(kind=jprb) :: vdtdz3 = 0.001_jprb
  real(kind=jprb) :: vdtdz4 = 0.0028_jprb
  real(kind=jprb) :: vdtdz5 = 0.0_jprb
  real(kind=jprb) :: vdtdz6 = -0.0028_jprb
  real(kind=jprb) :: vdtdz7 = -0.002_jprb
  real(kind=jprb) :: vdtdz8 = 0.0_jprb
  real(kind=jprb) :: vdtdz9 = 0.0_jprb
  real(kind=jprb) :: vztrop = 11000.0_jprb
  real(kind=jprb) :: vzstra = 20000.0_jprb
  real(kind=jprb) :: vzstr2 = 32000.0_jprb
  real(kind=jprb) :: vzstpo = 47000.0_jprb
  real(kind=jprb) :: vzmeso = 51000.0_jprb
  real(kind=jprb) :: vzmes2 = 71000.0_jprb
  real(kind=jprb) :: vzmepo = 84852.0_jprb
  real(kind=jprb) :: vzabov = 90000.0_jprb
  real(kind=jprb) :: vtsur = 288.15_jprb
  real(kind=jprb) :: vttrop = 216.65_jprb
  real(kind=jprb) :: vtstra = 216.65_jprb
  real(kind=jprb) :: vtstr2 = 228.65_jprb
  real(kind=jprb) :: vtstpo = 270.65_jprb
  real(kind=jprb) :: vtmeso = 270.65_jprb
  real(kind=jprb) :: vtmes2 = 214.65_jprb
  real(kind=jprb) :: vtmepo = 186.946_jprb
  real(kind=jprb) :: vtabov = 186.946_jprb
  real(kind=jprb) :: vptrop = 22632.06_jprb
  real(kind=jprb) :: vpstra = 5474.889_jprb
  real(kind=jprb) :: vpstr2 = 868.0187_jprb
  real(kind=jprb) :: vpstpo = 110.9063_jprb
  real(kind=jprb) :: vpmeso = 66.93887_jprb
  real(kind=jprb) :: vpmes2 = 3.956420_jprb
  real(kind=jprb) :: vpmepo = 0.3734_jprb
  real(kind=jprb) :: vpabov = 0.1_jprb
  real(kind=jprb) :: hextrap = 0.0_jprb
  integer(kind=jpim) :: nlextrap = 1
end module yomsta

module yomct0
  implicit none
  logical :: loldpp = .false.
  logical :: lecmwf = .true.
end module yomct0

module type_gflflds
  use parkind1, only: jpim
  implicit none
  type type_igflfld
    integer(kind=jpim) :: iq = 1
    integer(kind=jpim) :: il = 2
    integer(kind=jpim) :: ii = 3
    integer(kind=jpim) :: irr = 4
    integer(kind=jpim) :: is = 5
    integer(kind=jpim) :: ig = 6
    integer(kind=jpim) :: ih = 7
  end type type_igflfld
end module type_gflflds

module yomcver
  use parkind1, only: jpim, jprb
  implicit none
  type tcver
    logical :: laprxpk = .true.
    integer(kind=jpim) :: ndlnpr = 0
    real(kind=jprb) :: rhydr0 = 0.6931471805599453_jprb
    logical :: lregeta = .false.
    logical :: lvfe_regeta = .false.
    integer(kind=jpim) :: nvsch = 0
    integer(kind=jpim) :: nvfe_type = 0
    integer(kind=jpim) :: nvfe_order = 0
    integer(kind=jpim) :: nvfe_intbc = 0
    integer(kind=jpim) :: nvfe_derbc = 0
    integer(kind=jpim) :: nvfe_internals = 0
    integer(kind=jpim) :: nvfe_bc = 0
    logical :: lvertfe = .false.
    logical :: lvfe_lapl = .false.
    logical :: lvfe_lapl_bc = .false.
    logical :: lvfe_lapl_tbc = .false.
    logical :: lvfe_lapl_bbc = .false.
    logical :: lvfe_lapl2pi = .false.
    real(kind=jprb) :: rlapl2pi = 0.0_jprb
    logical :: lvfe_x_term = .false.
    logical :: lvfe_z_term = .false.
    logical :: lvfe_gw = .false.
    logical :: lvfe_delnhpre = .false.
    logical :: lvfe_gwmpa = .false.
    logical :: lvfe_centri = .false.
    logical :: lvfe_cheb = .false.
    real(kind=jprb) :: rvfe_centri = 0.0_jprb
    real(kind=jprb) :: rvfe_alpha = 0.0_jprb
    real(kind=jprb) :: rvfe_beta = 0.0_jprb
    real(kind=jprb) :: rvfe_knot_stretch = 0.0_jprb
    logical :: lvfe_approx = .false.
    logical :: lvfe_ecmwf = .false.
    logical :: lvfe_lapl_half = .false.
    logical :: lvfe_fix_order = .false.
    logical :: lvfe_gw_half = .false.
    logical :: lvfe_maximas = .false.
    logical :: lvfe_verbose = .false.
    logical :: lvfe_normalize = .false.
    logical :: ldyn_analysis_stability = .false.
    real(kind=jprb) :: rmindeta = 0.0_jprb
    real(kind=jprb) :: rfac1 = 0.0_jprb
    real(kind=jprb) :: rfac2 = 0.0_jprb
  end type tcver
end module yomcver

module intdyn_mod
  use parkind1, only: jpim
  implicit none
  type txyb
    integer(kind=jpim) :: m_delp = 1
    integer(kind=jpim) :: m_rdelp = 2
    integer(kind=jpim) :: m_lnpr = 3
    integer(kind=jpim) :: m_alph = 4
    integer(kind=jpim) :: m_rtgr = 5
    integer(kind=jpim) :: m_rpre = 6
    integer(kind=jpim) :: m_rpp = 7
    integer(kind=jpim) :: ndim = 7
  end type txyb
  type(txyb) :: yytxyb
end module intdyn_mod

module yomppvi
  use parkind1, only: jprb
  implicit none
  logical :: lescale = .false.
  logical :: lescale_t = .false.
  logical :: lescale_q = .false.
  logical :: lescale_u = .false.
  logical :: lescale_pd = .false.
  logical :: lescale_gfl = .false.
  logical :: lrppuv_cstext = .false.
  logical :: lrppuv_callitpq = .true.
  logical :: lppvivx = .false.
  logical :: lnots_t = .false.
  real(kind=jprb) :: rppvivx = 0.0_jprb
  real(kind=jprb) :: rppvivp = 0.0_jprb
end module yomppvi

module yomvert
  use parkind1, only: jprb
  use yomcver, only: tcver
  implicit none
  real(kind=jprb) :: vp00 = 100000.0_jprb
  real(kind=jprb) :: toppres = 0.1_jprb
  type tvab
    real(kind=jprb), allocatable :: valh(:)
    real(kind=jprb), allocatable :: vbh(:)
    real(kind=jprb), allocatable :: vah(:)
    real(kind=jprb), allocatable :: vc(:)
    real(kind=jprb), allocatable :: vaf(:)
    real(kind=jprb), allocatable :: vbf(:)
    real(kind=jprb), allocatable :: vdela(:)
    real(kind=jprb), allocatable :: vdelb(:)
  end type tvab
  type tveta
    real(kind=jprb), allocatable :: vfe_rdetah(:)
  end type tveta
  type tvfe
    integer :: unused = 0
  end type tvfe
  type tvertical_geom
    logical :: lnonhyd_geom = .false.
    type(tvab) :: yrvab
    type(tveta) :: yrveta
    type(tvfe) :: yrvfe
    type(tcver) :: yrcver
  end type tvertical_geom
end module yomvert

subroutine verdisint(ydvfe, ydcver, cdmethod, cdbc, kproma, kstart, kprof, kflev, pin, pout, kbc, kchunk)
  use parkind1, only: jpim, jprb
  use yomcver, only: tcver
  use yomvert, only: tvfe
  implicit none
  type(tvfe), intent(in) :: ydvfe
  type(tcver), intent(in) :: ydcver
  character(len=*), intent(in) :: cdmethod, cdbc
  integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
  real(kind=jprb), intent(in) :: pin(kproma,0:kflev+1)
  real(kind=jprb), intent(out) :: pout(kproma,0:kflev)
  integer(kind=jpim), optional, intent(in) :: kbc, kchunk
  pout(kstart:kprof,0:kflev) = pin(kstart:kprof,0:kflev)
end subroutine verdisint

subroutine abor1(cdtext)
  implicit none
  character(len=*), intent(in) :: cdtext
  error stop cdtext
end subroutine abor1

subroutine gprh(ldwmorh, kproma, kstart, kprof, kflev, prhmax, prhmin, pq, pt, presf, pes, prh, ldsonntag)
  use parkind1, only: jpim, jprb
  use yomcst, only: retv, rtt
  use yomphy, only: yrphy
  use yomct0, only: lecmwf
  implicit none

  logical, intent(in) :: ldwmorh
  integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
  real(kind=jprb), intent(in) :: prhmax, prhmin
  real(kind=jprb), intent(in) :: pq(kproma,kflev), pt(kproma,kflev), presf(kproma,kflev)
  real(kind=jprb), intent(inout) :: pes(kproma,kflev)
  real(kind=jprb), intent(out) :: prh(kproma,kflev)
  logical, optional, intent(in) :: ldsonntag

  integer(kind=jpim) :: jl, jlev
  real(kind=jprb) :: zdelta, zes

  include "fcttre.func.h"
  include "fcttrm.func.h"

  do jlev = 1, kflev
    do jl = kstart, kprof
      if (lecmwf) then
        if (ldwmorh) then
          pes(jl,jlev) = foewmo(pt(jl,jlev))
        else
          pes(jl,jlev) = foewm(pt(jl,jlev))
        endif
      else
        if (yrphy%lneige) then
          zdelta = max(0.0_jprb, sign(1.0_jprb, rtt - pt(jl,jlev)))
        else
          zdelta = 0.0_jprb
        endif
        pes(jl,jlev) = foew(pt(jl,jlev), zdelta)
      endif
      zes = pes(jl,jlev)
      prh(jl,jlev) = max(prhmin, min((presf(jl,jlev) * pq(jl,jlev) * (retv + 1.0_jprb)) / &
        & ((1.0_jprb + retv * pq(jl,jlev)) * zes), prhmax))
    enddo
  enddo
end subroutine gprh

real(kind=jprb) function foew(pt, zdelta)
  use parkind1, only: jprb
  use yomcst, only: ralpw, ralpd, ralps, rbetw, rbetd, rbets, rgamw, rgamd, rgams, rtt
  implicit none
  real(kind=jprb), intent(in) :: pt, zdelta
  real(kind=jprb) :: ztmp

  ztmp = max(pt, 1.0_jprb)
  if (zdelta > 0.5_jprb) then
    foew = exp(ralps + rbets / ztmp + rgams * log(ztmp))
  else
    foew = exp(ralpw + rbetw / ztmp + rgamw * log(ztmp))
  endif
end function foew

real(kind=jprb) function foewm(pt)
  use parkind1, only: jprb
  use yomcst, only: ralpw, rbetw, rgamw
  implicit none
  real(kind=jprb), intent(in) :: pt
  foewm = exp(ralpw + rbetw / max(pt, 1.0_jprb) + rgamw * log(max(pt, 1.0_jprb)))
end function foewm

real(kind=jprb) function foewmo(pt)
  use parkind1, only: jprb
  use yomcst, only: ralpw, rbetw, rgamw
  implicit none
  real(kind=jprb), intent(in) :: pt
  foewmo = exp(ralpw + rbetw / max(pt, 1.0_jprb) + rgamw * log(max(pt, 1.0_jprb)))
end function foewmo

real(kind=jprb) function foelson(pt)
  use parkind1, only: jprb
  use yomcst, only: ralps, rbets, rgams
  implicit none
  real(kind=jprb), intent(in) :: pt
  foelson = exp(ralps + rbets / max(pt, 1.0_jprb) + rgams * log(max(pt, 1.0_jprb)))
end function foelson

real(kind=jprb) function foqs(es)
  use parkind1, only: jprb
  use yomcst, only: retv
  implicit none
  real(kind=jprb), intent(in) :: es
  foqs = es / max(1.0_jprb, 1.0_jprb + retv)
end function foqs

real(kind=jprb) function fodqs(qs, es, desdt)
  use parkind1, only: jprb
  implicit none
  real(kind=jprb), intent(in) :: qs, es, desdt
  fodqs = 0.0_jprb
end function fodqs

real(kind=jprb) function fodlew(pt, zdelta)
  use parkind1, only: jprb
  implicit none
  real(kind=jprb), intent(in) :: pt, zdelta
  fodlew = 0.0_jprb
end function fodlew
