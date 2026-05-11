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
  real(kind=jprb) :: ratm = 100000.0_jprb
  real(kind=jprb) :: rg = 9.80665_jprb
  real(kind=jprb) :: rd = 287.0596736665907_jprb
  real(kind=jprb) :: rv = 461.5249933083879_jprb
  real(kind=jprb) :: rcpd = 3.5_jprb * 287.0596736665907_jprb
  real(kind=jprb) :: rcvd = 2.5_jprb * 287.0596736665907_jprb
  real(kind=jprb) :: rcpv = 4.0_jprb * 461.5249933083879_jprb
  real(kind=jprb) :: rcw = 4218.0_jprb
  real(kind=jprb) :: rcs = 2106.0_jprb
end module yomcst

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
end module yomct0

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
end module yomvert

subroutine abor1(cdtext)
  implicit none
  character(len=*), intent(in) :: cdtext
  error stop cdtext
end subroutine abor1
