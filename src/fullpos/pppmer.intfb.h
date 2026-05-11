interface
  subroutine pppmer(kproma, kstart, kprof, prpress, porog, ptstar, pt0, pmslppp, ldqnh)
    use parkind1, only: jpim, jprb
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof
    real(kind=jprb), intent(in) :: prpress(kproma), porog(kproma), ptstar(kproma), pt0(kproma)
    real(kind=jprb), intent(out) :: pmslppp(kproma)
    logical, optional, intent(in) :: ldqnh
  end subroutine pppmer
end interface
