interface
  subroutine fpview(kproma, kst, knd, koplev, psp1, psp2, pres2, pdpblc, palfa)
    use parkind1, only: jpim, jprb
    integer(kind=jpim), intent(in) :: kproma, kst, knd, koplev
    real(kind=jprb), intent(in) :: psp1(kproma), psp2(kproma), pres2(kproma,koplev), pdpblc
    real(kind=jprb), intent(out) :: palfa(kproma,koplev)
  end subroutine fpview
end interface
