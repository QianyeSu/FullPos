interface
  subroutine ppinit(kproma, kstart, kprof, kflev, kppm, prpresh, prpresf, prxp, prxpd)
    use parkind1, only: jpim, jprb
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm
    real(kind=jprb), intent(in) :: prpresh(kproma,0:kflev), prpresf(kproma,kflev)
    real(kind=jprb), intent(out) :: prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
  end subroutine ppinit
end interface
