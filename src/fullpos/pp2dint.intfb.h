interface
  subroutine pp2dint(kproma, kstart, kprof, kflev, kppm, klevp, klolev, khilev, klevb, kslct, &
                     ldbelo, ldblow, prpresc, prxp, prxpd, pfldi, pfbdo)
    use parkind1, only: jpim, jprb
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, kppm, klevp, klolev, khilev, kslct
    integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
    logical, intent(in) :: ldbelo(kproma,klevp), ldblow(klevp)
    real(kind=jprb), intent(in) :: prpresc(kproma,klevp), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    real(kind=jprb), intent(in) :: pfldi(kproma,0:kflev)
    real(kind=jprb), intent(out) :: pfbdo(kproma,klevp)
  end subroutine pp2dint
end interface
