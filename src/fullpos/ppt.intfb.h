interface
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
    real(kind=jprb), optional, intent(in) :: psttf(kproma,0:kflev), pr2(kproma,kflev)
    real(kind=jprb), optional, intent(out) :: prtpp(kproma,klevp)
  end subroutine ppt
end interface
