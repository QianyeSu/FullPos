interface
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
