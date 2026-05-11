interface
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
