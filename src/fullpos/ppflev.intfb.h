interface
  subroutine ppflev(kproma, kstart, kprof, kflev, klevp, kppm, prpres, prpresh, prpresf, &
                    klevb, ldbelo, ldbels, ldblow, ldbles)
    use parkind1, only: jpim, jprb
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev, klevp, kppm
    real(kind=jprb), intent(in) :: prpres(kproma,klevp), prpresh(kproma,0:kflev), prpresf(kproma,kflev)
    integer(kind=jpim), intent(inout) :: klevb(kproma,klevp,kppm)
    logical, intent(out) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp)
  end subroutine ppflev
end interface
