interface
  subroutine ppgeop(kproma,kstart,kprof,kflev,klevp,klolev,kppm,klevb,pres,plnpres,ldbelo,ldbels,ldblow,ldbles,ldextr,porog,pvert_geom,&
                    prxp,prxpd,ptstar,pr0,plnpr,palph,ptf,pzpp,psttf,prred)
    use parkind1, only: jpim, jprb
    use yomvert, only: tvertical_geom
    integer(kind=jpim), intent(in) :: kproma,kstart,kprof,kflev,klevp,klolev,kppm
    integer(kind=jpim), intent(in) :: klevb(kproma,klevp,kppm)
    real(kind=jprb), intent(in) :: pres(kproma,klevp), plnpres(kproma,klevp)
    logical, intent(in) :: ldbelo(kproma,klevp), ldbels(kproma,klevp), ldblow(klevp), ldbles(klevp), ldextr
    real(kind=jprb), intent(in) :: porog(kproma), prxp(kproma,0:kflev,kppm), prxpd(kproma,0:kflev,kppm)
    real(kind=jprb), intent(in) :: ptstar(kproma), pr0(kproma,kflev), plnpr(kproma,kflev), palph(kproma,kflev), ptf(kproma,kflev)
    type(tvertical_geom), intent(in) :: pvert_geom
    real(kind=jprb), intent(inout) :: pzpp(kproma,klevp)
    real(kind=jprb), optional, intent(in) :: psttf(kproma,0:kflev), prred(kproma,kflev)
  end subroutine ppgeop
end interface
