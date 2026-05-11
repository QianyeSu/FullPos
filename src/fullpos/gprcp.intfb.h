interface
  subroutine gprcp(kproma, kstart, kprof, kflev, pq, pqi, pql, pqr, pqs, pqg, pcp, pr, pkap, pgfl, kgfltyp, ldthermact)
    use parkind1, only: jpim, jprb
    use yom_ygfl, only: ygfl
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
    real(kind=jprb), optional, intent(in) :: pq(kproma,kflev)
    real(kind=jprb), optional, intent(in) :: pqi(kproma,kflev), pql(kproma,kflev)
    real(kind=jprb), optional, intent(in) :: pqr(kproma,kflev), pqs(kproma,kflev)
    real(kind=jprb), optional, intent(in) :: pqg(kproma,kflev)
    real(kind=jprb), optional, intent(in) :: pgfl(kproma,kflev,ygfl%ndim)
    real(kind=jprb), optional, intent(out) :: pcp(kproma,kflev), pr(kproma,kflev), pkap(kproma,kflev)
    integer(kind=jpim), optional, intent(in) :: kgfltyp
    logical, optional, intent(in) :: ldthermact
  end subroutine gprcp
end interface
