interface
  subroutine gpgeo(kproma, kstart, kprof, kflev, phi, phif, pt, pr, plnpr, palph, pvgeom)
    use parkind1, only: jpim, jprb
    use yomvert, only: tvertical_geom
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
    real(kind=jprb), intent(inout) :: phi(kproma,0:kflev)
    real(kind=jprb), intent(out) :: phif(kproma,kflev)
    real(kind=jprb), intent(in) :: pt(kproma,kflev), pr(kproma,kflev)
    real(kind=jprb), intent(in) :: plnpr(kproma,kflev), palph(kproma,kflev)
    type(tvertical_geom), intent(in) :: pvgeom
  end subroutine gpgeo
end interface
