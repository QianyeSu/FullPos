interface
  subroutine gprh(ldwmorh, kproma, kstart, kprof, kflev, prhmax, prhmin, pq, pt, presf, pes, prh, ldsonntag)
    use parkind1, only: jpim, jprb
    logical, intent(in) :: ldwmorh
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
    real(kind=jprb), intent(in) :: prhmax, prhmin
    real(kind=jprb), intent(in) :: pq(kproma,kflev), pt(kproma,kflev), presf(kproma,kflev)
    real(kind=jprb), intent(inout) :: pes(kproma,kflev)
    real(kind=jprb), intent(out) :: prh(kproma,kflev)
    logical, optional, intent(in) :: ldsonntag
  end subroutine gprh
end interface
