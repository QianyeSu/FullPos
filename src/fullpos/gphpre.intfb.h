interface
  subroutine gphpre(kproma, kflev, kstart, kprof, ydvab, ydcver, presh, pxyb, presf)
    use parkind1, only: jpim, jprb
    use yomcver, only: tcver
    use yomvert, only: tvab
    use intdyn_mod, only: yytxyb
    integer(kind=jpim), intent(in) :: kproma, kflev, kstart, kprof
    type(tvab), intent(in) :: ydvab
    type(tcver), intent(in) :: ydcver
    real(kind=jprb), intent(inout) :: presh(kproma,0:kflev)
    real(kind=jprb), optional, intent(out) :: pxyb(kproma,kflev,yytxyb%ndim)
    real(kind=jprb), optional, intent(out) :: presf(kproma,kflev)
  end subroutine gphpre
end interface
