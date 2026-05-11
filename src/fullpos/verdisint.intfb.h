interface
  subroutine verdisint(ydvfe, ydcver, cdmethod, cdbc, kproma, kstart, kprof, kflev, pin, pout, kbc, kchunk)
    use parkind1, only: jpim, jprb
    use yomcver, only: tcver
    use yomvert, only: tvfe
    type(tvfe), intent(in) :: ydvfe
    type(tcver), intent(in) :: ydcver
    character(len=*), intent(in) :: cdmethod, cdbc
    integer(kind=jpim), intent(in) :: kproma, kstart, kprof, kflev
    real(kind=jprb), intent(in) :: pin(kproma,0:kflev+1)
    real(kind=jprb), intent(out) :: pout(kproma,0:kflev)
    integer(kind=jpim), optional, intent(in) :: kbc, kchunk
  end subroutine verdisint
end interface
