interface
  real(kind=jprb) function foew(pt, zdelta)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: pt, zdelta
  end function foew
  real(kind=jprb) function foewm(pt)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: pt
  end function foewm
  real(kind=jprb) function foewmo(pt)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: pt
  end function foewmo
  real(kind=jprb) function foelson(pt)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: pt
  end function foelson
end interface
