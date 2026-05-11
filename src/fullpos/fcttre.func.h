interface
  real(kind=jprb) function foqs(es)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: es
  end function foqs
  real(kind=jprb) function fodqs(qs, es, desdt)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: qs, es, desdt
  end function fodqs
  real(kind=jprb) function fodlew(pt, zdelta)
    use parkind1, only: jprb
    real(kind=jprb), intent(in) :: pt, zdelta
  end function fodlew
end interface
