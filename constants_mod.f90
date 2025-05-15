module constants_mod
    implicit none
    ! Define precisão dupla
    integer, parameter :: dp = kind(1.0d0)
    
    ! Constante gravitacional em (km/s)^2 kpc / M_sun
    real(kind=dp), parameter :: G = 4.30091e-6_dp 
    
    ! Outras constantes úteis
    real(kind=dp), parameter :: pi = acos(-1.0_dp)
    real(kind=dp), parameter :: tiny_val = 1.0e-10_dp  ! Valor pequeno para evitar divisão por zero
end module constants_mod

