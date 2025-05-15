module halo_models_mod
    use constants_mod
    use math_utils_mod
    implicit none
contains
    ! Modelo de halo de Burkert (corrigido)
    function burkert_halo(params, r) result(vh)
        real(kind=dp), intent(in) :: params(2)  ! params(1) = rs, params(2) = rho_s
        real(kind=dp), intent(in) :: r(:)
        real(kind=dp) :: vh(size(r))
        real(kind=dp) :: rs, rho_s, r_val, x, term1, term2, term3, numerator
        integer :: i, n
        
        rs = params(1)
        rho_s = params(2)
        n = size(r)
        
        do i = 1, n
            r_val = r(i)
            if (r_val > tiny_val) then
                x = r_val / rs
                term1 = log(1.0_dp + x)
                term2 = 0.5_dp * log(1.0_dp + x**2)
                term3 = atan(x)
                numerator = 4.0_dp * pi * G * rho_s * rs**3 * (term1 + term2 - term3)
                vh(i) = sqrt(max(numerator / r_val, 0.0_dp))  ! Garante valor não-negativo
            else
                vh(i) = 0.0_dp
            end if
        end do
    end function burkert_halo

    ! Modelo de halo NFW 
    function nfw_halo(params, r) result(vh)
        real(kind=dp), intent(in) :: params(2)  ! params(1) = rs, params(2) = rho_s
        real(kind=dp), intent(in) :: r(:)
        real(kind=dp) :: vh(size(r))
        real(kind=dp) :: rs, rho_s, r_val, x, term
        integer :: i, n
        
        rs = params(1)
        rho_s = params(2)
        n = size(r)
        
        do i = 1, n
            r_val = r(i)
            if (r_val > tiny_val) then
                x = r_val / rs
                term = log(1.0_dp + x) - x/(1.0_dp + x)
                vh(i) = sqrt(4.0_dp * pi * G * rho_s * rs**3 * term / r_val)
            else
                vh(i) = 0.0_dp
            end if
        end do
    end function nfw_halo
end module halo_models_mod
