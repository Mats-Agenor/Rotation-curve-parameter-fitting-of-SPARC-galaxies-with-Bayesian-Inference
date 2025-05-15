module stats_mod
    use constants_mod
    use math_utils_mod
    use halo_models_mod
    implicit none
    
    ! Interface para função IEEE_IS_NAN se disponível
    interface
        elemental function ieee_is_nan(x) result(y)
            use, intrinsic :: ieee_arithmetic
            real(kind=8), intent(in) :: x
            logical :: y
        end function ieee_is_nan
    end interface

contains
    ! Função para verificar NaN de forma portável
    function has_nan(array) result(res)
        real(kind=dp), intent(in) :: array(:)
        logical :: res
        integer :: i
        
        res = .false.
        do i = 1, size(array)
            ! Verificação de NaN portável
            if (array(i) /= array(i)) then
                res = .true.
                exit
            end if
        end do
    end function has_nan

    ! Função chi-quadrado com termo de prior e verificações
    function chi_squared(params, yd, r, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, profile) result(chi2)
        real(kind=dp), intent(in) :: params(2)  ! params(1) = rs, params(2) = rho_s
        real(kind=dp), intent(in) :: yd
        real(kind=dp), intent(in) :: r(:), Vobs(:), Vgas(:), Vdisk_orig(:), e_Vobs(:)
        real(kind=dp), intent(in) :: mu_yd, sigma_yd
        character(len=*), intent(in) :: profile
        real(kind=dp) :: chi2
        
        real(kind=dp) :: Vdisk(size(r)), Vh(size(r)), Vtotal(size(r))
        real(kind=dp) :: chi2_data, chi2_prior
        integer :: i, n
        
        ! Verificação inicial de parâmetros físicos
        if (params(1) <= 0.0_dp .or. params(2) <= 0.0_dp .or. yd <= 0.0_dp) then
            chi2 = huge(1.0_dp)
            return
        end if
        
        n = size(r)
        
        ! Verifica consistência dos tamanhos dos arrays
        if (size(Vobs) /= n .or. size(Vgas) /= n .or. size(Vdisk_orig) /= n .or. size(e_Vobs) /= n) then
            chi2 = huge(1.0_dp)
            return
        end if
        
        ! Calcula componentes da curva de rotação
        Vdisk = yd * Vdisk_orig
        
        select case(trim(profile))
        case("Burk")
            Vh = burkert_halo(params, r)
        case("NFW")
            Vh = nfw_halo(params, r)
        case default
            chi2 = huge(1.0_dp)
            return
        end select
        
        ! Verifica valores inválidos usando a função portável
        if (has_nan(Vh)) then
            chi2 = huge(1.0_dp)
            return
        end if
        
        Vtotal = sqrt(Vdisk**2 + Vh**2 + Vgas**2)
        
        ! Calcula chi-quadrado dos dados
        chi2_data = 0.0_dp
        do i = 1, n
            if (e_Vobs(i) > tiny_val .and. Vtotal(i) > 0.0_dp) then
                chi2_data = chi2_data + ((Vobs(i) - Vtotal(i))**2) / (e_Vobs(i)**2)
            else
                chi2_data = huge(1.0_dp)
                exit
            end if
        end do
        
        ! Adiciona termo de prior para yd
        if (sigma_yd > tiny_val) then
            chi2_prior = ((yd - mu_yd)**2) / (sigma_yd**2)
        else
            if (abs(yd - mu_yd) > tiny_val) then
                chi2_prior = huge(1.0_dp)
            else
                chi2_prior = 0.0_dp
            end if
        end if
        
        chi2 = chi2_data + chi2_prior
    end function chi_squared

    ! Função de densidade de probabilidade não normalizada
    function unnorm_pdf(params, r, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, profile) result(pdf_val)
        real(kind=dp), intent(in) :: params(3)  ! params(1) = YD, params(2) = rs, params(3) = rho_s
        real(kind=dp), intent(in) :: r(:), Vobs(:), Vgas(:), Vdisk_orig(:), e_Vobs(:)
        real(kind=dp), intent(in) :: mu_yd, sigma_yd
        character(len=*), intent(in) :: profile
        real(kind=dp) :: pdf_val
        
        real(kind=dp) :: chi2_val
        real(kind=dp) :: halo_params(2)
        
        halo_params(1) = params(2)  ! rs
        halo_params(2) = params(3)  ! rho_s
        
        chi2_val = chi_squared(halo_params, params(1), r, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, profile)
        
        ! Usa log para evitar underflow, depois converte de volta
        pdf_val = exp(-0.5_dp * chi2_val)
    end function unnorm_pdf

    ! Calcula incertezas (percentis 16 e 84)
    function calculate_uncertainty(values, pdf) result(uncertainty)
        real(kind=dp), intent(in) :: values(:), pdf(:)
        real(kind=dp) :: uncertainty
        real(kind=dp), allocatable :: pdf_cumsum(:)
        real(kind=dp) :: total_pdf, lower_bound, upper_bound
        integer :: n
        
        n = size(values)
        if (size(pdf) /= n) then
            error stop "Erro em calculate_uncertainty: tamanhos de arrays diferentes"
        end if
        
        if (n < 2) then
            uncertainty = 0.0_dp
            return
        end if
        
        allocate(pdf_cumsum(n))
        
        ! Soma cumulativa e normalização
        call cumsum(pdf, pdf_cumsum)
        total_pdf = pdf_cumsum(n)
        
        if (total_pdf <= tiny_val) then
            uncertainty = 0.0_dp
            deallocate(pdf_cumsum)
            return
        end if
        
        pdf_cumsum = pdf_cumsum / total_pdf
        
        ! Encontra os percentis 16 e 84
        lower_bound = interp(0.16_dp, pdf_cumsum, values)
        upper_bound = interp(0.84_dp, pdf_cumsum, values)
        
        uncertainty = (upper_bound - lower_bound) / 2.0_dp
        deallocate(pdf_cumsum)
    end function calculate_uncertainty
end module stats_mod
