module optimization_mod
    use constants_mod
    use math_utils_mod
    use halo_models_mod
    use stats_mod
    implicit none
contains
    ! Busca em grade melhorada com refinamento
    subroutine grid_search(r, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, profile, &
                          best_params, best_chi2)
        real(kind=dp), intent(in) :: r(:), Vobs(:), Vgas(:), Vdisk_orig(:), e_Vobs(:)
        real(kind=dp), intent(in) :: mu_yd, sigma_yd
        character(len=*), intent(in) :: profile
        real(kind=dp), intent(out) :: best_params(3)  ! rs, rho_s, yd
        real(kind=dp), intent(out) :: best_chi2
        
        ! Parâmetros para busca em grade
        real(kind=dp) :: rs_start, rs_stop, rs_step
        real(kind=dp) :: rho_s_start, rho_s_stop, rho_s_step
        real(kind=dp) :: yd_start, yd_stop, yd_step
        
        real(kind=dp) :: current_rs, current_rho_s, current_yd
        real(kind=dp) :: current_chi2, refined_chi2
        real(kind=dp) :: halo_params(2), refined_params(3)
        integer :: i, j, k
        integer :: rs_steps, rho_s_steps, yd_steps
        
        ! Configurações da grade
        rs_start = 1.0_dp       ! kpc
        rs_stop = 300.0_dp
        rs_step = 1.0_dp
        
        rho_s_start = 1.0e6_dp   ! M_sun/kpc^3
        rho_s_stop = 1.0e15_dp
        rho_s_step = 1.0e12_dp   ! Passo menor para melhor resolução
        
        yd_start = 0.1_dp
        yd_stop = 1.0_dp
        yd_step = 0.01_dp
        
        ! Inicializa com valores ruins
        best_chi2 = huge(1.0_dp)
        best_params = -1.0_dp
        
        ! Calcula número de passos para cada parâmetro
        rs_steps = nint((rs_stop - rs_start) / rs_step)
        rho_s_steps = nint((rho_s_stop - rho_s_start) / rho_s_step)
        yd_steps = nint((yd_stop - yd_start) / yd_step)
        
        ! Loop sobre todos os parâmetros (busca em grade principal)
        do i = 0, rs_steps
            current_rs = rs_start + i * rs_step
            
            do j = 0, rho_s_steps
                current_rho_s = rho_s_start + j * rho_s_step
                
                do k = 0, yd_steps
                    current_yd = yd_start + k * yd_step
                    
                    ! Calcula chi-quadrado
                    current_chi2 = chi_squared([current_rs, current_rho_s], current_yd, r, Vobs, Vgas, Vdisk_orig, &
                                              e_Vobs, mu_yd, sigma_yd, profile)
                    
                    ! Atualiza melhor solução
                    if (current_chi2 < best_chi2) then
                        best_chi2 = current_chi2
                        best_params = [current_rs, current_rho_s, current_yd]
                    end if
                end do
            end do
        end do
        
        ! Refinamento local (busca em grade mais fina ao redor do melhor ponto)
        if (best_params(1) > 0.0_dp) then
            rs_start = max(1.0_dp, best_params(1) - 5.0_dp)
            rs_stop = best_params(1) + 5.0_dp
            rs_step = 0.1_dp
            
            rho_s_start = max(1.0e6_dp, best_params(2) - 5.0e12_dp)
            rho_s_stop = best_params(2) + 5.0e12_dp
            rho_s_step = 1.0e10_dp
            
            yd_start = max(0.1_dp, best_params(3) - 0.05_dp)
            yd_stop = min(1.0_dp, best_params(3) + 0.05_dp)
            yd_step = 0.001_dp
            
            ! Executa busca refinada
            do i = 0, nint((rs_stop - rs_start)/rs_step)
                current_rs = rs_start + i * rs_step
                
                do j = 0, nint((rho_s_stop - rho_s_start)/rho_s_step)
                    current_rho_s = rho_s_start + j * rho_s_step
                    
                    do k = 0, nint((yd_stop - yd_start)/yd_step)
                        current_yd = yd_start + k * yd_step
                        
                        current_chi2 = chi_squared([current_rs, current_rho_s], current_yd, r, Vobs, Vgas, Vdisk_orig, &
                                                  e_Vobs, mu_yd, sigma_yd, profile)
                        
                        if (current_chi2 < best_chi2) then
                            best_chi2 = current_chi2
                            best_params = [current_rs, current_rho_s, current_yd]
                        end if
                    end do
                end do
            end do
        end if
        
        ! Verifica se encontrou uma solução válida
        if (best_params(1) < 0.0_dp) then
            print *, "Aviso: busca em grade não encontrou solução válida"
        end if
    end subroutine grid_search
end module optimization_mod
