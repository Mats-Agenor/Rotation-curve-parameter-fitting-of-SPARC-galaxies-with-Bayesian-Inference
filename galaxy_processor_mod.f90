module galaxy_processor_mod
    use constants_mod
    use math_utils_mod
    use halo_models_mod
    use data_reader_mod
    use stats_mod
    use optimization_mod
    implicit none

contains
    ! Processa uma única galáxia
    subroutine process_galaxy(file_path, output_prefix, mu_yd, sigma_yd, rank)
        character(len=*), intent(in) :: file_path, output_prefix
        real(kind=dp), intent(in) :: mu_yd, sigma_yd
        integer, intent(in) :: rank
        
        ! Dados observacionais
        real(kind=dp), allocatable :: R(:), Vobs(:), e_Vobs(:), Vgas(:), Vdisk_orig(:)
        logical :: read_success
        
        ! Parâmetros para busca em grade
        real(kind=dp) :: best_params(3), best_chi2
        real(kind=dp) :: best_rs, best_rho_s, best_yd
        
        ! Curvas de rotação do melhor ajuste
        real(kind=dp), allocatable :: Vdisk_best(:), Vh_best(:), Vtotal_best(:)
        
        ! PDFs para cálculo de incertezas
        integer, parameter :: n_pdf = 100
        real(kind=dp) :: yd_values(n_pdf), rs_values(n_pdf), rho_s_values(n_pdf)
        real(kind=dp) :: pdf_yd(n_pdf), pdf_rs(n_pdf), pdf_rho_s(n_pdf)
        real(kind=dp) :: unc_yd, unc_rs, unc_rho_s
        
        ! Variáveis para arquivos
        integer :: unit_num, i, j
        character(len=256) :: output_file
        character(len=4), dimension(2) :: profiles = ["Burk", "NFW "]
        character(len=4) :: current_profile
        
        ! Lê dados da galáxia
        call read_galaxy_data(file_path, R, Vobs, e_Vobs, Vgas, Vdisk_orig, read_success)
        if (.not. read_success) then
            if (rank == 0) print *, "Erro ao ler dados da galáxia: ", trim(file_path)
            return
        end if
        
        ! Processa cada perfil de halo
        do i = 1, size(profiles)
            current_profile = profiles(i)
            if (rank == 0) print *, "Processando perfil: ", trim(current_profile), " para ", trim(file_path)
            
            ! Executa busca em grade (chamada corrigida)
            call grid_search(R, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, trim(current_profile), &
                            best_params, best_chi2)
            
            if (best_chi2 == huge(1.0_dp)) then
                if (rank == 0) print *, "Falha na busca em grade para perfil ", trim(current_profile)
                cycle
            end if
            
            best_rs = best_params(1)
            best_rho_s = best_params(2)
            best_yd = best_params(3)
            
            ! Calcula curvas de rotação para o melhor ajuste
            allocate(Vdisk_best(size(R)), Vh_best(size(R)), Vtotal_best(size(R)))
            Vdisk_best = best_yd * Vdisk_orig
            
            select case(trim(current_profile))
            case("Burk")
                Vh_best = burkert_halo([best_rs, best_rho_s], R)
            case("NFW")
                Vh_best = nfw_halo([best_rs, best_rho_s], R)
            end select
            
            Vtotal_best = sqrt(Vdisk_best**2 + Vh_best**2 + Vgas**2)
            
            ! Calcula PDFs para estimar incertezas
            ! PDF para YD (fixa rs e rho_s)
            call linspace(max(0.01_dp, best_yd - 0.2_dp), min(1.0_dp, best_yd + 0.2_dp), n_pdf, yd_values)
            do j = 1, n_pdf
                pdf_yd(j) = unnorm_pdf([yd_values(j), best_rs, best_rho_s], &
                                    R, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, &
                                    trim(current_profile))
            end do
            unc_yd = calculate_uncertainty(yd_values, pdf_yd)
            
            ! PDF para rs (fixa YD e rho_s)
            call linspace(max(0.1_dp, best_rs - 20.0_dp), best_rs + 20.0_dp, n_pdf, rs_values)
            do j = 1, n_pdf
                pdf_rs(j) = unnorm_pdf([best_yd, rs_values(j), best_rho_s], &
                                   R, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, &
                                   trim(current_profile))
            end do
            unc_rs = calculate_uncertainty(rs_values, pdf_rs)
            
            ! PDF para rho_s (fixa YD e rs)
            call logspace(max(1.0e7_dp, best_rho_s/10.0_dp), min(1.0e16_dp, best_rho_s*10.0_dp), n_pdf, rho_s_values)
            do j = 1, n_pdf
                pdf_rho_s(j) = unnorm_pdf([best_yd, best_rs, rho_s_values(j)], &
                                      R, Vobs, Vgas, Vdisk_orig, e_Vobs, mu_yd, sigma_yd, &
                                      trim(current_profile))
            end do
            unc_rho_s = calculate_uncertainty(rho_s_values, pdf_rho_s)
            
            ! Salva resultados em arquivo
            output_file = trim(output_prefix) // "_" // trim(current_profile) // ".txt"
            open(newunit=unit_num, file=trim(output_file), status='replace', action='write')
            
            write(unit_num, '(A)') "# Resultados do ajuste para " // trim(file_path)
            write(unit_num, '(A, A)') "# Perfil: ", trim(current_profile)
            write(unit_num, '(A, F10.4, A, F10.4)') "# rs = ", best_rs, " ± ", unc_rs
            write(unit_num, '(A, ES14.6, A, ES14.6)') "# rho_s = ", best_rho_s, " ± ", unc_rho_s
            write(unit_num, '(A, F6.4, A, F6.4)') "# yd = ", best_yd, " ± ", unc_yd
            write(unit_num, '(A, F12.4)') "# Chi2 mínimo = ", best_chi2
            write(unit_num, '(A)') "#"
            write(unit_num, '(A)') "# R[kpc]  Vobs[km/s]  e_Vobs[km/s]  Vgas[km/s]  Vdisk[km/s]  Vhalo[km/s]  Vtotal[km/s]"
            
            do j = 1, size(R)
                write(unit_num, '(7F12.4)') R(j), Vobs(j), e_Vobs(j), Vgas(j), Vdisk_best(j), Vh_best(j), Vtotal_best(j)
            end do
            
            close(unit_num)
            
            if (rank == 0) print *, "Resultados salvos em: ", trim(output_file)
            deallocate(Vdisk_best, Vh_best, Vtotal_best)
        end do
        
        ! Libera memória
        deallocate(R, Vobs, e_Vobs, Vgas, Vdisk_orig)
    end subroutine process_galaxy
end module galaxy_processor_mod
