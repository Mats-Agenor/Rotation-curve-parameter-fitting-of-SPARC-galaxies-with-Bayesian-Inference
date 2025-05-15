! Módulo para leitura de dados
module data_reader_mod
    use constants_mod
    implicit none
contains
    ! Lê dados de um arquivo de galáxia, pulando as duas primeiras linhas
    subroutine read_galaxy_data(file_path, R, Vobs, e_Vobs, Vgas, Vdisk, success)
        character(len=*), intent(in) :: file_path
        real(kind=dp), allocatable, intent(out) :: R(:), Vobs(:), e_Vobs(:), Vgas(:), Vdisk(:)
        logical, intent(out) :: success
        integer :: unit_num, io_stat, n_lines, i
        character(len=256) :: line
        
        ! Inicializa variáveis
        success = .false.
        if (allocated(R)) deallocate(R)
        if (allocated(Vobs)) deallocate(Vobs)
        if (allocated(e_Vobs)) deallocate(e_Vobs)
        if (allocated(Vgas)) deallocate(Vgas)
        if (allocated(Vdisk)) deallocate(Vdisk)
        
        ! Abre o arquivo
        open(newunit=unit_num, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "Erro ao abrir arquivo: ", trim(file_path)
            return
        end if
        
        ! Conta o número de linhas de dados (pula cabeçalho)
        n_lines = 0
        read(unit_num, '(A)', iostat=io_stat) ! Pula primeira linha (cabeçalho 1)
        read(unit_num, '(A)', iostat=io_stat) ! Pula segunda linha (cabeçalho 2)
        
        do
            read(unit_num, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (len_trim(line) > 0) n_lines = n_lines + 1
        end do
        
        ! Verifica se há dados
        if (n_lines <= 0) then
            print *, "Arquivo sem dados: ", trim(file_path)
            close(unit_num)
            return
        end if
        
        ! Aloca arrays
        allocate(R(n_lines), Vobs(n_lines), e_Vobs(n_lines), Vgas(n_lines), Vdisk(n_lines))
        
        ! Lê os dados
        rewind(unit_num)
        read(unit_num, '(A)') ! Pula cabeçalho 1
        read(unit_num, '(A)') ! Pula cabeçalho 2
        
        do i = 1, n_lines
            ! Formato: Name D R Vobs Err Vgas Vdisk Vbul SBdisk SBbul
            read(unit_num, *, iostat=io_stat) line, line, R(i), Vobs(i), e_Vobs(i), Vgas(i), Vdisk(i)
            
            if (io_stat /= 0) then
                print *, "Erro ao ler linha ", i, " do arquivo ", trim(file_path)
                deallocate(R, Vobs, e_Vobs, Vgas, Vdisk)
                close(unit_num)
                return
            end if
        end do
        
        close(unit_num)
        success = .true.
    end subroutine read_galaxy_data
end module data_reader_mod
