program galaxy_fitter
    use mpi_f08
    use constants_mod
    use galaxy_processor_mod
    implicit none
    
    integer :: rank, num_procs, ierr
    integer :: i, file_idx
    integer :: num_files, files_per_proc, remainder
    integer, allocatable :: counts(:), displs(:)
    
    real(kind=dp), parameter :: mu_yd = 0.5_dp
    real(kind=dp), parameter :: sigma_yd = 0.1_dp
    character(len=256) :: data_dir = "./data"
    character(len=256) :: output_dir = "./results"
    
    ! Lista de arquivos de galáxias
    character(len=12), dimension(4) :: galaxy_files = [ &
        "NGC0055.txt  ", &
        "NGC2366.txt  ", &
        "NGC3109.txt  ", &
        "UGC05721.txt "  &
    ]
    character(len=256) :: input_file, output_prefix
    character(len=256) :: filename
    
    real(kind=dp) :: start_time, end_time
    
    ! Inicializa MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    
    num_files = size(galaxy_files)
    
    if (rank == 0) then
        print *, "Iniciando ajuste de galáxias com MPI"
        print *, "Número de processos: ", num_procs
        print *, "Arquivos para processar: "
        do i = 1, num_files
            print *, trim(galaxy_files(i))
        end do
        start_time = MPI_Wtime()
    end if
    
    ! Distribuição balanceada de trabalho usando MPI_SCATTERV
    allocate(counts(0:num_procs-1), displs(0:num_procs-1))
    
    files_per_proc = num_files / num_procs
    remainder = mod(num_files, num_procs)
    
    do i = 0, num_procs-1
        if (i < remainder) then
            counts(i) = files_per_proc + 1
        else
            counts(i) = files_per_proc
        end if
    end do
    
    displs(0) = 0
    do i = 1, num_procs-1
        displs(i) = displs(i-1) + counts(i-1)
    end do
    
    ! Cada processo processa seus arquivos atribuídos
    do file_idx = displs(rank)+1, displs(rank)+counts(rank)
        if (file_idx > num_files) cycle
        
        ! Constrói caminhos de entrada e saída
        input_file = trim(data_dir) // "/" // trim(adjustl(galaxy_files(file_idx)))
        
        ! Remove extensão .txt para o prefixo de saída
        filename = trim(adjustl(galaxy_files(file_idx)))
        i = index(filename, ".txt")
        if (i > 0) filename = filename(:i-1)
        output_prefix = trim(output_dir) // "/" // trim(filename)
        
        ! Processa esta galáxia
        call process_galaxy(input_file, output_prefix, mu_yd, sigma_yd, rank)
    end do
    
    ! Barreira para sincronização
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    ! Finalização
    if (rank == 0) then
        end_time = MPI_Wtime()
        print *, "Processamento concluído"
        print *, "Galáxias processadas: ", num_files
        print *, "Tempo total: ", end_time - start_time, " segundos"
    end if
    
    deallocate(counts, displs)
    call MPI_Finalize(ierr)
end program galaxy_fitter
