module math_utils_mod
    use constants_mod
    implicit none
contains
    ! Soma cumulativa de um array
    subroutine cumsum(arr, result)
        real(kind=dp), intent(in) :: arr(:)
        real(kind=dp), intent(out) :: result(:)
        integer :: i, n
        
        n = size(arr)
        if (size(result) /= n) then
            error stop "Erro em cumsum: tamanhos de array incompatíveis"
        end if
        
        if (n > 0) then
            result(1) = arr(1)
            do i = 2, n
                result(i) = result(i-1) + arr(i)
            end do
        end if
    end subroutine cumsum

    ! Interpolação linear com verificação de ordenação
    function interp(x_target, x_arr, y_arr) result(y_target)
        real(kind=dp), intent(in) :: x_target
        real(kind=dp), intent(in) :: x_arr(:), y_arr(:)
        real(kind=dp) :: y_target
        integer :: i, n
        logical :: is_ascending
        
        n = size(x_arr)
        
        ! Verificações de entrada
        if (size(y_arr) /= n) then
            error stop "Erro em interp: tamanhos de x_arr e y_arr diferentes"
        end if
        
        if (n == 0) then
            error stop "Erro em interp: arrays vazios"
        end if
        
        ! Caso especial: array com um único elemento
        if (n == 1) then
            y_target = y_arr(1)
            return
        end if
        
        ! Verifica se o array está ordenado
        is_ascending = x_arr(n) > x_arr(1)
        
        ! Verifica limites
        if ((is_ascending .and. x_target <= x_arr(1)) .or. &
            (.not. is_ascending .and. x_target >= x_arr(1))) then
            y_target = y_arr(1)
            return
        end if
        
        if ((is_ascending .and. x_target >= x_arr(n)) .or. &
            (.not. is_ascending .and. x_target <= x_arr(n))) then
            y_target = y_arr(n)
            return
        end if
        
        ! Encontra o intervalo e interpola
        if (is_ascending) then
            do i = 1, n-1
                if (x_target >= x_arr(i) .and. x_target <= x_arr(i+1)) then
                    y_target = y_arr(i) + (y_arr(i+1) - y_arr(i)) * &
                              (x_target - x_arr(i)) / (x_arr(i+1) - x_arr(i))
                    return
                end if
            end do
        else
            do i = 1, n-1
                if (x_target <= x_arr(i) .and. x_target >= x_arr(i+1)) then
                    y_target = y_arr(i) + (y_arr(i+1) - y_arr(i)) * &
                              (x_target - x_arr(i)) / (x_arr(i+1) - x_arr(i))
                    return
                end if
            end do
        end if
        
        ! Se chegou aqui, algo deu errado
        error stop "Erro em interp: não foi possível encontrar o intervalo (array não ordenado?)"
    end function interp

    ! Integração trapezoidal com verificação de monotonicidade
    function trapz(y, x) result(integral_val)
        real(kind=dp), intent(in) :: y(:), x(:)
        real(kind=dp) :: integral_val
        integer :: i, n
        
        n = size(y)
        if (size(x) /= n) then
            error stop "Erro em trapz: tamanhos de x e y diferentes"
        end if
        
        integral_val = 0.0_dp
        if (n > 1) then
            do i = 1, n-1
                ! Verifica se os pontos estão em ordem crescente ou decrescente
                if ((x(i+1) > x(i)) .or. (x(i+1) < x(i))) then
                    integral_val = integral_val + 0.5_dp * (y(i) + y(i+1)) * (x(i+1) - x(i))
                else
                    error stop "Erro em trapz: valores de x não são monotônicos"
                end if
            end do
        end if
    end function trapz

    ! Gera array logspace (valores igualmente espaçados em escala log)
    subroutine logspace(start, stop, num, result)
        real(kind=dp), intent(in) :: start, stop
        integer, intent(in) :: num
        real(kind=dp), intent(out) :: result(:)
        real(kind=dp) :: step, exponent
        integer :: i
        
        if (size(result) /= num) then
            error stop "Erro em logspace: tamanho do array de saída incorreto"
        end if
        
        if (num == 1) then
            result(1) = 10.0_dp**start
            return
        end if
        
        step = (stop - start) / real(num - 1, dp)
        do i = 1, num
            exponent = start + real(i - 1, dp) * step
            result(i) = 10.0_dp**exponent
        end do
    end subroutine logspace

    ! Gera array com valores igualmente espaçados
    subroutine linspace(start, stop, num, result)
        real(kind=dp), intent(in) :: start, stop
        integer, intent(in) :: num
        real(kind=dp), intent(out) :: result(:)
        real(kind=dp) :: step
        integer :: i
        
        if (size(result) /= num) error stop "Erro em linspace: tamanho do array incorreto"
        
        if (num == 1) then
            result(1) = start
            return
        end if
        
        step = (stop - start) / real(num - 1, dp)
        do i = 1, num
            result(i) = start + (i - 1) * step
        end do
    end subroutine linspace
end module math_utils_mod
