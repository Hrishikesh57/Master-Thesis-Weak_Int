program print_polynomial_table
    implicit none

    ! Declare variables
    integer :: x
    real :: resultG , resultF


    ! Print the header of the table
    print *, "  x     G(x)    F(x)"
    print *, "--------------
    ! Loop over x from 1 to 100
    do x = 1, 100
        resultG = G(real(x))  ! Calculate the value of G(x)
        resultF = F(real(x))    ! Calculate the value of F(x)
        print '(I3, 5X, F7.2,5X,F10.2)', x, resultG , resultF ! Print x and G(x)
    end do

contains

    ! Define the polynomial function G(x)
    real function G(x)
        real, intent(in) :: x
        G = x**2 + 2.0*x + 1.0
    end function G
    ! Define the polynomial function F(x)
    real function F(x)
        real, intent(in) :: x
        F = 8.0*x**3 + 7.0*x**2 + 2.0*x + 1.0
    end function F
end program print_polynomial_table



