module bspqed_fourier_density
    use bspqed_function, only: function_type
    implicit none
    private

    public :: coulomb_density_fourier
    public :: sphere_density_fourier

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: sqrt2picube = sqrt(2.0d0 * pi)**3
    double precision, parameter :: pihalf = 0.5d0 * pi

    type, extends(function_type) :: coulomb_density_fourier
        double precision :: zalpha
    contains
        procedure, pass(self) :: eval => coulomb_density_fourier_eval
    end type coulomb_density_fourier

    type, extends(function_type) :: sphere_density_fourier
        double precision :: zalpha
        double precision :: rad
    contains
        procedure, pass(self) :: eval => sphere_density_fourier_eval
    end type sphere_density_fourier

contains

    double precision function coulomb_density_fourier_eval(x, self)
        double precision, intent(in) :: x
        class(coulomb_density_fourier), intent(in) :: self

        double precision :: zalpha

        zalpha = self%zalpha
        coulomb_density_fourier_eval = -zalpha / sqrt2picube
        return
    end function coulomb_density_fourier_eval

    double precision function sphere_density_fourier_eval(x, self)
        double precision, intent(in) :: x
        class(sphere_density_fourier), intent(in) :: self

        double precision :: zalpha, rad, qrad

        zalpha = self%zalpha
        rad = self%rad
        qrad = x * rad
        sphere_density_fourier_eval = -zalpha / sqrt2picube * 3.0d0 * sj1dx(qrad)
        return
    end function sphere_density_fourier_eval

    double precision function sj1dx(x)
        double precision, intent(in) :: x

        double precision, dimension(1) :: j1
        double precision :: y
        integer :: nz

        if (x < 0.1d0) then
            y = x**2
            sj1dx = (1.0d0 - y * (1.0d0 - y * (1.0d0 - y / 54.0d0) / 28.0d0) / 10.0d0) / 3.0d0
            return
        end if
        call dbesj(x, 1.5d0, 1, j1, nz)
        j1 = sqrt(pihalf/x) * j1
        sj1dx = j1(1) / x
        return
    end function sj1dx

end module bspqed_fourier_density
