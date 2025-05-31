module bspse_zeropotential
    use bspqed_momentum, only: momentumgf
    implicit none
    private

    public :: zeropot

contains

    subroutine zeropot(knots, order, kappa, epsa, d1, d2, res, abserr, epsabs, epsrel)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, intent(in) :: epsa
        double precision, dimension(:, :), intent(in) :: d1
        double precision, dimension(:, :), intent(in) :: d2

        double precision, intent(out) :: res
        double precision, intent(out) :: abserr

        double precision, optional, intent(in) :: epsabs
        double precision, optional, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        integer, dimension(limit) :: iwork
        double precision, dimension(lenw) :: work

        double precision :: epsab, epsrl
        integer :: neval, ier, last

        epsab = 1.0d-10
        if (present(epsabs)) epsab = epsabs
        epsrl = 5.0d-07
        if (present(epsrel)) epsrl = epsrel
        call dqagi(f, 0.0d0, 1, epsab, epsrl, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        return

    contains

        double precision function f(p)
            double precision, intent(in) :: p

            double precision :: rho, lnrho, arho, brho, gp, fp

            rho = 1.0d0 - epsa**2 + p**2
            lnrho = log(rho)
            arho = 2.0d0 * (1.0d0 + 2.0d0 * rho * lnrho / (1.0d0 - rho))
            brho = -(2.0d0 - rho) * (1.0d0 + rho * lnrho / (1.0d0 - rho)) / (1.0d0 - rho)
            call momentumgf(knots, order, kappa, d1, d2, p, gp, fp)
            f = 0.25d0 * (arho * (gp**2 - fp**2) + brho * (epsa * (gp**2 + fp**2) + 2.0d0 * p * gp * fp))
            return
        end function f
    end subroutine zeropot

end module bspse_zeropotential
