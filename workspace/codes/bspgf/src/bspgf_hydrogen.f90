module bspgf_hydrogen
    use bspqed_function, only: function_type
    use bspqed_integrals, only: intlv
    implicit none
    private

    public :: gfactor

contains

    double precision function gfactor(knots, order, kappa, coefa, coefb, magnpot, nodes, weights)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefa
        double precision, dimension(:), intent(in) :: coefb
        class(function_type), intent(in) :: magnpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: angu, lvab
        integer :: twoj

        twoj = 2 * abs(kappa) - 1
        angu = 4.0d0 * kappa / (twoj * (twoj + 2.0d0))
        lvab = intlv(knots, order, kappa, coefa, kappa, coefb, magnpot, nodes, weights)
        gfactor = -angu * lvab
        return
    end function gfactor

end module bspgf_hydrogen
