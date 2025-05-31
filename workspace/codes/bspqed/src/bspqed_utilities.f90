module bspqed_utilities
    use bspqed_function, only: function_type
    use bspqed_integrals, only: intlw
    implicit none
    private

    public :: gam
    public :: energy
    public :: meanr
    public :: ptsknts
    public :: brkpts1
    public :: kntseq1
    public :: perturbation

contains

    double precision function gam(zalpha)
        double precision, intent(in) :: zalpha
        gam = sqrt(1.0d0 - zalpha**2)
        return
    end function gam

    double precision function energy(zalpha, kappa, nmain)
        double precision, intent(in) :: zalpha
        integer, intent(in) :: kappa
        integer, intent(in) :: nmain

        double precision :: gamfct, radicn

        gamfct = sqrt(kappa**2 - zalpha**2)
        radicn = 1.0d0 + (zalpha / (nmain - abs(kappa) + gamfct))**2
        energy = 1.0d0 / sqrt(radicn)
        return
    end function energy

    double precision function meanr(zalpha, kappa, nmain)
        double precision, intent(in) :: zalpha
        integer, intent(in) :: kappa
        integer, intent(in) :: nmain

        double precision :: eps, sm1, sm2

        eps = energy(zalpha, kappa, nmain)
        sm1 = 1.5d0 * zalpha * eps / (1.0d0 - eps**2)
        sm2 = -0.5d0 * kappa * (1.0d0 + kappa * eps) / zalpha
        meanr = sm1 + sm2
        return
    end function meanr

    subroutine ptsknts(points, order, knots)
        double precision, dimension(:), intent(in) :: points
        integer, intent(in) :: order

        double precision, dimension(:), intent(out) :: knots

        integer :: npoints, nknots

        npoints = size(points)
        nknots = npoints + 2 * (order - 1)
        knots(1:order) = points(1)
        knots((order+1):(nknots-order)) = points(2:(npoints-1))
        knots((nknots-order+1):nknots) = points(npoints)
        return
    end subroutine ptsknts

    subroutine brkpts1(rnucl, rcavity, ninside, points)
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity
        integer, intent(in) :: ninside

        double precision, dimension(:), intent(out) :: points

        double precision :: ratio
        integer :: npoints, noutside, i

        npoints = size(points)
        noutside = npoints - ninside
        ratio = rcavity / rnucl
        points(1) = 0.0d0
        do i = 2, ninside
            points(i) = rnucl * (i-1) / dble(ninside-1)
        end do
        do i = 1, noutside
            points(ninside + i) = rnucl * ratio**(i / dble(noutside))
        end do
        return
    end subroutine brkpts1

    subroutine kntseq1(rnucl, rcavity, ninside, order, knots)
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity
        integer, intent(in) :: ninside
        integer, intent(in) :: order

        double precision, dimension(:), intent(out) :: knots

        double precision, dimension(size(knots)-2*(order-1)) :: points

        call brkpts1(rnucl, rcavity, ninside, points)
        call ptsknts(points, order, knots)
        return
    end subroutine kntseq1

    subroutine perturbation(knots, order, na, kappa, spectrum, coefmat, &
    & potential, nodes, weights, pertcoef)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: na
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: spectrum
        double precision, dimension(:, :), intent(in) :: coefmat
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: pertcoef

        double precision, dimension(size(coefmat, 1)) :: coefa
        double precision :: epsa, epsn, diff, vna
        integer :: nspectrum, nn

        pertcoef = 0.0d0
        nspectrum = size(spectrum)
        coefa = coefmat(:, na)
        epsa = spectrum(na)
        do nn = 1, nspectrum
            epsn = spectrum(nn)
            diff = epsa - epsn
            if (abs(diff) == 0.0d0 .or. nn == na) cycle
            vna = intlw(knots, order, kappa, coefmat(:, nn), kappa, coefa, &
            & potential, nodes, weights)
            pertcoef = pertcoef + vna * coefmat(:, nn) / diff
        end do
        return
    end subroutine perturbation

end module bspqed_utilities
