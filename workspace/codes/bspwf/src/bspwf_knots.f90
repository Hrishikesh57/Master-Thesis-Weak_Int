module bspwf_knots
    implicit none
    private

    public :: meanr_nonrel
    public :: knots_01
    public :: knots_exp

!! TODO(Halil): add knots with linear and sine distribution
!! TODO(Halil): analyze knot sequences

contains

    double precision function meanr_nonrel(orbital, nmain)
        integer, intent(in) :: orbital
        integer, intent(in) :: nmain

        meanr_nonrel = 0.5d0 * (3.0d0 * nmain**2 - orbital * (orbital + 1))
        return
    end function meanr_nonrel

    subroutine knots_01(rnucl, rcavity, ninside, order, knots)
        double precision, intent(in) :: rnucl
        double precision, intent(in) :: rcavity
        integer, intent(in) :: ninside
        integer, intent(in) :: order

        double precision, dimension(:), intent(out) :: knots

        double precision, dimension(size(knots)-2*(order-1)) :: points

        call points_01(rnucl, rcavity, ninside, points)
        call ptsknts(points, order, knots)
        return
    end subroutine knots_01

    subroutine knots_exp(scale, rcavity, order, knots)
        double precision, intent(in) :: scale
        double precision, intent(in) :: rcavity
        integer, intent(in) :: order

        double precision, dimension(:), intent(out) :: knots

        double precision, dimension(size(knots)-2*(order-1)) :: points

        call points_exp(scale, rcavity, points)
        call ptsknts(points, order, knots)
        return
    end subroutine knots_exp

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

    subroutine points_01(rnucl, rcavity, ninside, points)
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
    end subroutine points_01

    subroutine points_exp(scale, rcavity, points)
        double precision, intent(in) :: scale
        double precision, intent(in) :: rcavity

        double precision, dimension(:), intent(out) :: points

        integer :: npoints, i

        npoints = size(points)
        do i = 1, npoints
            points(i) = rcavity * (exp(scale * (i-1) / (npoints-1)) - 1.0d0) / (exp(scale) - 1.0d0)
        end do
        return
    end subroutine points_exp

end module bspwf_knots
