module bspqed_coefficients
    implicit none
    private

    public :: bspppf
    public :: ppfpbn
    public :: bsppbn
    public :: bspgfc

contains

    subroutine bspppf(knots, order, coefs, c1, c2)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:) :: coefs

        double precision, dimension(0:, :), intent(out) :: c1
        double precision, dimension(0:, :), intent(out) :: c2

        double precision, dimension(0:(order-1), size(knots)-2*(order-1)-1) &
        & :: c
        double precision, dimension(size(coefs)/2) :: bcoefs
        double precision, dimension(size(knots)-2*(order-1)) :: xi
        double precision, dimension(order*(size(knots)-order+3)) :: work
        double precision, dimension(0:(order-1)) :: factorial
        integer :: nknots, ordm1, npoints, nintervals, nbspl, nxi, j

        nknots = size(knots)
        ordm1 = order - 1
        npoints = nknots - 2 * ordm1
        nintervals = npoints - 1
        nbspl = nknots - order
        factorial(0) = 1.0d0
        do j = 1, ordm1
            factorial(j) = j * factorial(j-1)
        end do
        bcoefs = coefs(:nbspl)
        call dbsppp(knots, bcoefs, nbspl, order, order, c, xi, nxi, work)
        do j = 1, nintervals
            c1(:, j) = c(:, j) / factorial
        end do
        bcoefs = coefs((nbspl+1):)
        call dbsppp(knots, bcoefs, nbspl, order, order, c, xi, nxi, work)
        do j = 1, nintervals
            c2(:, j) = c(:, j) / factorial
        end do
        return
    end subroutine

    subroutine ppfpbn(knots, order, c1, c2, d1, d2)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(0:, :), intent(in) :: c1
        double precision, dimension(0:, :), intent(in) :: c2

        double precision, dimension(0:, :), intent(out) :: d1
        double precision, dimension(0:, :), intent(out) :: d2

        double precision, external :: dbinom
        double precision :: minusa, sum1, sum2, bnmc, mapw
        integer :: ordm1, nintervals, jleft, j, m, k

        ordm1 = order - 1
        nintervals = size(c1, 2)
        do j = 1, nintervals
            jleft = order - 1 + j
            minusa = -knots(jleft)
            do m = 0, ordm1
                sum1 = 0.0d0
                sum2 = 0.0d0
                do k = m, ordm1
                    bnmc = dbinom(k, m)
                    mapw = minusa**(k-m)
                    sum1 = sum1 + bnmc * mapw * c1(k, j)
                    sum2 = sum2 + bnmc * mapw * c2(k, j)
                end do
                d1(m, j) = sum1
                d2(m, j) = sum2
            end do
        end do
        return
    end subroutine ppfpbn

    subroutine bsppbn(knots, order, coefs, d1, d2)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:) :: coefs

        double precision, dimension(0:, :), intent(out) :: d1
        double precision, dimension(0:, :), intent(out) :: d2

        double precision, dimension(order, size(knots)-2*(order-1)-1) :: c1, c2

        call bspppf(knots, order, coefs, c1, c2)
        call ppfpbn(knots, order, c1, c2, d1, d2)
        return
    end subroutine bsppbn

    subroutine bspgfc(knots, order, kappa, coefs, gcf, fcf)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:) :: coefs

        double precision, dimension(0:, :), intent(out) :: gcf
        double precision, dimension(0:, :), intent(out) :: fcf

        double precision, dimension(0:(order-1), size(knots)-2*(order-1)-1) :: &
        & d1, d2
        double precision :: a, b, gmi, fmi
        integer :: nintervals, i, ileft, m

        call bsppbn(knots, order, coefs, d1, d2)
        nintervals = size(knots) - 2 * (order - 1) - 1
        do i = 1, nintervals
            ileft = order - 1 + i
            a = knots(ileft)
            b = knots(ileft+1)
            do m = 0, order
                if (m == 0) then
                    gmi = -0.5d0 * kappa * d2(0, i)
                    fmi =  0.5d0 * kappa * d1(0, i)
                else if (m == order) then
                    gmi = d1(order-1, i)
                    fmi = d2(order-1, i)
                else
                    gmi = d1(m-1, i) + 0.5d0 * (m - kappa) * d2(m, i)
                    fmi = d2(m-1, i) + 0.5d0 * (m + kappa) * d1(m, i)
                end if
                gcf(m, i) = gmi
                fcf(m, i) = fmi
            end do
        end do
        return
    end subroutine bspgfc

end module bspqed_coefficients
