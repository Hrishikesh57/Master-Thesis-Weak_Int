module bspqed_discretization
    use bspqed_dkb, only: radialgf
    use bspqed_function, only: function_type
    implicit none
    private

    public :: nodwgt
    public :: storewf
    public :: storefunc
    public :: matelems

contains

    subroutine nodwgt(points, nodes, weights, xs, ws)
        double precision, dimension(:), intent(in) :: points
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:), intent(out) :: xs
        double precision, dimension(:), intent(out) :: ws

        double precision :: left, rght, hlfd
        integer :: ngl, npoints, nintervals, i, j, k

        ngl = size(nodes)
        npoints = size(points)
        nintervals = npoints - 1
        do i = 1, nintervals
            left = points(i)
            rght = points(i+1)
            hlfd = 0.5d0 * (rght - left)
            do j = 1, ngl
                k = (i - 1) * ngl + j
                xs(k) = hlfd * (nodes(j) + 1.0d0) + left
                ws(k) = hlfd * weights(j)
            end do
        end do
        return
    end subroutine nodwgt

    subroutine storewf(xs, knots, order, kappa, coefs, gs, fs)
        double precision, dimension(:), intent(in) :: xs
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs

        double precision, dimension(:), intent(out) :: gs
        double precision, dimension(:), intent(out) :: fs

        double precision :: x, gx, fx
        integer :: ngrid, inbv, i

        ngrid = size(xs)
        inbv = order
        do i = 1, ngrid
            x = xs(i)
            call radialgf(knots, order, kappa, coefs, x, gx, fx, inbv)
            gs(i) = gx
            fs(i) = fx
        end do
        return
    end subroutine storewf

    subroutine storefunc(f, xs, fs)
        class(function_type), intent(in) :: f
        double precision, dimension(:), intent(in) :: xs

        double precision, dimension(:), intent(out) :: fs

        double precision :: x, fx
        integer :: ngrid, i

        ngrid = size(xs)
        do i = 1, ngrid
            x = xs(i)
            fx = f%eval(x)
            fs(i) = fx
        end do
        return
    end subroutine storefunc

    subroutine matelems(w, f, g1, f1, g2, f2, f12)
        double precision, dimension(:), intent(in) :: w
        double precision, dimension(:), intent(in) :: f
        double precision, dimension(:, :), intent(in) :: g1
        double precision, dimension(:, :), intent(in) :: f1
        double precision, dimension(:, :), intent(in) :: g2
        double precision, dimension(:, :), intent(in) :: f2

        double precision, dimension(:, :), intent(out) :: f12

        double precision, dimension(size(g1, 1)) :: gb, fb, ga, fa
        integer :: nspeca, nspecb, na, nb

        nspeca = size(g1, 2)
        nspecb = size(g2, 2)
        do nb = 1, nspecb
            gb = g2(:, nb)
            fb = f2(:, nb)
            do na = 1, nspeca
                ga = g1(:, na)
                fa = f1(:, na)
                f12(na, nb) = sum(f * (ga * gb + fa * fb) * w)
            end do
        end do
        return
    end subroutine matelems

end module bspqed_discretization
