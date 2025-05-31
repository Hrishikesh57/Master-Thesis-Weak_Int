module bspqed_momentum
    implicit none
    private

    public :: momentumgf

    double precision, parameter :: pihalf = 2.0d0 * atan(1.0d0)
    double precision, parameter :: sqrtpihalf = sqrt(pihalf)

contains

    subroutine momentumgf(knots, order, kappa, d1, d2, p, gp, fp)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(0:, :), intent(in) :: d1
        double precision, dimension(0:, :), intent(in) :: d2
        double precision, intent(in) :: p

        double precision, intent(out) :: gp
        double precision, intent(out) :: fp

        double precision, dimension(0:order, 0:max((abs(2*kappa+1)-1)/2, (abs(2*kappa-1)-1)/2)) :: j0, j1
        double precision :: sum1, sum2, a, b, gmi, fmi, pm
        integer :: nknots, npoints, nintervals, ileft, l1, l2, i, m

        if (p == 0.0d0) then
            gp = 0.0d0
            fp = 0.0d0
            return
        end if

        nknots = size(knots)
        npoints = nknots - 2 * (order - 1)
        nintervals = npoints - 1

        l1 = (abs(2 * kappa + 1) - 1) / 2
        l2 = (abs(2 * kappa - 1) - 1) / 2

        sum1 = 0.0d0
        sum2 = 0.0d0
        j0 = 0.0d0
        do i = 1, nintervals
            ileft = order - 1 + i
            a = knots(ileft)
            b = knots(ileft+1)
            call sphjrec(p*b, j1)
            pm = 1.0d0
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
                sum1 = sum1 + gmi * (j1(m, l1) - j0(m, l1)) / pm
                sum2 = sum2 + fmi * (j1(m, l2) - j0(m, l2)) / pm
                pm = pm * p
            end do
            j0 = j1
        end do
        gp = sum1 / sqrtpihalf
        fp = -kappa / abs(kappa) * sum2 / sqrtpihalf
        return
    end subroutine momentumgf

    subroutine sphjrec(x, j)
        double precision, intent(in) :: x

        double precision, dimension(0:, 0:), intent(out) :: j

        double precision, dimension(-1:(size(j,1)-1)) :: s
        double precision, dimension( 0:(size(j,1)-1)) :: c
        double precision, dimension( 0:(size(j,2)-1)) :: sphj
        double precision :: sphjxk, l2p1, l2p3, l2p5, res
        integer :: kmax, lmax, k, l

        kmax = size(j, 1) - 1
        lmax = size(j, 2) - 1

        if (x < 0.5d0) then
            l2p1 = 1.0d0
            do l = 0, lmax
                l2p1 = l2p1 * (2 * l + 1)
                l2p3 = l2p1 * (2 * l + 3)
                l2p5 = l2p3 * (2 * l + 5)
                do k = 0, kmax
                    res = x**(k+l+1) / (l2p1 * (k + l + 1))
                    res = res - 0.5d0 * x**(k+l+3) / (l2p3 * (k + l + 3))
                    res = res + 0.125d0 * x**(k+l+5) / (l2p5 * (k + l + 5))
                    j(k, l) = res
                end do
            end do
            return
        end if

        s(-1) = sinint(x)
        call screc(x, s(0:), c)
        do k = 0, kmax
            j(k, 0) = s(k-1)
        end do
        call sphbesj(x, sphj)
        do l = 1, lmax
            j(0, l) = (j(1, l-1) - x * sphj(l)) / l
            sphjxk = x * sphj(l-1)
            do k = 1, kmax
                j(k, l) = (l + k - 1) * j(k-1, l-1) - sphjxk
                sphjxk = sphjxk * x
            end do
        end do
        return
    end subroutine sphjrec

    subroutine screc(x, s, c)
        double precision, intent(in) :: x

        double precision, dimension(0:), intent(out) :: s
        double precision, dimension(0:), intent(out) :: c

        double precision :: x1, x2, x3, x4, x5, x6, sinxk, cosxk
        integer :: n, k

        n = min(size(s), size(c)) - 1

        if (x < 0.05d0) then
            x1 = x
            x2 = x1 * x
            x3 = x2 * x
            x4 = x3 * x
            x5 = x4 * x
            x6 = x5 * x
            do k = 0, n
                s(k) = x2 / (k + 2) - x4 / (6 * (k + 4)) + x6 / (120 * (k + 6))
                c(k) = x1 / (k + 1) - x3 / (2 * (k + 3)) + x5 / (24 * (k + 5))
                x1 = x1 * x
                x2 = x2 * x
                x3 = x3 * x
                x4 = x4 * x
                x5 = x5 * x
                x6 = x6 * x
            end do
        else
            sinxk = sin(x)
            cosxk = cos(x)
            s(0) = 2.0d0 * sin(0.5d0 * x)**2
            c(0) = sinxk
            do k = 1, n
                sinxk = sinxk * x
                cosxk = cosxk * x
                s(k) = k * c(k-1) - cosxk
                c(k) = sinxk - k * s(k-1)
            end do
        end if
        return
    end subroutine screc

    subroutine sphbesj(x, sphj)
        double precision, intent(in) :: x

        double precision, dimension(0:), intent(out) :: sphj

        integer :: n, nz

        if (x == 0.0d0) then
            sphj(0) = 1.0d0
            sphj(1:) = 0.0d0
            return
        end if
        n = size(sphj)
        call dbesj(x, 0.5d0, n, sphj, nz)
        sphj = sqrt(pihalf/x) * sphj
        return
    end subroutine sphbesj

    double precision function sinint(x)
        double precision, intent(in) :: x

        double precision, dimension(8), parameter :: csup = [&
            &  1.00000000000000000d+00, &
            & -4.54393409816329991d-02, &
            &  1.15457225751016682d-03, &
            & -1.41018536821330254d-05, &
            &  9.43280809438713025d-08, &
            & -3.53201978997168357d-10, &
            &  7.08240282274875911d-13, &
            & -6.05338212010422477d-16  &
        &]

        double precision, dimension(7), parameter :: csdw = [&
            &  1.00000000000000000d+00, &
            &  1.01162145739225565d-02, &
            &  4.99175116169755106d-05, &
            &  1.55654986308745614d-07, &
            &  3.28067571055789734d-10, &
            &  4.50490975753865810d-13, &
            &  3.21107051193712168d-16  &
        &]

        double precision, dimension(11), parameter :: cfup = [&
            &  1.00000000000000000000d+00, &
            &  7.44437068161936700618d+02, &
            &  1.96396372895146869801d+05, &
            &  2.37750310125431834034d+07, &
            &  1.43073403821274636888d+09, &
            &  4.33736238870432522765d+10, &
            &  6.40533830574022022911d+11, &
            &  4.20968180571076940208d+12, &
            &  1.00795182980368574617d+13, &
            &  4.94816688199951963482d+12, &
            & -4.94701168645415959931d+11  &
        &]

        double precision, dimension(10), parameter :: cfdw = [&
            &  1.00000000000000000000d+00, &
            &  7.46437068161927678031d+02, &
            &  1.97865247031583951450d+05, &
            &  2.41535670165126845144d+07, &
            &  1.47478952192985464958d+09, &
            &  4.58595115847765779830d+10, &
            &  7.08501308149515401563d+11, &
            &  5.06084464593475076774d+12, &
            &  1.43468549171581016479d+13, &
            &  1.11535493509914254097d+13  &
        &]

        double precision, dimension(11), parameter :: cgup = [&
            &  1.00000000000000000d+00, &
            &  8.13595201151686150d+02, &
            &  2.35239181626478200d+05, &
            &  3.12557570795778731d+07, &
            &  2.06297595146763354d+09, &
            &  6.83052205423625007d+10, &
            &  1.09049528450362786d+12, &
            &  7.57664583257834349d+12, &
            &  1.81004487464664575d+13, &
            &  6.43291613143049485d+12, &
            & -1.36517137670871689d+12  &
        &]

        double precision, dimension(10), parameter :: cgdw = [&
            &  1.00000000000000000d+00, &
            &  8.19595201151451564d+02, &
            &  2.40036752835578777d+05, &
            &  3.26026661647090822d+07, &
            &  2.23355543278099360d+09, &
            &  7.87465017341829930d+10, &
            &  1.39866710696414565d+12, &
            &  1.17164723371736605d+13, &
            &  4.01839087307656620d+13, &
            &  3.99653257887490811d+13  &
        &]

        double precision :: sgn, xi, arg, pup, pdw, fvl, gvl

        sgn = 1.0d0
        if (x < 0.0d0) sgn = -1.0d0
        xi = sgn * x
        if (xi <= 4.0d0) then
            arg = xi**2
            pup = polyeval(csup, arg)
            pdw = polyeval(csdw, arg)
            sinint = xi * pup / pdw
        else
            arg = 1.0d0 / xi**2
            pup = polyeval(cfup, arg)
            pdw = polyeval(cfdw, arg)
            fvl = pup / pdw / xi
            pup = polyeval(cgup, arg)
            pdw = polyeval(cgdw, arg)
            gvl = pup / pdw / xi**2
            sinint = pihalf - fvl * cos(xi) - gvl * sin(xi)
        end if
        sinint = sgn * sinint
        return
    end function sinint

    double precision function polyeval(coef, x)
        double precision, dimension(:), intent(in) :: coef
        double precision, intent(in) :: x
        double precision :: y
        integer :: n, i
        n = size(coef)
        y = coef(n)
        do i = n-1, 1, -1
            y = y * x + coef(i)
        end do
        polyeval = y
        return
    end function polyeval

end module bspqed_momentum
