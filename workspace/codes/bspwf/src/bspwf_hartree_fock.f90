module bspwf_hartree_fock
    use iso_fortran_env, only: error_unit
    use bspwf_constants, only: ialpha
    implicit none
    private

    public :: shf_matrix
    public :: shf_dir_matrix
    public :: shf_exc_matrix
    public :: dhfc_matrix
    public :: dhfc_dir_matrix
    public :: dhfc_exc_matrix
    public :: dhfb_matrix
    public :: dhfcb_matrix

    double precision, parameter :: speed_of_light_in_au = ialpha

contains

    subroutine shf_matrix(knots, order, la, lb, coefb, nodes, weights, shfmat)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: la
        integer, intent(in) :: lb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: shfmat

        double precision, dimension(abs(la-lb):(la+lb)) :: c3j
        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: vdir
        double precision, dimension(abs(la-lb):(la+lb), size(knots)-order, size(nodes), size(knots)-2*(order-1)-1)&
        & :: h1, h2
        double precision, dimension(size(knots)-order-2, size(knots)-order-2) :: shfdir, shfexc
        double precision, dimension(order, 0:0) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: dla, dlb, l1min, l1max, left, right, halfd, x, w, v, b0j, b0i, fnc, u, hj, lambda
        integer :: ndim, ier, lmin, lmax, nknots, nintervals, nbsplines, ngl, l, ileft, imin, imax, m, j, k, i
        external :: drc3jj, dbspvd

        ndim = size(c3j)
        dla = dble(la)
        dlb = dble(lb)
        call drc3jj(dla, dlb, 0.0d0, 0.0d0, l1min, l1max, c3j, ndim, ier)
        if (ier /= 0) then
            write (error_unit, *) "* ERROR in 'shf_matrix':"
            write (error_unit, *) "* 3j coefficients could not be computed."
            stop
        end if
        lmin = nint(l1min)
        lmax = nint(l1max)

        call shf_dir_pot(knots, order, coefb, nodes, weights, vdir)
        call shf_exc_int(knots, order, coefb, lmin, lmax, nodes, weights, h1, h2)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)

        shfdir = 0.0d0
        shfexc = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 1, x, ileft, order, dbsx, work)
                v = vdir(m, l)
                do j = imin, imax
                    b0j = dbsx(j+order-ileft, 0)
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        fnc = v * b0i * b0j
                        shfdir(i-1, j-1) = shfdir(i-1, j-1) + fnc * w
                    end do
                end do
                u = 0.0d0
                do i = imin, imax
                    b0i = dbsx(i+order-ileft, 0)
                    u = u + coefb(i) * b0i
                end do
                do j = 2, nbsplines - 1
                    hj = 0.0d0
                    do k = lmin, lmax
                        lambda = 0.5d0 * c3j(k)**2
                        hj = hj - lambda * h1(k, j, m, l) / x**(k+1)
                        hj = hj - lambda * h2(k, j, m, l) * x**k
                    end do
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        fnc = b0i * u * hj
                        shfexc(i-1, j-1) = shfexc(i-1, j-1) + fnc * w
                    end do
                end do
            end do
        end do
        shfdir = 0.5d0 * (shfdir + transpose(shfdir))
        shfexc = 0.5d0 * (shfexc + transpose(shfexc))
        shfmat = shfdir + shfexc
        return
    end subroutine shf_matrix

    subroutine shf_dir_matrix(knots, order, coefb, nodes, weights, shfdir)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: shfdir

        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: vdir
        double precision, dimension(order, 0:0) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: left, right, halfd, x, w, v, b0j, b0i, fnc
        integer :: nknots, nintervals, nbsplines, ngl, ileft, imin, imax, l, m, j, i
        external :: dbspvd

        call shf_dir_pot(knots, order, coefb, nodes, weights, vdir)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)

        shfdir = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 1, x, ileft, order, dbsx, work)
                v = vdir(m, l)
                do j = imin, imax
                    b0j = dbsx(j+order-ileft, 0)
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        fnc = v * b0i * b0j
                        shfdir(i-1, j-1) = shfdir(i-1, j-1) + fnc * w
                    end do
                end do
            end do
        end do
        shfdir = 0.5d0 * (shfdir + transpose(shfdir))
        return
    end subroutine shf_dir_matrix

    subroutine shf_exc_matrix(knots, order, la, lb, coefb, nodes, weights, shfexc)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: la
        integer, intent(in) :: lb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: shfexc

        double precision, dimension(abs(la-lb):(la+lb)) :: c3j
        double precision, dimension(abs(la-lb):(la+lb), size(knots)-order, size(nodes), size(knots)-2*(order-1)-1)&
        & :: h1, h2
        double precision, dimension(order, 0:0) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: dla, dlb, l1min, l1max, left, right, halfd, x, w, u, b0i, hj, lambda, fnc
        integer :: ndim, ier, lmin, lmax, nknots, nintervals, nbsplines, ngl, l, ileft, imin, imax, m, j, k, i
        external :: drc3jj, dbspvd

        ndim = size(c3j)
        dla = dble(la)
        dlb = dble(lb)
        call drc3jj(dla, dlb, 0.0d0, 0.0d0, l1min, l1max, c3j, ndim, ier)
        if (ier /= 0) then
            write (error_unit, *) "* ERROR in 'shf_exc_matrix':"
            write (error_unit, *) "* 3j coefficients could not be computed."
            stop
        end if
        lmin = nint(l1min)
        lmax = nint(l1max)

        call shf_exc_int(knots, order, coefb, lmin, lmax, nodes, weights, h1, h2)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)

        shfexc = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 1, x, ileft, order, dbsx, work)
                u = 0.0d0
                do i = imin, imax
                    b0i = dbsx(i+order-ileft, 0)
                    u = u + coefb(i) * b0i
                end do
                do j = 2, nbsplines - 1
                    hj = 0.0d0
                    do k = lmin, lmax
                        lambda = 0.5d0 * c3j(k)**2
                        hj = hj - lambda * h1(k, j, m, l) / x**(k+1)
                        hj = hj - lambda * h2(k, j, m, l) * x**k
                    end do
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        fnc = b0i * u * hj
                        shfexc(i-1, j-1) = shfexc(i-1, j-1) + fnc * w
                    end do
                end do
            end do
        end do
        shfexc = 0.5d0 * (shfexc + transpose(shfexc))
        return
    end subroutine shf_exc_matrix

    subroutine dhfc_matrix(knots, order, kapa, kapb, coefb, nodes, weights, dhfmat, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: dhfmat

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1)) :: c3j
        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: vdir
        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1), size(knots)-order,&
        & size(nodes), size(knots)-2*(order-1)-1) :: h1p, h1m, h2p, h2m
        double precision, dimension(size(knots)-order) :: coefb1, coefb2
        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: dhfdir, dhfexc
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, ja, jb, l1min, l1max, left, right, halfd, x, w, v, b0j, b1j, b0i, b1i, fnc, k1, dk1&
        &, k2, dk2, gb, fb, u0, u1, v0, v1, hpj, hmj, lambda
        integer :: ndim, ier, lmin, lmax, la, lb, nknots, nintervals, nbsplines, ngl, s, l, ileft, imin, imax, m,&
        & j, i, k
        external :: drc3jj, dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        ndim = size(c3j)
        ja = 0.5d0 * (2 * abs(kapa) - 1)
        jb = 0.5d0 * (2 * abs(kapb) - 1)
        call drc3jj(ja, jb, -0.5d0, 0.5d0, l1min, l1max, c3j, ndim, ier)
        if (ier /= 0) then
            write (error_unit, *) "* ERROR in 'dhfc_matrix':"
            write (error_unit, *) "* 3j coefficients could not be computed."
            stop
        end if
        lmin = nint(l1min)
        lmax = nint(l1max)

        la = (abs(2 * kapa + 1) - 1) / 2
        lb = (abs(2 * kapb + 1) - 1) / 2

        call dhfc_dir_pot(knots, order, kapb, coefb, nodes, weights, vdir, speed_of_light=c)
        call dhfc_exc_int(knots, order, kapa, kapb, coefb, lmin, lmax, nodes, weights, h1p, h1m, h2p, h2m,&
        & speed_of_light=c)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)
        s = nbsplines - 2

        coefb1 = coefb(1:nbsplines)
        coefb2 = coefb((nbsplines+1):)

        dhfdir = 0.0d0
        dhfexc = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 2, x, ileft, order, dbsx, work)
                v = vdir(m, l)
                do j = imin, imax
                    b0j = dbsx(j+order-ileft, 0)
                    b1j = dbsx(j+order-ileft, 1)
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        b1i = dbsx(i+order-ileft, 1)

                        fnc = v * ((1.0d0 + (0.5d0 * kapa / (c * x))**2) * b0i * b0j &
                        & + 0.25d0 * (b1i * b1j + kapa / x * (b1i * b0j + b0i * b1j)) / c**2)

                        dhfdir(i-1, j-1) = dhfdir(i-1, j-1) + fnc * w

                        fnc = v * 0.5d0 * (b1i * b0j + b0i * b1j) / c
                        dhfdir(i-1, s+j-1) = dhfdir(i-1, s+j-1) + fnc * w
                        dhfdir(s+i-1, j-1) = dhfdir(s+i-1, j-1) + fnc * w

                        fnc = v * ((1.0d0 + (0.5d0 * kapa / (c * x))**2) * b0i * b0j &
                        & + 0.25d0 * (b1i * b1j - kapa / x * (b1i * b0j + b0i * b1j)) / c**2)
                        dhfdir(s+i-1, s+j-1) = dhfdir(s+i-1, s+j-1) + fnc * w
                    end do
                end do
                k1 = 0.0d0
                dk1 = 0.0d0
                k2 = 0.0d0
                dk2 = 0.0d0
                do i = imin, imax
                    b0i = dbsx(i+order-ileft, 0)
                    b1i = dbsx(i+order-ileft, 1)
                    k1 = k1 + coefb1(i) * b0i
                    dk1 = dk1 + coefb1(i) * b1i
                    k2 = k2 + coefb2(i) * b0i
                    dk2 = dk2 + coefb2(i) * b1i
                end do
                gb = k1 + 0.5d0 * (dk2 - kapb / x * k2) / c
                fb = k2 + 0.5d0 * (dk1 + kapb / x * k1) / c
                u0 = gb + 0.5d0 * kapa / (c * x) * fb
                u1 = 0.5d0 * fb / c
                v0 = fb - 0.5d0 * kapa / (c * x) * gb
                v1 = 0.5d0 * gb / c
                do j = 2, nbsplines - 1
                    hpj = 0.0d0
                    hmj = 0.0d0
                    do k = lmin, lmax
                        if (mod(la + lb + k, 2) /= 0) cycle
                        lambda = c3j(k)**2
                        hpj = hpj - lambda * h1p(k, j, m, l) / x**(k+1)
                        hpj = hpj - lambda * h2p(k, j, m, l) * x**k
                        hmj = hmj - lambda * h1m(k, j, m, l) / x**(k+1)
                        hmj = hmj - lambda * h2m(k, j, m, l) * x**k
                    end do
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        b1i = dbsx(i+order-ileft, 1)

                        fnc = (b0i * u0 + b1i * u1) * hpj
                        dhfexc(i-1, j-1) = dhfexc(i-1, j-1) + fnc * w

                        fnc = (b0i * u0 + b1i * u1) * hmj
                        dhfexc(i-1, s+j-1) = dhfexc(i-1, s+j-1) + fnc * w

                        fnc = (b0i * v0 + b1i * v1) * hpj
                        dhfexc(s+i-1, j-1) = dhfexc(s+i-1, j-1) + fnc * w

                        fnc = (b0i * v0 + b1i * v1) * hmj
                        dhfexc(s+i-1, s+j-1) = dhfexc(s+i-1, s+j-1) + fnc * w
                    end do
                end do
            end do
        end do
        dhfdir = 0.5d0 * (dhfdir + transpose(dhfdir))
        dhfexc = 0.5d0 * (dhfexc + transpose(dhfexc))
        dhfmat = dhfdir + dhfexc
        return
    end subroutine dhfc_matrix

    subroutine dhfc_dir_matrix(knots, order, kapa, kapb, coefb, nodes, weights, dhfdir, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: dhfdir

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: vdir
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, left, right, halfd, x, w, v, b0j, b1j, b0i, b1i, fnc
        integer :: nknots, nintervals, nbsplines, nspectrum, ngl, s, ileft, imin, imax, l, m, j, i
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        nspectrum = 2 * (nbsplines - 2)
        ngl = size(nodes)
        s = nbsplines - 2

        call dhfc_dir_pot(knots, order, kapb, coefb, nodes, weights, vdir, speed_of_light=c)

        dhfdir = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 2, x, ileft, order, dbsx, work)
                v = vdir(m, l)
                do j = imin, imax
                    b0j = dbsx(j+order-ileft, 0)
                    b1j = dbsx(j+order-ileft, 1)
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        b1i = dbsx(i+order-ileft, 1)

                        fnc = v * ((1.0d0 + (0.5d0 * kapa / (c * x))**2) * b0i * b0j &
                        & + 0.25d0 * (b1i * b1j + kapa / x * (b1i * b0j + b0i * b1j)) / c**2)
                        dhfdir(i-1, j-1) = dhfdir(i-1, j-1) + fnc * w

                        fnc = v * 0.5d0 * (b1i * b0j + b0i * b1j) / c
                        dhfdir(i-1, s+j-1) = dhfdir(i-1, s+j-1) + fnc * w
                        dhfdir(s+i-1, j-1) = dhfdir(s+i-1, j-1) + fnc * w

                        fnc = v * ((1.0d0 + (0.5d0 * kapa / (c * x))**2) * b0i * b0j &
                        & + 0.25d0 * (b1i * b1j - kapa / x * (b1i * b0j + b0i * b1j)) / c**2)
                        dhfdir(s+i-1, s+j-1) = dhfdir(s+i-1, s+j-1) + fnc * w
                    end do
                end do
            end do
        end do
        dhfdir = 0.5d0 * (dhfdir + transpose(dhfdir))
        return
    end subroutine dhfc_dir_matrix

    subroutine dhfc_exc_matrix(knots, order, kapa, kapb, coefb, nodes, weights, dhfexc, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: dhfexc

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1)) :: c3j
        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1), size(knots)-order,&
        & size(nodes), size(knots)-2*(order-1)-1) :: h1p, h1m, h2p, h2m
        double precision, dimension(size(knots)-order) :: coefb1, coefb2
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, ja, jb, l1min, l1max, left, right, halfd, x, w, k1, dk1, k2, dk2, b0i, b1i, gb, fb,&
        & u0, u1, v0, v1, hpj, hmj, lambda, fnc
        integer :: ndim, ier, lmin, lmax, la, lb, nknots, nintervals, nbsplines, ngl, s, l, ileft, imin, imax, m,&
        & i, j, k
        external :: drc3jj, dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        ndim = size(c3j)
        ja = 0.5d0 * (2 * abs(kapa) - 1)
        jb = 0.5d0 * (2 * abs(kapb) - 1)
        call drc3jj(ja, jb, -0.5d0, 0.5d0, l1min, l1max, c3j, ndim, ier)
        if (ier /= 0) then
            write (error_unit, *) "* ERROR in 'dhfc_exc_matrix':"
            write (error_unit, *) "* 3j coefficients could not be computed."
            stop
        end if
        lmin = nint(l1min)
        lmax = nint(l1max)

        la = (abs(2 * kapa + 1) - 1) / 2
        lb = (abs(2 * kapb + 1) - 1) / 2

        call dhfc_exc_int(knots, order, kapa, kapb, coefb, lmin, lmax, nodes, weights, h1p, h1m, h2p, h2m,&
        & speed_of_light=c)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)
        s = nbsplines - 2

        coefb1 = coefb(1:nbsplines)
        coefb2 = coefb((nbsplines+1):)

        dhfexc = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 2, x, ileft, order, dbsx, work)
                k1 = 0.0d0
                dk1 = 0.0d0
                k2 = 0.0d0
                dk2 = 0.0d0
                do i = imin, imax
                    b0i = dbsx(i+order-ileft, 0)
                    b1i = dbsx(i+order-ileft, 1)
                    k1 = k1 + coefb1(i) * b0i
                    dk1 = dk1 + coefb1(i) * b1i
                    k2 = k2 + coefb2(i) * b0i
                    dk2 = dk2 + coefb2(i) * b1i
                end do
                gb = k1 + 0.5d0 * (dk2 - kapb / x * k2) / c
                fb = k2 + 0.5d0 * (dk1 + kapb / x * k1) / c
                u0 = gb + 0.5d0 * kapa / (c * x) * fb
                u1 = 0.5d0 * fb / c
                v0 = fb - 0.5d0 * kapa / (c * x) * gb
                v1 = 0.5d0 * gb / c
                do j = 2, nbsplines - 1
                    hpj = 0.0d0
                    hmj = 0.0d0
                    do k = lmin, lmax
                        if (mod(la + lb + k, 2) /= 0) cycle
                        lambda = c3j(k)**2
                        hpj = hpj - lambda * h1p(k, j, m, l) / x**(k+1)
                        hpj = hpj - lambda * h2p(k, j, m, l) * x**k
                        hmj = hmj - lambda * h1m(k, j, m, l) / x**(k+1)
                        hmj = hmj - lambda * h2m(k, j, m, l) * x**k
                    end do
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        b1i = dbsx(i+order-ileft, 1)

                        fnc = (b0i * u0 + b1i * u1) * hpj
                        dhfexc(i-1, j-1) = dhfexc(i-1, j-1) + fnc * w

                        fnc = (b0i * u0 + b1i * u1) * hmj
                        dhfexc(i-1, s+j-1) = dhfexc(i-1, s+j-1) + fnc * w

                        fnc = (b0i * v0 + b1i * v1) * hpj
                        dhfexc(s+i-1, j-1) = dhfexc(s+i-1, j-1) + fnc * w

                        fnc = (b0i * v0 + b1i * v1) * hmj
                        dhfexc(s+i-1, s+j-1) = dhfexc(s+i-1, s+j-1) + fnc * w
                    end do
                end do
            end do
        end do
        dhfexc = 0.5d0 * (dhfexc + transpose(dhfexc))
        return
    end subroutine dhfc_exc_matrix

    subroutine dhfb_matrix(knots, order, kapa, kapb, coefb, nodes, weights, bhfmat, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: bhfmat

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1)) :: c3j
        double precision, dimension(abs(abs(kapa)-abs(kapb)):(abs(kapa)+abs(kapb)-1), size(knots)-order,&
        & size(nodes), size(knots)-2*(order-1)-1) :: h1p0, h1m0, h2p0, h2m0, h1p2, h1m2, h2p2, h2m2, h1p3, h1m3,&
        & h2p3, h2m3, h1p4, h1m4, h2p4, h2m4, h1p5, h1m5, h2p5, h2m5
        double precision, dimension(size(knots)-order) :: coefb1, coefb2
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision :: c, ja, jb, l1min, l1max, left, right, halfd, x, w, k1, dk1, k2, dk2, b0i, b1i, gb, fb,&
        & u0_hpj, u1_hpj, u0_hmj, u1_hmj, v0_hpj, v1_hpj, v0_hmj, v1_hmj, pb, qb, u0, u1, v0, v1, factor, fhpx,&
        & fhmx, fnc
        integer :: ndim, ier, lmin, lmax, la, lb, nknots, nintervals, nbsplines, ngl, s, l, ileft, imin, imax, m,&
        & i, j, k
        external :: drc3jj, dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        ndim = size(c3j)
        ja = 0.5d0 * (2 * abs(kapa) - 1)
        jb = 0.5d0 * (2 * abs(kapb) - 1)
        call drc3jj(ja, jb, -0.5d0, 0.5d0, l1min, l1max, c3j, ndim, ier)
        if (ier /= 0) then
            write (error_unit, *) "* ERROR in 'dhfb_matrix':"
            write (error_unit, *) "* 3j coefficients could not be computed."
            stop
        end if
        lmin = nint(l1min)
        lmax = nint(l1max)

        la = (abs(2 * kapa + 1) - 1) / 2
        lb = (abs(2 * kapb + 1) - 1) / 2

        call dhfb_exc_int(knots, order, kapa, kapb, coefb, lmin, lmax, nodes, weights, h1p0, h1m0, h2p0, h2m0,&
        & h1p2, h1m2, h2p2, h2m2, h1p3, h1m3, h2p3, h2m3, h1p4, h1m4, h2p4, h2m4, h1p5, h1m5, h2p5, h2m5,&
        & speed_of_light=c)

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)
        s = nbsplines - 2

        coefb1 = coefb(1:nbsplines)
        coefb2 = coefb((nbsplines+1):)

        bhfmat = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            left = knots(ileft)
            right = knots(ileft+1)
            halfd = 0.5d0 * (right - left)
            imin = ileft - order + 1
            if (imin == 1) imin = 2
            imax = ileft
            if (imax == nbsplines) imax = nbsplines - 1
            do m = 1, ngl
                x = halfd * (nodes(m) + 1.0d0) + left
                w = halfd * weights(m)
                call dbspvd(knots, order, 2, x, ileft, order, dbsx, work)
                k1 = 0.0d0
                dk1 = 0.0d0
                k2 = 0.0d0
                dk2 = 0.0d0
                do i = imin, imax
                    b0i = dbsx(i+order-ileft, 0)
                    b1i = dbsx(i+order-ileft, 1)
                    k1 = k1 + coefb1(i) * b0i
                    dk1 = dk1 + coefb1(i) * b1i
                    k2 = k2 + coefb2(i) * b0i
                    dk2 = dk2 + coefb2(i) * b1i
                end do
                gb = k1 + 0.5d0 * (dk2 - kapb / x * k2) / c
                fb = k2 + 0.5d0 * (dk1 + kapb / x * k1) / c
                do j = 2, nbsplines - 1
                    u0_hpj = 0.0d0
                    u1_hpj = 0.0d0
                    u0_hmj = 0.0d0
                    u1_hmj = 0.0d0
                    v0_hpj = 0.0d0
                    v1_hpj = 0.0d0
                    v0_hmj = 0.0d0
                    v1_hmj = 0.0d0
                    do k = lmin, lmax
                        if (k == 0) cycle
                        if (mod(la + lb + k, 2) == 0) then
                            factor = c3j(k)**2 * dble(k * (k + 1)) / (2 * k + 1)

                            pb = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb
                            qb = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = factor / (2 * k + 3) * (h1p2(k, j, m, l) / x**(k+2) &
                            & + h2p2(k, j, m, l) * x**(k+1))
                            fhmx = factor / (2 * k + 3) * (h1m2(k, j, m, l) / x**(k+2) &
                            & + h2m2(k, j, m, l) * x**(k+1))
                            u0_hpj = u0_hpj + u0 * fhpx
                            u1_hpj = u1_hpj + u1 * fhpx
                            u0_hmj = u0_hmj + u0 * fhmx
                            u1_hmj = u1_hmj + u1 * fhmx
                            v0_hpj = v0_hpj + v0 * fhpx
                            v1_hpj = v1_hpj + v1 * fhpx
                            v0_hmj = v0_hmj + v0 * fhmx
                            v1_hmj = v1_hmj + v1 * fhmx

                            pb = (dble(kapb - kapa) / k + 1.0d0) * fb
                            qb = (dble(kapb - kapa) / k - 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = factor / (2 * k - 1) * (h1p3(k, j, m, l) / x**k &
                            & + h2p3(k, j, m, l) * x**(k-1))
                            fhmx = factor / (2 * k - 1) * (h1m3(k, j, m, l) / x**k &
                            & + h2m3(k, j, m, l) * x**(k-1))
                            u0_hpj = u0_hpj + u0 * fhpx
                            u1_hpj = u1_hpj + u1 * fhpx
                            u0_hmj = u0_hmj + u0 * fhmx
                            u1_hmj = u1_hmj + u1 * fhmx
                            v0_hpj = v0_hpj + v0 * fhpx
                            v1_hpj = v1_hpj + v1 * fhpx
                            v0_hmj = v0_hmj + v0 * fhmx
                            v1_hmj = v1_hmj + v1 * fhmx

                            pb = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb
                            qb = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = 0.5d0 * factor * h1p4(k, j, m, l) / x**k
                            fhmx = 0.5d0 * factor * h1m4(k, j, m, l) / x**k
                            u0_hpj = u0_hpj - u0 * fhpx
                            u1_hpj = u1_hpj - u1 * fhpx
                            u0_hmj = u0_hmj - u0 * fhmx
                            u1_hmj = u1_hmj - u1 * fhmx
                            v0_hpj = v0_hpj - v0 * fhpx
                            v1_hpj = v1_hpj - v1 * fhpx
                            v0_hmj = v0_hmj - v0 * fhmx
                            v1_hmj = v1_hmj - v1 * fhmx

                            pb = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb
                            qb = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = 0.5d0 * factor * h1p5(k, j, m, l) / x**(k+2)
                            fhmx = 0.5d0 * factor * h1m5(k, j, m, l) / x**(k+2)
                            u0_hpj = u0_hpj + u0 * fhpx
                            u1_hpj = u1_hpj + u1 * fhpx
                            u0_hmj = u0_hmj + u0 * fhmx
                            u1_hmj = u1_hmj + u1 * fhmx
                            v0_hpj = v0_hpj + v0 * fhpx
                            v1_hpj = v1_hpj + v1 * fhpx
                            v0_hmj = v0_hmj + v0 * fhmx
                            v1_hmj = v1_hmj + v1 * fhmx

                            pb = (dble(kapb - kapa) / k + 1.0d0) * fb
                            qb = (dble(kapb - kapa) / k - 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = 0.5d0 * factor * h2p4(k, j, m, l) * x**(k-1)
                            fhmx = 0.5d0 * factor * h2m4(k, j, m, l) * x**(k-1)
                            u0_hpj = u0_hpj - u0 * fhpx
                            u1_hpj = u1_hpj - u1 * fhpx
                            u0_hmj = u0_hmj - u0 * fhmx
                            u1_hmj = u1_hmj - u1 * fhmx
                            v0_hpj = v0_hpj - v0 * fhpx
                            v1_hpj = v1_hpj - v1 * fhpx
                            v0_hmj = v0_hmj - v0 * fhmx
                            v1_hmj = v1_hmj - v1 * fhmx

                            pb = (dble(kapb - kapa) / k + 1.0d0) * fb
                            qb = (dble(kapb - kapa) / k - 1.0d0) * gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = 0.5d0 * factor * h2p5(k, j, m, l) * x**(k+1)
                            fhmx = 0.5d0 * factor * h2m5(k, j, m, l) * x**(k+1)
                            u0_hpj = u0_hpj + u0 * fhpx
                            u1_hpj = u1_hpj + u1 * fhpx
                            u0_hmj = u0_hmj + u0 * fhmx
                            u1_hmj = u1_hmj + u1 * fhmx
                            v0_hpj = v0_hpj + v0 * fhpx
                            v1_hpj = v1_hpj + v1 * fhpx
                            v0_hmj = v0_hmj + v0 * fhmx
                            v1_hmj = v1_hmj + v1 * fhmx
                        else
                            factor = c3j(k)**2 * dble(kapa + kapb)**2 / (k *&
                            & (k + 1))

                            pb = fb
                            qb = gb
                            u0 = pb + 0.5d0 * kapa / (c * x) * qb
                            u1 = 0.5d0 * qb / c
                            v0 = qb - 0.5d0 * kapa / (c * x) * pb
                            v1 = 0.5d0 * pb / c
                            fhpx = factor * (h1p0(k, j, m, l) / x**(k+1) + h2p0(k, j, m, l) * x**k)
                            fhmx = factor * (h1m0(k, j, m, l) / x**(k+1) + h2m0(k, j, m, l) * x**k)
                            u0_hpj = u0_hpj + u0 * fhpx
                            u1_hpj = u1_hpj + u1 * fhpx
                            u0_hmj = u0_hmj + u0 * fhmx
                            u1_hmj = u1_hmj + u1 * fhmx
                            v0_hpj = v0_hpj + v0 * fhpx
                            v1_hpj = v1_hpj + v1 * fhpx
                            v0_hmj = v0_hmj + v0 * fhmx
                            v1_hmj = v1_hmj + v1 * fhmx
                        end if
                    end do
                    do i = imin, imax
                        b0i = dbsx(i+order-ileft, 0)
                        b1i = dbsx(i+order-ileft, 1)

                        fnc = b0i * u0_hpj + b1i * u1_hpj
                        bhfmat(i-1, j-1) = bhfmat(i-1, j-1) + fnc * w

                        fnc = b0i * u0_hmj + b1i * u1_hmj
                        bhfmat(i-1, s+j-1) = bhfmat(i-1, s+j-1) + fnc * w

                        fnc = b0i * v0_hpj + b1i * v1_hpj
                        bhfmat(s+i-1, j-1) = bhfmat(s+i-1, j-1) + fnc * w

                        fnc = b0i * v0_hmj + b1i * v1_hmj
                        bhfmat(s+i-1, s+j-1) = bhfmat(s+i-1, s+j-1) + fnc * w
                    end do
                end do
            end do
        end do
        bhfmat = 0.5d0 * (bhfmat + transpose(bhfmat))
        return
    end subroutine dhfb_matrix

    subroutine dhfcb_matrix(knots, order, kapa, kapb, coefb, nodes, weights, dbhfmat, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: dbhfmat

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(2*(size(knots)-order-2), 2*(size(knots)-order-2)) :: dhfmat, bhfmat
        double precision :: c

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        call dhfc_matrix(knots, order, kapa, kapb, coefb, nodes, weights, dhfmat, speed_of_light=c)
        call dhfb_matrix(knots, order, kapa, kapb, coefb, nodes, weights, bhfmat, speed_of_light=c)
        dbhfmat = dhfmat + bhfmat
        return
    end subroutine dhfcb_matrix

    subroutine shf_dir_pot(knots, order, coefs, nodes, weights, vdir)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefs
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: vdir

        double precision, dimension(order, 0:0) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: v1, v2
        double precision :: i0, x0, a1, b1, d1, x1, a2, b2, d2, s2, x2, w2, u2, y2, i1
        integer :: nknots, nintervals, ngl, l, ileft, m, n, jleft, jmin, jmax
        external :: dbspvd

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        ngl = size(nodes)

        x0 = 0.0d0
        i0 = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 1, x2, jleft, order, dbsx, work)
                    u2 = sum(coefs(jmin:jmax) * dbsx(:, 0))
                    y2 = u2**2
                    s2 = s2 + y2 * w2
                end do
                i1 = i0 + s2
                v1(m, l) = i1 / x1
                i0 = i1
                x0 = x1
            end do
        end do

        x0 = knots(order + nintervals)
        i0 = 0.0d0
        do l = nintervals, 1, -1
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 1, x2, jleft, order, dbsx, work)
                    u2 = sum(coefs(jmin:jmax) * dbsx(:, 0))
                    y2 = u2**2 / x2
                    s2 = s2 + y2 * w2
                end do
                i1 = i0 + s2
                v2(m, l) = i1
                i0 = i1
                x0 = x1
            end do
        end do

        vdir = v1 + v2
        return
    end subroutine shf_dir_pot

    subroutine shf_exc_int(knots, order, coefs, lmin, lmax, nodes, weights, h1, h2)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        double precision, dimension(:), intent(in) :: coefs
        integer, intent(in) :: lmin
        integer, intent(in) :: lmax
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :, :, :), intent(out) :: h1
        double precision, dimension(:, :, :, :), intent(out) :: h2

        double precision, dimension(order, 0:0) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(lmin:lmax, size(knots)-order) :: r0, s2, r1
        double precision :: x0, a1, b1, d1, x1, a2, b2, d2, x2, w2, u2, b0j, y2
        integer :: nknots, nintervals, ngl, l, ileft, m, n, jleft, jmin, jmax, j, k
        external :: dbspvd

        h1 = 0.0d0
        h2 = 0.0d0

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        ngl = size(nodes)

        x0 = 0.0d0
        r0 = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 1, x2, jleft, order, dbsx, work)
                    u2 = sum(coefs(jmin:jmax) * dbsx(:, 0))
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        do k = lmin, lmax
                            y2 = u2 * b0j * x2**k
                            s2(k, j) = s2(k, j) + y2 * w2
                        end do
                    end do
                end do
                r1 = r0 + s2
                h1(:, :, m, l) = r1
                r0 = r1
                x0 = x1
            end do
        end do

        x0 = knots(order + nintervals)
        r0 = 0.0d0
        do l = nintervals, 1, -1
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 1, x2, jleft, order, dbsx, work)
                    u2 = sum(coefs(jmin:jmax) * dbsx(:, 0))
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        do k = lmin, lmax
                            y2 = u2 * b0j / x2**(k+1)
                            s2(k, j) = s2(k, j) + y2 * w2
                        end do
                    end do
                end do
                r1 = r0 + s2
                h2(:, :, m, l) = r1
                r0 = r1
                x0 = x1
            end do
        end do
        return
    end subroutine shf_exc_int

    subroutine dhfc_dir_pot(knots, order, kappa, coefs, nodes, weights, vdir, speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kappa
        double precision, dimension(:), intent(in) :: coefs
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :), intent(out) :: vdir

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(size(knots)-order) :: coefs1, coefs2
        double precision, dimension(size(nodes), (size(knots)-2*(order-1)-1)) :: v1, v2
        double precision :: c, i0, x0, a1, b1, d1, x1, a2, b2, d2, s2, x2, w2, k1, dk1, k2, dk2, g2, f2, y2, i1
        integer :: nknots, nbsplines, nintervals, ngl, l, ileft, m, n, jleft, jmin, jmax
        external :: dbspvd

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = size(knots) - order
        ngl = size(nodes)

        coefs1 = coefs(1:nbsplines)
        coefs2 = coefs((nbsplines+1):)

        x0 = 0.0d0
        i0 = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefs1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefs1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefs2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefs2(jmin:jmax) * dbsx(:, 1))
                    g2 = k1 + 0.5d0 * (dk2 - kappa / x2 * k2) / c
                    f2 = k2 + 0.5d0 * (dk1 + kappa / x2 * k1) / c
                    y2 = g2**2 + f2**2
                    s2 = s2 + y2 * w2
                end do
                i1 = i0 + s2
                v1(m, l) = i1 / x1
                i0 = i1
                x0 = x1
            end do
        end do

        x0 = knots(order + nintervals)
        i0 = 0.0d0
        do l = nintervals, 1, -1
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefs1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefs1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefs2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefs2(jmin:jmax) * dbsx(:, 1))
                    g2 = k1 + 0.5d0 * (dk2 - kappa / x2 * k2) / c
                    f2 = k2 + 0.5d0 * (dk1 + kappa / x2 * k1) / c
                    y2 = (g2**2 + f2**2) / x2
                    s2 = s2 + y2 * w2
                end do
                i1 = i0 + s2
                v2(m, l) = i1
                i0 = i1
                x0 = x1
            end do
        end do

        vdir = v1 + v2
        return
    end subroutine dhfc_dir_pot

    subroutine dhfc_exc_int(knots, order, kapa, kapb, coefb, lmin, lmax, nodes, weights, h1p, h1m, h2p, h2m,&
    & speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: lmin
        integer, intent(in) :: lmax
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :, :, :), intent(out) :: h1p
        double precision, dimension(:, :, :, :), intent(out) :: h1m
        double precision, dimension(:, :, :, :), intent(out) :: h2p
        double precision, dimension(:, :, :, :), intent(out) :: h2m

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(size(knots)-order) :: coefb1, coefb2
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(lmin:lmax, size(knots)-order) :: ru0, rv0, su2, sv2, ru1, rv1
        double precision :: c, x0, a1, b1, d1, x1, a2, b2, d2, x2, w2, k1, dk1, k2, dk2, gb2, fb2, u02, u12, v02,&
        & v12, b0j, b1j, y2
        integer :: nknots, nintervals, nbsplines, ngl, l, ileft, m, n, jleft, jmin, jmax, j, k
        external :: dbspvd

        h1p = 0.0d0
        h1m = 0.0d0
        h2p = 0.0d0
        h2m = 0.0d0

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)

        coefb1 = coefb(1:nbsplines)
        coefb2 = coefb((nbsplines+1):)

        x0 = 0.0d0
        ru0 = 0.0d0
        rv0 = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                su2 = 0.0d0
                sv2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefb1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefb1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefb2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefb2(jmin:jmax) * dbsx(:, 1))
                    gb2 = k1 + 0.5d0 * (dk2 - kapb / x2 * k2) / c
                    fb2 = k2 + 0.5d0 * (dk1 + kapb / x2 * k1) / c
                    u02 = gb2 + 0.5d0 * kapa / (c * x2) * fb2
                    u12 = 0.5d0 * fb2 / c
                    v02 = fb2 - 0.5d0 * kapa / (c * x2) * gb2
                    v12 = 0.5d0 * gb2 / c
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        b1j = dbsx(j+order-jleft, 1)
                        do k = lmin, lmax
                            y2 = (u02 * b0j + u12 * b1j) * x2**k
                            su2(k, j) = su2(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**k
                            sv2(k, j) = sv2(k, j) + y2 * w2
                        end do
                    end do
                end do
                ru1 = ru0 + su2
                h1p(:, :, m, l) = ru1
                ru0 = ru1
                rv1 = rv0 + sv2
                h1m(:, :, m, l) = rv1
                rv0 = rv1
                x0 = x1
            end do
        end do

        x0 = knots(order + nintervals)
        ru0 = 0.0d0
        rv0 = 0.0d0
        do l = nintervals, 1, -1
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                su2 = 0.0d0
                sv2 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefb1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefb1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefb2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefb2(jmin:jmax) * dbsx(:, 1))
                    gb2 = k1 + 0.5d0 * (dk2 - kapb / x2 * k2) / c
                    fb2 = k2 + 0.5d0 * (dk1 + kapb / x2 * k1) / c
                    u02 = gb2 + 0.5d0 * kapa / (c * x2) * fb2
                    u12 = 0.5d0 * fb2 / c
                    v02 = fb2 - 0.5d0 * kapa / (c * x2) * gb2
                    v12 = 0.5d0 * gb2 / c
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        b1j = dbsx(j+order-jleft, 1)
                        do k = lmin, lmax
                            y2 = (u02 * b0j + u12 * b1j) / x2**(k+1)
                            su2(k, j) = su2(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**(k+1)
                            sv2(k, j) = sv2(k, j) + y2 * w2
                        end do
                    end do
                end do
                ru1 = ru0 + su2
                h2p(:, :, m, l) = ru1
                ru0 = ru1
                rv1 = rv0 + sv2
                h2m(:, :, m, l) = rv1
                rv0 = rv1
                x0 = x1
            end do
        end do
        return
    end subroutine dhfc_exc_int

    subroutine dhfb_exc_int(knots, order, kapa, kapb, coefb, lmin, lmax, nodes, weights, h1p0, h1m0, h2p0, h2m0,&
    & h1p2, h1m2, h2p2, h2m2, h1p3, h1m3, h2p3, h2m3, h1p4, h1m4, h2p4, h2m4, h1p5, h1m5, h2p5, h2m5,&
    & speed_of_light)
        double precision, dimension(:), intent(in) :: knots
        integer, intent(in) :: order
        integer, intent(in) :: kapa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: coefb
        integer, intent(in) :: lmin
        integer, intent(in) :: lmax
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, dimension(:, :, :, :), intent(out) :: h1p0
        double precision, dimension(:, :, :, :), intent(out) :: h1m0
        double precision, dimension(:, :, :, :), intent(out) :: h2p0
        double precision, dimension(:, :, :, :), intent(out) :: h2m0

        double precision, dimension(:, :, :, :), intent(out) :: h1p2
        double precision, dimension(:, :, :, :), intent(out) :: h1m2
        double precision, dimension(:, :, :, :), intent(out) :: h2p2
        double precision, dimension(:, :, :, :), intent(out) :: h2m2

        double precision, dimension(:, :, :, :), intent(out) :: h1p3
        double precision, dimension(:, :, :, :), intent(out) :: h1m3
        double precision, dimension(:, :, :, :), intent(out) :: h2p3
        double precision, dimension(:, :, :, :), intent(out) :: h2m3

        double precision, dimension(:, :, :, :), intent(out) :: h1p4
        double precision, dimension(:, :, :, :), intent(out) :: h1m4
        double precision, dimension(:, :, :, :), intent(out) :: h2p4
        double precision, dimension(:, :, :, :), intent(out) :: h2m4

        double precision, dimension(:, :, :, :), intent(out) :: h1p5
        double precision, dimension(:, :, :, :), intent(out) :: h1m5
        double precision, dimension(:, :, :, :), intent(out) :: h2p5
        double precision, dimension(:, :, :, :), intent(out) :: h2m5

        double precision, intent(in), optional :: speed_of_light

        double precision, dimension(size(knots)-order) :: coefb1, coefb2
        double precision, dimension(order, 0:1) :: dbsx
        double precision, dimension((order+2)*(order+1)/2) :: work
        double precision, dimension(lmin:lmax, size(knots)-order) :: ru00, rv00, ru20, rv20, ru30, rv30, ru40,&
        & rv40, ru50, rv50, su02, sv02, su22, sv22, su32, sv32, su42, sv42, su52, sv52, ru01, rv01, ru21, rv21,&
        & ru31, rv31, ru41, rv41, ru51, rv51
        double precision :: c, x0, a1, b1, d1, x1, a2, b2, d2, x2, w2, k1, dk1, k2, dk2, gb2, fb2, b0j, b1j, pb2,&
        & qb2, u02, u12, v02, v12, y2
        integer :: nknots, nintervals, nbsplines, ngl, l, ileft, m, n, jleft, jmin, jmax, j, k
        external :: dbspvd

        h1p0 = 0.0d0
        h1m0 = 0.0d0
        h2p0 = 0.0d0
        h2m0 = 0.0d0

        h1p2 = 0.0d0
        h1m2 = 0.0d0
        h2p2 = 0.0d0
        h2m2 = 0.0d0

        h1p3 = 0.0d0
        h1m3 = 0.0d0
        h2p3 = 0.0d0
        h2m3 = 0.0d0

        h1p4 = 0.0d0
        h1m4 = 0.0d0
        h2p4 = 0.0d0
        h2m4 = 0.0d0

        h1p5 = 0.0d0
        h1m5 = 0.0d0
        h2p5 = 0.0d0
        h2m5 = 0.0d0

        c = speed_of_light_in_au
        if (present(speed_of_light)) c = speed_of_light

        nknots = size(knots)
        nintervals = nknots - 2 * (order - 1) - 1
        nbsplines = nknots - order
        ngl = size(nodes)

        coefb1 = coefb(1:nbsplines)
        coefb2 = coefb((nbsplines+1):)

        x0 = 0.0d0
        ru00 = 0.0d0
        rv00 = 0.0d0
        ru20 = 0.0d0
        rv20 = 0.0d0
        ru30 = 0.0d0
        rv30 = 0.0d0
        ru40 = 0.0d0
        rv40 = 0.0d0
        ru50 = 0.0d0
        rv50 = 0.0d0
        do l = 1, nintervals
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = 1, ngl
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                su02 = 0.0d0
                sv02 = 0.0d0
                su22 = 0.0d0
                sv22 = 0.0d0
                su32 = 0.0d0
                sv32 = 0.0d0
                su42 = 0.0d0
                sv42 = 0.0d0
                su52 = 0.0d0
                sv52 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l > 1 .and. m == 1 .and. x2 < a1) jleft = ileft - 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefb1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefb1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefb2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefb2(jmin:jmax) * dbsx(:, 1))
                    gb2 = k1 + 0.5d0 * (dk2 - kapb / x2 * k2) / c
                    fb2 = k2 + 0.5d0 * (dk1 + kapb / x2 * k1) / c
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        b1j = dbsx(j+order-jleft, 1)
                        do k = lmin, lmax
                            pb2 = fb2
                            qb2 = gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) * x2**k
                            su02(k, j) = su02(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**k
                            sv02(k, j) = sv02(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) * x2**(k+1)
                            su22(k, j) = su22(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**(k+1)
                            sv22(k, j) = sv22(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / k + 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / k - 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) * x2**(k-1)
                            su32(k, j) = su32(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**(k-1)
                            sv32(k, j) = sv32(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / k + 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / k - 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) * x2**(k-1)
                            su42(k, j) = su42(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**(k-1)
                            sv42(k, j) = sv42(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / k + 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / k - 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) * x2**(k+1)
                            su52(k, j) = su52(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) * x2**(k+1)
                            sv52(k, j) = sv52(k, j) + y2 * w2
                        end do
                    end do
                end do
                ru01 = ru00 + su02
                h1p0(:, :, m, l) = ru01
                ru00 = ru01
                rv01 = rv00 + sv02
                h1m0(:, :, m, l) = rv01
                rv00 = rv01

                ru21 = ru20 + su22
                h1p2(:, :, m, l) = ru21
                ru20 = ru21
                rv21 = rv20 + sv22
                h1m2(:, :, m, l) = rv21
                rv20 = rv21

                ru31 = ru30 + su32
                h1p3(:, :, m, l) = ru31
                ru30 = ru31
                rv31 = rv30 + sv32
                h1m3(:, :, m, l) = rv31
                rv30 = rv31

                ru41 = ru40 + su42
                h1p4(:, :, m, l) = ru41
                ru40 = ru41
                rv41 = rv40 + sv42
                h1m4(:, :, m, l) = rv41
                rv40 = rv41

                ru51 = ru50 + su52
                h1p5(:, :, m, l) = ru51
                ru50 = ru51
                rv51 = rv50 + sv52
                h1m5(:, :, m, l) = rv51
                rv50 = rv51

                x0 = x1
            end do
        end do

        x0 = knots(order + nintervals)
        ru00 = 0.0d0
        rv00 = 0.0d0
        ru20 = 0.0d0
        rv20 = 0.0d0
        ru30 = 0.0d0
        rv30 = 0.0d0
        ru40 = 0.0d0
        rv40 = 0.0d0
        ru50 = 0.0d0
        rv50 = 0.0d0
        do l = nintervals, 1, -1
            ileft = order - 1 + l
            a1 = knots(ileft)
            b1 = knots(ileft+1)
            d1 = 0.5d0 * (b1 - a1)
            do m = ngl, 1, -1
                x1 = d1 * (nodes(m) + 1.0d0) + a1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                su02 = 0.0d0
                sv02 = 0.0d0
                su22 = 0.0d0
                sv22 = 0.0d0
                su32 = 0.0d0
                sv32 = 0.0d0
                su42 = 0.0d0
                sv42 = 0.0d0
                su52 = 0.0d0
                sv52 = 0.0d0
                do n = ngl, 1, -1
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    w2 = d2 * weights(n)
                    jleft = ileft
                    if (l < nintervals .and. m == ngl .and. x2 >= b1) jleft = ileft + 1
                    jmin = jleft - order + 1
                    jmax = jleft
                    call dbspvd(knots, order, 2, x2, jleft, order, dbsx, work)
                    k1 = sum(coefb1(jmin:jmax) * dbsx(:, 0))
                    dk1 = sum(coefb1(jmin:jmax) * dbsx(:, 1))
                    k2 = sum(coefb2(jmin:jmax) * dbsx(:, 0))
                    dk2 = sum(coefb2(jmin:jmax) * dbsx(:, 1))
                    gb2 = k1 + 0.5d0 * (dk2 - kapb / x2 * k2) / c
                    fb2 = k2 + 0.5d0 * (dk1 + kapb / x2 * k1) / c
                    do j = jmin, jmax
                        b0j = dbsx(j+order-jleft, 0)
                        b1j = dbsx(j+order-jleft, 1)
                        do k = lmin, lmax
                            pb2 = fb2
                            qb2 = gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) / x2**(k+1)
                            su02(k, j) = su02(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**(k+1)
                            sv02(k, j) = sv02(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) / x2**(k+2)
                            su22(k, j) = su22(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**(k+2)
                            sv22(k, j) = sv22(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / k + 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / k - 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) / x2**k
                            su32(k, j) = su32(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**k
                            sv32(k, j) = sv32(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) / x2**k
                            su42(k, j) = su42(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**k
                            sv42(k, j) = sv42(k, j) + y2 * w2

                            pb2 = (dble(kapb - kapa) / (k + 1) - 1.0d0) * fb2
                            qb2 = (dble(kapb - kapa) / (k + 1) + 1.0d0) * gb2
                            u02 = pb2 + 0.5d0 * kapa / (c * x2) * qb2
                            u12 = 0.5d0 * qb2 / c
                            v02 = qb2 - 0.5d0 * kapa / (c * x2) * pb2
                            v12 = 0.5d0 * pb2 / c
                            y2 = (u02 * b0j + u12 * b1j) / x2**(k+2)
                            su52(k, j) = su52(k, j) + y2 * w2
                            y2 = (v02 * b0j + v12 * b1j) / x2**(k+2)
                            sv52(k, j) = sv52(k, j) + y2 * w2
                        end do
                    end do
                end do
                ru01 = ru00 + su02
                h2p0(:, :, m, l) = ru01
                ru00 = ru01
                rv01 = rv00 + sv02
                h2m0(:, :, m, l) = rv01
                rv00 = rv01

                ru21 = ru20 + su22
                h2p2(:, :, m, l) = ru21
                ru20 = ru21
                rv21 = rv20 + sv22
                h2m2(:, :, m, l) = rv21
                rv20 = rv21

                ru31 = ru30 + su32
                h2p3(:, :, m, l) = ru31
                ru30 = ru31
                rv31 = rv30 + sv32
                h2m3(:, :, m, l) = rv31
                rv30 = rv31

                ru41 = ru40 + su42
                h2p4(:, :, m, l) = ru41
                ru40 = ru41
                rv41 = rv40 + sv42
                h2m4(:, :, m, l) = rv41
                rv40 = rv41

                ru51 = ru50 + su52
                h2p5(:, :, m, l) = ru51
                ru50 = ru51
                rv51 = rv50 + sv52
                h2m5(:, :, m, l) = rv51
                rv50 = rv51

                x0 = x1
            end do
        end do
        return
    end subroutine dhfb_exc_int

end module bspwf_hartree_fock
