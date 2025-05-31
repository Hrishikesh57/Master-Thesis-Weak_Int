module bspqed_integrals
    use bspqed_dkb, only: vfunc, wfunc, pfunc, qfunc
    use bspqed_function, only: function_type
    implicit none
    private

    public :: intlv
    public :: intlw
    public :: intlvx

    public :: intif_re
    public :: intic_re
    public :: inti1_re
    public :: inti1_st
    public :: inti1_dy_re
    public :: inti2_re
    public :: inti2_st
    public :: inti2_dy_re
    public :: inti3_re
    public :: inti3_st
    public :: inti3_dy_re
    public :: inti4_re
    public :: inti4_st
    public :: inti4_dy_re
    public :: intji_re
    public :: intji_st
    public :: intji_dy_re

    public :: intif_re_drv
    public :: intic_re_drv
    public :: inti1_re_drv
    public :: inti1_dy_re_drv
    public :: inti2_re_drv
    public :: inti2_dy_re_drv
    public :: inti3_re_drv
    public :: inti3_dy_re_drv
    public :: inti4_re_drv
    public :: inti4_dy_re_drv
    public :: intji_re_drv
    public :: intji_dy_re_drv

    public :: intif_re_drv_2
    public :: intic_re_drv_2
    public :: inti1_re_drv_2
    public :: inti1_st_drv_2
    public :: inti1_dy_re_drv_2
    public :: inti2_re_drv_2
    public :: inti2_st_drv_2
    public :: inti2_dy_re_drv_2
    public :: inti3_re_drv_2
    public :: inti3_st_drv_2
    public :: inti3_dy_re_drv_2
    public :: inti4_re_drv_2
    public :: inti4_st_drv_2
    public :: inti4_dy_re_drv_2
    public :: intji_re_drv_2
    public :: intji_st_drv_2
    public :: intji_dy_re_drv_2

    public :: intif_ia
    public :: inti1_ia
    public :: inti1_dy_ia
    public :: inti2_ia
    public :: inti2_dy_ia
    public :: inti3_ia
    public :: inti3_dy_ia
    public :: intji_ia
    public :: intji_dy_ia

    public :: intif_ia_drv
    public :: inti1_ia_drv
    public :: inti1_dy_ia_drv
    public :: inti2_ia_drv
    public :: inti2_dy_ia_drv
    public :: inti3_ia_drv
    public :: inti3_dy_ia_drv
    public :: intji_ia_drv
    public :: intji_dy_ia_drv

    double precision, parameter :: pihalf = 2.0d0 * atan(1.0d0)

    contains

    double precision function intlv(knt, ord, kapa, cfa, kapb, cfb, potential, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x, s, t, a, b, d, pot, vab
        integer :: ngl, npoints, nintervals, il, iv, i, j

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1
        s = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a = knt(il)
            b = knt(il+1)
            d = 0.5d0 * (b - a)
            t = 0.0d0
            do j = 1, ngl
                x = d * (nodes(j) + 1.0d0) + a
                pot = potential%eval(x)
                vab = vfunc(knt, ord, kapa, cfa, kapb, cfb, x, iv)
                t = t + pot * vab * weights(j)
            enddo
            t = d * t
            s = s + t
        enddo
        intlv = s
        return
    end function intlv

    double precision function intlw(knt, ord, kapa, cfa, kapb, cfb, potential, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        class(function_type), intent(in) :: potential
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x, s, t, a, b, d, pot, wab
        integer :: ngl, npoints, nintervals, il, iv, i, j

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1
        s = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a = knt(il)
            b = knt(il+1)
            d = 0.5d0 * (b - a)
            t = 0.0d0
            do j = 1, ngl
                x = d * (nodes(j) + 1.0d0) + a
                pot = potential%eval(x)
                wab = wfunc(knt, ord, kapa, cfa, kapb, cfb, x, iv)
                t = t + pot * wab * weights(j)
            enddo
            t = d * t
            s = s + t
        enddo
        intlw = s
        return
    end function intlw

    double precision function intlvx(knt, ord, kapa, cfa, kapb, cfb, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x, s, t, a, b, d, vab
        integer :: ngl, npoints, nintervals, il, iv, i, j

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1
        s = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a = knt(il)
            b = knt(il+1)
            d = 0.5d0 * (b - a)
            t = 0.0d0
            do j = 1, ngl
                x = d * (nodes(j) + 1.0d0) + a
                vab = vfunc(knt, ord, kapa, cfa, kapb, cfb, x, iv)
                t = t + x * vab * weights(j)
            enddo
            t = d * t
            s = s + t
        enddo
        intlvx = s
        return
    end function intlvx

    double precision function intif_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3

        i1 = inti1_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        i2 = inti2_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        intif_re = i1 + (k + 1) * i2 / (2 * k + 3)
        if (k /= 0) then
            i3 = inti3_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            intif_re = intif_re + k * i3 / (2 * k - 1)
        end if
        return
    end function intif_re

    double precision function intic_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3, i4, i5

        i1 = inti1_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        intic_re = i1
        if (k /= 0) then
            i2 = inti2_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i3 = inti3_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i4 = inti4_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i5 = inti4_re(knt, ord, kapb, cfb, kapa, cfa, kapd, cfd, kapc, cfc, k, omg, nodes, weights)
            intic_re = intic_re + k * (k + 1) * (i2 / (2 * k + 3) + i3 / (2 * k - 1) - i4 - i5) / (2 * k + 1)
        end if
        return
    end function intic_re

    double precision function inti1_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti1_re = inti1_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti1_re = inti1_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti1_re

    double precision function inti1_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f1 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f2 = f2 * x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f1 = f1 * x1**k
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f2 = f2 / x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti1_st = i1 + i2
        return
    end function inti1_st

    double precision function inti1_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesy_sph(k, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesy_sph(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti1_dy_re = -omg * (2 * k + 1) * (i1 + i2)
        return
    end function inti1_dy_re

    double precision function inti2_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti2_re = inti2_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti2_re = inti2_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti2_re

    double precision function inti2_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f1 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f2 = f2 * x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f1 = f1 * x1**(k+1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f2 = f2 / x2**(k+2)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti2_st = i1 + i2
        return
    end function inti2_st

    double precision function inti2_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesy_sph(k+1, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k+1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesy_sph(k+1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti2_dy_re = -omg * (2 * k + 3) * (i1 + i2)
        return
    end function inti2_dy_re

    double precision function inti3_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti3_re = inti3_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti3_re = inti3_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti3_re

    double precision function inti3_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f1 = f1 / x1**k
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f2 = f2 * x2**(k-1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f1 = f1 * x1**(k-1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f2 = f2 / x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti3_st = i1 + i2
        return
    end function inti3_st

    double precision function inti3_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesy_sph(k-1, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k-1, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesy_sph(k-1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti3_dy_re = -omg * (2 * k - 1) * (i1 + i2)
        return
    end function inti3_dy_re

    double precision function inti4_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg < 1.0d-08) then
            inti4_re = inti4_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti4_re = inti4_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti4_re

    double precision function inti4_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, x0, h01, h02, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22, x2&
        &, f2, f21, f22, h11, h12
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        i2 = 0.0d0
        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 / x1**k
                f12 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 * x2**(k-1)
                    f22 = f2 * x2**(k+1)
                    s21 = s21 + f21 * weights(n)
                    s22 = s22 + f22 * weights(n)
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                s11 = s11 + f11 * h11 * weights(j)
                s12 = s12 + f12 * h12 * weights(j)
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i1 = i1 + s11
            i2 = i2 + s12
        end do

        inti4_st = 0.5d0 * (i1 - i2)
        return
    end function inti4_st

    double precision function inti4_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, x0, h01, h02, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21, s22&
        &, x2, p2, q2, f21, f22, h11, h12, i3
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        i2 = 0.0d0
        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesy_sph(k+1, omg * x1)
                f11 = p1 * q1
                f12 = p1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    f22 = p2 * x2**(k-1)
                    s21 = s21 + f21 * weights(n)
                    s22 = s22 + f22 * weights(n)
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                s11 = s11 + f11 * h11 * weights(j)
                s12 = s12 + f12 * h12 * weights(j)
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i1 = i1 + s11
            i2 = i2 + s12
        end do

        i3 = 0.0d0
        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f11 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    s21 = s21 + f21 * weights(n)
                end do
                s21 = d2 * s21
                h11 = h01 + s21
                s11 = s11 + f11 * h11 * weights(j)
                x0 = x1
                h01 = h11
            end do
            s11 = d1 * s11
            i3 = i3 + s11
        end do

        inti4_dy_re = -omg * (i1 + i3) - (2 * k + 1) * i2 / omg**2
        return
    end function inti4_dy_re

    double precision function intji_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            intji_re = intji_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            intji_re = intji_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function intji_re

    double precision function intji_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, f1, a2, b2, d2, s2, x2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f1 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f2 = f2 * x2**k
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f1 = f1 * x1**k
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f2 = f2 / x2**(k+1)
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        intji_st = i1 + i2
        return
    end function intji_st

    double precision function intji_dy_re(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesy_sph(k, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesy_sph(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        intji_dy_re = -omg * (2 * k + 1) * (i1 + i2)
        return
    end function intji_dy_re

    double precision function intif_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3

        i1 = inti1_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        i2 = inti2_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        intif_re_drv = i1 + (k + 1) * i2 / (2 * k + 3)
        if (k /= 0) then
            i3 = inti3_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            intif_re_drv = intif_re_drv + k * i3 / (2 * k - 1)
        end if
        return
    end function intif_re_drv

    double precision function intic_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i2, i3, i4, i5

        intic_re_drv = 0.0d0
        if (k /= 0) then
            i2 = inti2_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i3 = inti3_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i4 = inti4_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i5 = inti4_re_drv(knt, ord, kapb, cfb, kapa, cfa, kapd, cfd, kapc, cfc, k, omg, nodes, weights)
            intic_re_drv = intic_re_drv + k * (k + 1) * (i2 / (2 * k + 3) + i3 / (2 * k - 1) - i4 - i5) / (2 * k +&
            & 1)
        end if
        return
    end function intic_re_drv

    double precision function inti1_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti1_re_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti1_re_drv = inti1_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti1_re_drv

    double precision function inti1_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+1, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti1_dy_re_drv = -(2 * k + 1) * (i1 - i2)
        return
    end function inti1_dy_re_drv

    double precision function inti2_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti2_re_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti2_re_drv = inti2_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti2_re_drv

    double precision function inti2_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesy_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k+1, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k+1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+2, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+2, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti2_dy_re_drv = -(2 * k + 3) * (i1 - i2)
        return
    end function inti2_dy_re_drv

    double precision function inti3_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti3_re_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti3_re_drv = inti3_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti3_re_drv

    double precision function inti3_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-2, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k-1, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-2, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k-1, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti3_dy_re_drv = -(2 * k - 1) * (i1 - i2)
        return
    end function inti3_dy_re_drv

    double precision function inti4_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti4_re_drv = 0.0d0
        if (omg >= 1.0d-08) then
            inti4_re_drv = inti4_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti4_re_drv

    double precision function inti4_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, i1, i2, i3, i4, a1, b1, d1, s11, s12, s13, s14, x1, p1, q1,&
        & f11, f12, f13, f14, a2, b2, d2, s21, s22, s23, s24, x2, p2, q2, f21, f22, f23, f24, wgt, h11, h12, h13,&
        & h14, i5, i6, i7
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        i1 = 0.0d0
        i2 = 0.0d0
        i3 = 0.0d0
        i4 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesy_sph(k+1, omg * x1)
                f11 = p1 * q1
                f12 = f11
                q1 = omg * x1 * dbesy_sph(k+2, omg * x1)
                f13 = p1 * q1
                f14 = p1 / (omg * x1)**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k, omg * x2)
                    f22 = p2 * q2
                    f23 = f21
                    f24 = p2 * (omg * x2)**(k-1)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            i1 = i1 + s11
            i2 = i2 + s12
            i3 = i3 + s13
            i4 = i4 + s14
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        i5 = 0.0d0
        i6 = 0.0d0
        i7 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+2, omg * x1)
                f12 = p1 * q1
                f13 = f11
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    f22 = f21
                    q2 = omg * x2 * dbesy_sph(k, omg * x2)
                    f23 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            i5 = i5 + s11
            i6 = i6 + s12
            i7 = i7 + s13
        end do

        inti4_dy_re_drv = -(2 * k + 1) * (i1 - 2.0d0 * i4 + i5) + i2 + i3 + i6 + i7
        return
    end function inti4_dy_re_drv

    double precision function intji_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        intji_re_drv = 0.0d0
        if (omg /= 0.0d0) then
            intji_re_drv = intji_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function intji_re_drv

    double precision function intji_dy_re_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+1, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k, omg * x2)
                    f22 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        intji_dy_re_drv = -(2 * k + 1) * (i1 - i2)
        return
    end function intji_dy_re_drv

    double precision function intif_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3

        i1 = inti1_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        i2 = inti2_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        intif_re_drv_2 = i1 + (k + 1) * i2 / (2 * k + 3)
        if (k /= 0) then
            i3 = inti3_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            intif_re_drv_2 = intif_re_drv_2 + k * i3 / (2 * k - 1)
        end if
        return
    end function intif_re_drv_2

    double precision function intic_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i2, i3, i4, i5

        intic_re_drv_2 = 0.0d0
        if (k /= 0) then
            i2 = inti2_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i3 = inti3_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i4 = inti4_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            i5 = inti4_re_drv_2(knt, ord, kapb, cfb, kapa, cfa, kapd, cfd, kapc, cfc, k, omg, nodes, weights)
            intic_re_drv_2 = intic_re_drv_2 + k * (k + 1) * (i2 / (2 * k + 3) + i3 / (2 * k - 1) - i4 - i5) / (2 *&
            & k + 1)
        end if
        return
    end function intic_re_drv_2

    double precision function inti1_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti1_re_drv_2 = inti1_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
            & weights)
        else
            inti1_re_drv_2 = inti1_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes&
            &, weights)
        end if
        return
    end function inti1_re_drv_2

    double precision function inti1_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22,&
        & x2, f2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f11 = f1 / x1**(k-1)
                f12 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f21 = f2 * x2**k
                    f22 = f2 * x2**(k+2)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f11 = f1 * x1**k
                f12 = f1 * x1**(k+2)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f21 = f2 / x2**(k-1)
                    f22 = f2 / x2**(k+1)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti1_st_drv_2 = i1 / (2 * k - 1) - i2 / (2 * k + 3)
        return
    end function inti1_st_drv_2

    double precision function inti1_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, h05, i11, i21, i31, i41, i51, a1, b1, d1, s11, s12, s13, s14,&
        & s15, x1, p1, q1, f11, f12, f13, f14, f15, a2, b2, d2, s21, s22, s23, s24, s25, x2, p2, q2, f21, f22, f23&
        &, f24, f25, wgt, h11, h12, h13, h14, h15, i12, i22, i32, i42, i52, i1, i2, i3, i4, i5
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        i31 = 0.0d0
        i41 = 0.0d0
        i51 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k, omg * x1)
                f12 = p1 * q1
                f13 = f12
                f14 = f12 * (omg * x1)**2
                q1 = omg * x1 * dbesy_sph(k+1, omg * x1)
                f15 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    f23 = f21 * (omg * x2)**2
                    f24 = f21
                    f25 = f22
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i11 = i11 + s11
            i21 = i21 + s12
            i31 = i31 + s13
            i41 = i41 + s14
            i51 = i51 + s15
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        i32 = 0.0d0
        i42 = 0.0d0
        i52 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+1, omg * x1)
                f12 = p1 * q1
                f13 = f11 * (omg * x1)**2
                f14 = f11
                f15 = f12
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k, omg * x2)
                    f22 = p2 * q2
                    f23 = f22
                    f24 = f22 * (omg * x2)**2
                    q2 = omg * x2 * dbesy_sph(k+1, omg * x2)
                    f25 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i12 = i12 + s11
            i22 = i22 + s12
            i32 = i32 + s13
            i42 = i42 + s14
            i52 = i52 + s15
        end do

        i1 = i11 + i12
        i2 = i21 + i22
        i3 = i31 + i32
        i4 = i41 + i42
        i5 = i51 + i52

        inti1_dy_re_drv_2 = -(2 * k + 1) * (2 * k * (i1 - i2) - i3 - i4 + 2.0d0 * i5) / omg
        return
    end function inti1_dy_re_drv_2

    double precision function inti2_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti2_re_drv_2 = inti2_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
            & weights)
        else
            inti2_re_drv_2 = inti2_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes&
            &, weights)
        end if
        return
    end function inti2_re_drv_2

    double precision function inti2_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22,&
        & x2, f2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 / x1**k
                f12 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 * x2**(k+1)
                    f22 = f2 * x2**(k+3)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 * x1**(k+1)
                f12 = f1 * x1**(k+3)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 / x2**k
                    f22 = f2 / x2**(k+2)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti2_st_drv_2 = i1 / (2 * k + 1) - i2 / (2 * k + 5)
        return
    end function inti2_st_drv_2

    double precision function inti2_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, h05, i11, i21, i31, i41, i51, a1, b1, d1, s11, s12, s13, s14,&
        & s15, x1, p1, q1, f11, f12, f13, f14, f15, a2, b2, d2, s21, s22, s23, s24, s25, x2, p2, q2, f21, f22, f23&
        &, f24, f25, wgt, h11, h12, h13, h14, h15, i12, i22, i32, i42, i52, i1, i2, i3, i4, i5
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        i31 = 0.0d0
        i41 = 0.0d0
        i51 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesy_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k+1, omg * x1)
                f12 = p1 * q1
                f13 = f12
                f14 = f12 * (omg * x1)**2
                q1 = omg * x1 * dbesy_sph(k+2, omg * x1)
                f15 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k+1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+2, omg * x2)
                    f22 = p2 * q2
                    f23 = f21 * (omg * x2)**2
                    f24 = f21
                    f25 = f22
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i11 = i11 + s11
            i21 = i21 + s12
            i31 = i31 + s13
            i41 = i41 + s14
            i51 = i51 + s15
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        i32 = 0.0d0
        i42 = 0.0d0
        i52 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+2, omg * x1)
                f12 = p1 * q1
                f13 = f11 * (omg * x1)**2
                f14 = f11
                f15 = f12
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    f23 = f22
                    f24 = f22 * (omg * x2)**2
                    q2 = omg * x2 * dbesy_sph(k+2, omg * x2)
                    f25 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i12 = i12 + s11
            i22 = i22 + s12
            i32 = i32 + s13
            i42 = i42 + s14
            i52 = i52 + s15
        end do

        i1 = i11 + i12
        i2 = i21 + i22
        i3 = i31 + i32
        i4 = i41 + i42
        i5 = i51 + i52

        inti2_dy_re_drv_2 = -(2 * k + 3) * (2 * (k + 1) * (i1 - i2) - i3 - i4 + 2.0d0 * i5) / omg
        return
    end function inti2_dy_re_drv_2

    double precision function inti3_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti3_re_drv_2 = inti3_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
            & weights)
        else
            inti3_re_drv_2 = inti3_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes&
            &, weights)
        end if
        return
    end function inti3_re_drv_2

    double precision function inti3_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22,&
        & x2, f2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 / x1**(k-2)
                f12 = f1 / x1**k
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 * x2**(k-1)
                    f22 = f2 * x2**(k+1)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 * x1**(k-1)
                f12 = f1 * x1**(k+1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 / x2**(k-2)
                    f22 = f2 / x2**k
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti3_st_drv_2 = i1 / (2 * k - 3) - i2 / (2 * k + 1)
        return
    end function inti3_st_drv_2

    double precision function inti3_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, h05, i11, i21, i31, i41, i51, a1, b1, d1, s11, s12, s13, s14,&
        & s15, x1, p1, q1, f11, f12, f13, f14, f15, a2, b2, d2, s21, s22, s23, s24, s25, x2, p2, q2, f21, f22, f23&
        &, f24, f25, wgt, h11, h12, h13, h14, h15, i12, i22, i32, i42, i52, i1, i2, i3, i4, i5
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        i31 = 0.0d0
        i41 = 0.0d0
        i51 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-2, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k-1, omg * x1)
                f12 = p1 * q1
                f13 = f12
                f14 = f12 * (omg * x1)**2
                q1 = omg * x1 * dbesy_sph(k, omg * x1)
                f15 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k, omg * x2)
                    f22 = p2 * q2
                    f23 = f21 * (omg * x2)**2
                    f24 = f21
                    f25 = f22
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i11 = i11 + s11
            i21 = i21 + s12
            i31 = i31 + s13
            i41 = i41 + s14
            i51 = i51 + s15
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        i32 = 0.0d0
        i42 = 0.0d0
        i52 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k, omg * x1)
                f12 = p1 * q1
                f13 = f11 * (omg * x1)**2
                f14 = f11
                f15 = f12
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-2, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k-1, omg * x2)
                    f22 = p2 * q2
                    f23 = f22
                    f24 = f22 * (omg * x2)**2
                    q2 = omg * x2 * dbesy_sph(k, omg * x2)
                    f25 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i12 = i12 + s11
            i22 = i22 + s12
            i32 = i32 + s13
            i42 = i42 + s14
            i52 = i52 + s15
        end do

        i1 = i11 + i12
        i2 = i21 + i22
        i3 = i31 + i32
        i4 = i41 + i42
        i5 = i51 + i52

        inti3_dy_re_drv_2 = -(2 * k - 1) * (2 * (k - 1) * (i1 - i2) - i3 - i4 + 2.0d0 * i5) / omg
        return
    end function inti3_dy_re_drv_2

    double precision function inti4_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg < 1.0d-08) then
            inti4_re_drv_2 = inti4_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
            & weights)
        else
            inti4_re_drv_2 = inti4_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes&
            &, weights)
        end if
        return
    end function inti4_re_drv_2

    double precision function inti4_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, i1, i2, i3, a1, b1, d1, s11, s12, s13, x1, f1, f11, f12, f13, a2,&
        & b2, d2, s21, s22, s23, x2, f2, f21, f22, f23, wgt, h11, h12, h13, i4
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        i1 = 0.0d0
        i2 = 0.0d0
        i3 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 / x1**(k-2)
                f12 = f1 / x1**k
                f13 = f1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 * x2**(k-1)
                    f22 = f2 * x2**(k+1)
                    f23 = f2 * x2**(k+3)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            i1 = i1 + s11
            i2 = i2 + s12
            i3 = i3 + s13
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        i4 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                f11 = f1 * x1**(k+1)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    f21 = f2 / x2**k
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                end do
                s21 = d2 * s21
                h11 = h01 + s21
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                x0 = x1
                h01 = h11
            end do
            s11 = d1 * s11
            i4 = i4 + s11
        end do

        inti4_st_drv_2 = 0.25d0 * (i1 / (2 * k - 1) - 2 * i2 / (2 * k + 1) + i3 / (2 * k + 3) + i4 / ((2 * k + 3) &
        &* (2 * k + 1) * (2 * k - 1)))
        return
    end function inti4_st_drv_2

    double precision function inti4_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, h05, h06, h07, i1, i2, i3, i4, i5, i6, i7, a1, b1, d1, s11,&
        & s12, s13, s14, s15, s16, s17, x1, p1, q1, f11, f12, f13, f14, f15, f16, f17, a2, b2, d2, s21, s22, s23,&
        & s24, s25, s26, s27, x2, p2, q2, f21, f22, f23, f24, f25, f26, f27, wgt, h11, h12, h13, h14, h15, h16,&
        & h17, s, t
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        h06 = 0.0d0
        h07 = 0.0d0
        i1 = 0.0d0
        i2 = 0.0d0
        i3 = 0.0d0
        i4 = 0.0d0
        i5 = 0.0d0
        i6 = 0.0d0
        i7 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            s16 = 0.0d0
            s17 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesy_sph(k+1, omg * x1)
                f11 = p1 * q1
                f12 = f11
                q1 = omg * x1 * dbesy_sph(k+2, omg * x1)
                f13 = p1 * q1
                f14 = f11
                f15 = f11 * (omg * x1)**2
                f16 = f13
                f17 = p1 / x1**(k+2)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                s26 = 0.0d0
                s27 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesj_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k, omg * x2)
                    f22 = p2 * q2
                    f23 = f21
                    f24 = f21 * (omg * x2)**2
                    f25 = f21
                    f26 = f22
                    f27 = p2 * x2**(k-1)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                    s26 = s26 + f26 * wgt
                    s27 = s27 + f27 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                s26 = d2 * s26
                s27 = d2 * s27
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                h16 = h06 + s26
                h17 = h07 + s27
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                s16 = s16 + f16 * h16 * wgt
                s17 = s17 + f17 * h17 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
                h06 = h16
                h07 = h17
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            s16 = d1 * s16
            s17 = d1 * s17
            i1 = i1 + s11
            i2 = i2 + s12
            i3 = i3 + s13
            i4 = i4 + s14
            i5 = i5 + s15
            i6 = i6 + s16
            i7 = i7 + s17
        end do
        s = -((2 * k + 1) * 2 * k * i1 - 2 * (k + 1) * i2 - 2 * (k - 1) * i3 - i4 - i5 + 2.0d0 * i6) / omg - 6.0d0&
        & * (2 * k + 1) * i7 / omg**4

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        h06 = 0.0d0
        i1 = 0.0d0
        i2 = 0.0d0
        i3 = 0.0d0
        i4 = 0.0d0
        i5 = 0.0d0
        i6 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            s16 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesj_sph(k+1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+2, omg * x1)
                f12 = p1 * q1
                f13 = f11
                f14 = f11 * (omg * x1)**2
                f15 = f11
                f16 = f12
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                s26 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    f22 = f21
                    q2 = omg * x2 * dbesy_sph(k, omg * x2)
                    f23 = p2 * q2
                    f24 = f21
                    f25 = f21 * (omg * x2)**2
                    f26 = f23
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                    s26 = s26 + f26 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                s26 = d2 * s26
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                h16 = h06 + s26
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                s16 = s16 + f16 * h16 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
                h06 = h16
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            s16 = d1 * s16
            i1 = i1 + s11
            i2 = i2 + s12
            i3 = i3 + s13
            i4 = i4 + s14
            i5 = i5 + s15
            i6 = i6 + s16
        end do
        t = -((2 * k + 1) * 2 * k * i1 - 2 * (k - 1) * i2 - 2 * (k + 1) * i3 - i4 - i5 + 2.0d0 * i6) / omg

        inti4_dy_re_drv_2 = s + t
        return
    end function inti4_dy_re_drv_2

    double precision function intji_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            intji_re_drv_2 = intji_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
            & weights)
        else
            intji_re_drv_2 = intji_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes&
            &, weights)
        end if
        return
    end function intji_re_drv_2

    double precision function intji_st_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, f1, f11, f12, a2, b2, d2, s21, s22,&
        & x2, f2, f21, f22, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f11 = f1 / x1**(k-1)
                f12 = f1 / x1**(k+1)
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f21 = f2 * x2**k
                    f22 = f2 * x2**(k+2)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                f1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                f11 = f1 * x1**k
                f12 = f1 * x1**(k+2)
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    f2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    f21 = f2 / x2**(k-1)
                    f22 = f2 / x2**(k+1)
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                h11 = h01 + s21
                h12 = h02 + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        intji_st_drv_2 = i1 / (2 * k - 1) - i2 / (2 * k + 3)
        return
    end function intji_st_drv_2

    double precision function intji_dy_re_drv_2(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg,&
    & nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, h03, h04, h05, i11, i21, i31, i41, i51, a1, b1, d1, s11, s12, s13, s14,&
        & s15, x1, p1, q1, f11, f12, f13, f14, f15, a2, b2, d2, s21, s22, s23, s24, s25, x2, p2, q2, f21, f22, f23&
        &, f24, f25, wgt, h11, h12, h13, h14, h15, i12, i22, i32, i42, i52, i1, i2, i3, i4, i5
        integer :: ngl, npoints, nintervals, i, j, il, iv, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        i31 = 0.0d0
        i41 = 0.0d0
        i51 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesy_sph(k-1, omg * x1)
                f11 = p1 * q1
                q1 = dbesy_sph(k, omg * x1)
                f12 = p1 * q1
                f13 = f12
                f14 = f12 * (omg * x1)**2
                q1 = omg * x1 * dbesy_sph(k+1, omg * x1)
                f15 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesj_sph(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesj_sph(k+1, omg * x2)
                    f22 = p2 * q2
                    f23 = f21 * (omg * x2)**2
                    f24 = f21
                    f25 = f22
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i11 = i11 + s11
            i21 = i21 + s12
            i31 = i31 + s13
            i41 = i41 + s14
            i51 = i51 + s15
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        h03 = 0.0d0
        h04 = 0.0d0
        h05 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        i32 = 0.0d0
        i42 = 0.0d0
        i52 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            s13 = 0.0d0
            s14 = 0.0d0
            s15 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesj_sph(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesj_sph(k+1, omg * x1)
                f12 = p1 * q1
                f13 = f11 * (omg * x1)**2
                f14 = f11
                f15 = f12
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                s23 = 0.0d0
                s24 = 0.0d0
                s25 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesy_sph(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesy_sph(k, omg * x2)
                    f22 = p2 * q2
                    f23 = f22
                    f24 = f22 * (omg * x2)**2
                    q2 = omg * x2 * dbesy_sph(k+1, omg * x2)
                    f25 = p2 * q2
                    wgt = weights(n)
                    s21 = s21 + f21 * wgt
                    s22 = s22 + f22 * wgt
                    s23 = s23 + f23 * wgt
                    s24 = s24 + f24 * wgt
                    s25 = s25 + f25 * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                s23 = d2 * s23
                s24 = d2 * s24
                s25 = d2 * s25
                h11 = h01 + s21
                h12 = h02 + s22
                h13 = h03 + s23
                h14 = h04 + s24
                h15 = h05 + s25
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                s13 = s13 + f13 * h13 * wgt
                s14 = s14 + f14 * h14 * wgt
                s15 = s15 + f15 * h15 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
                h03 = h13
                h04 = h14
                h05 = h15
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            s13 = d1 * s13
            s14 = d1 * s14
            s15 = d1 * s15
            i12 = i12 + s11
            i22 = i22 + s12
            i32 = i32 + s13
            i42 = i42 + s14
            i52 = i52 + s15
        end do

        i1 = i11 + i12
        i2 = i21 + i22
        i3 = i31 + i32
        i4 = i41 + i42
        i5 = i51 + i52

        intji_dy_re_drv_2 = -(2 * k + 1) * (2 * k * (i1 - i2) - i3 - i4 + 2.0d0 * i5) / omg
        return
    end function intji_dy_re_drv_2

    double precision function intif_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3

        i1 = inti1_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        i2 = inti2_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        intif_ia = i1 + (k + 1) * i2 / (2 * k + 3)
        if (k /= 0) then
            i3 = inti3_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            intif_ia = intif_ia + k * i3 / (2 * k - 1)
        end if
        return
    end function intif_ia

    double precision function inti1_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti1_ia = inti1_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti1_ia = inti1_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti1_ia

    double precision function inti1_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesk_sph_scl(k, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesi1_sph_scl(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x2 - x1)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x0 - x1)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesi1_sph_scl(k, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesk_sph_scl(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x1 - x2)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x1 - x0)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti1_dy_ia = omg / pihalf * (2 * k + 1) * (i1 + i2)
        return
    end function inti1_dy_ia

    double precision function inti2_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti2_ia = inti2_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti2_ia = inti2_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti2_ia

    double precision function inti2_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesk_sph_scl(k+1, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesi1_sph_scl(k+1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x2 - x1)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x0 - x1)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesi1_sph_scl(k+1, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesk_sph_scl(k+1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x1 - x2)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x1 - x0)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti2_dy_ia = omg / pihalf * (2 * k + 3) * (i1 + i2)
        return
    end function inti2_dy_ia

    double precision function inti3_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            inti3_ia = inti3_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            inti3_ia = inti3_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function inti3_ia

    double precision function inti3_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesk_sph_scl(k-1, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesi1_sph_scl(k-1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x2 - x1)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x0 - x1)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesi1_sph_scl(k-1, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesk_sph_scl(k-1, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x1 - x2)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x1 - x0)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        inti3_dy_ia = omg / pihalf * (2 * k - 1) * (i1 + i2)
        return
    end function inti3_dy_ia

    double precision function intji_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        if (omg == 0.0d0) then
            intji_ia = intji_st(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, nodes, weights)
        else
            intji_ia = intji_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        end if
        return
    end function intji_ia

    double precision function intji_dy_ia(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, x0, h0, a1, b1, d1, s1, x1, p1, q1, f1, a2, b2, d2, s2, x2, p2, q2, f2, h1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        i1 = 0.0d0
        x0 = 0.0d0
        h0 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesk_sph_scl(k, omg * x1)
                f1 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesi1_sph_scl(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x2 - x1)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x0 - x1)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i1 = i1 + s1
        end do

        i2 = 0.0d0
        x0 = knt(ord + nintervals)
        h0 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s1 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesi1_sph_scl(k, omg * x1)
                f1 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s2 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesk_sph_scl(k, omg * x2)
                    f2 = p2 * q2
                    s2 = s2 + f2 * exp(omg * (x1 - x2)) * weights(n)
                end do
                s2 = d2 * s2
                h1 = h0 * exp(omg * (x1 - x0)) + s2
                s1 = s1 + f1 * h1 * weights(j)
                x0 = x1
                h0 = h1
            end do
            s1 = d1 * s1
            i2 = i2 + s1
        end do

        intji_dy_ia = omg / pihalf * (2 * k + 1) * (i1 + i2)
        return
    end function intji_dy_ia

    double precision function intif_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: i1, i2, i3

        i1 = inti1_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        i2 = inti2_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
        intif_ia_drv = i1 + (k + 1) * i2 / (2 * k + 3)
        if (k /= 0) then
            i3 = inti3_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes, weights)
            intif_ia_drv = intif_ia_drv + k * i3 / (2 * k - 1)
        end if
        return
    end function intif_ia_drv

    double precision function inti1_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti1_ia_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti1_ia_drv = inti1_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti1_ia_drv

    double precision function inti1_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, scl, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesk_sph_scl(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesk_sph_scl(k-1, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesi1_sph_scl(k+1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesi1_sph_scl(k, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x2 - x1))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x0 - x1))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = wfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesi1_sph_scl(k+1, omg * x1)
                f11 = p1 * q1
                q1 = dbesi1_sph_scl(k, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = wfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesk_sph_scl(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesk_sph_scl(k-1, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x1 - x2))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x1 - x0))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti1_dy_ia_drv = (2 * k + 1) * (i1 - i2) / pihalf
        return
    end function inti1_dy_ia_drv

    double precision function inti2_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti2_ia_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti2_ia_drv = inti2_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti2_ia_drv

    double precision function inti2_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, scl, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesk_sph_scl(k+1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesk_sph_scl(k, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesi1_sph_scl(k+2, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesi1_sph_scl(k+1, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x2 - x1))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x0 - x1))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = qfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesi1_sph_scl(k+2, omg * x1)
                f11 = p1 * q1
                q1 = dbesi1_sph_scl(k+1, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = qfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesk_sph_scl(k+1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesk_sph_scl(k, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x1 - x2))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x1 - x0))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti2_dy_ia_drv = (2 * k + 3) * (i1 - i2) / pihalf
        return
    end function inti2_dy_ia_drv

    double precision function inti3_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        inti3_ia_drv = 0.0d0
        if (omg /= 0.0d0) then
            inti3_ia_drv = inti3_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function inti3_ia_drv

    double precision function inti3_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, scl, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = dbesk_sph_scl(k-1, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesk_sph_scl(k-2, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = omg * x2 * dbesi1_sph_scl(k, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesi1_sph_scl(k-1, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x2 - x1))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x0 - x1))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = pfunc(knt, ord, kapa, cfa, kapc, cfc, k, x1, iv)
                q1 = omg * x1 * dbesi1_sph_scl(k, omg * x1)
                f11 = p1 * q1
                q1 = dbesi1_sph_scl(k-1, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = pfunc(knt, ord, kapb, cfb, kapd, cfd, k, x2, iv)
                    q2 = dbesk_sph_scl(k-1, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesk_sph_scl(k-2, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x1 - x2))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x1 - x0))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        inti3_dy_ia_drv = (2 * k - 1) * (i1 - i2) / pihalf
        return
    end function inti3_dy_ia_drv

    double precision function intji_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        intji_ia_drv = 0.0d0
        if (omg /= 0.0d0) then
            intji_ia_drv = intji_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
            & weights)
        end if
        return
    end function intji_ia_drv

    double precision function intji_dy_ia_drv(knt, ord, kapa, cfa, kapb, cfb, kapc, cfc, kapd, cfd, k, omg, nodes,&
    & weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        integer, intent(in) :: kapc
        double precision, dimension(:), intent(in) :: cfc
        integer, intent(in) :: kapd
        double precision, dimension(:), intent(in) :: cfd
        integer, intent(in) :: k
        double precision, intent(in) :: omg
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x0, h01, h02, i11, i21, a1, b1, d1, s11, s12, x1, p1, q1, f11, f12, a2, b2, d2, s21,&
        & s22, x2, p2, q2, f21, f22, scl, wgt, h11, h12, i12, i22, i1, i2
        integer :: ngl, npoints, nintervals, il, iv, i, j, n

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1

        x0 = 0.0d0
        h01 = 0.0d0
        h02 = 0.0d0
        i11 = 0.0d0
        i21 = 0.0d0
        do i = 1, nintervals
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = 1, ngl
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = dbesk_sph_scl(k, omg * x1)
                f11 = p1 * q1
                q1 = omg * x1 * dbesk_sph_scl(k-1, omg * x1)
                f12 = p1 * q1
                a2 = x0
                b2 = x1
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = omg * x2 * dbesi1_sph_scl(k+1, omg * x2)
                    f21 = p2 * q2
                    q2 = dbesi1_sph_scl(k, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x2 - x1))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x0 - x1))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i11 = i11 + s11
            i21 = i21 + s12
        end do

        x0 = knt(ord + nintervals)
        h01 = 0.0d0
        h02 = 0.0d0
        i12 = 0.0d0
        i22 = 0.0d0
        do i = nintervals, 1, -1
            il = ord - 1 + i
            iv = il
            a1 = knt(il)
            b1 = knt(il+1)
            d1 = 0.5d0 * (b1 - a1)
            s11 = 0.0d0
            s12 = 0.0d0
            do j = ngl, 1, -1
                x1 = d1 * (nodes(j) + 1.0d0) + a1
                p1 = vfunc(knt, ord, kapa, cfa, kapc, cfc, x1, iv)
                q1 = omg * x1 * dbesi1_sph_scl(k+1, omg * x1)
                f11 = p1 * q1
                q1 = dbesi1_sph_scl(k, omg * x1)
                f12 = p1 * q1
                a2 = x1
                b2 = x0
                d2 = 0.5d0 * (b2 - a2)
                s21 = 0.0d0
                s22 = 0.0d0
                do n = 1, ngl
                    x2 = d2 * (nodes(n) + 1.0d0) + a2
                    p2 = vfunc(knt, ord, kapb, cfb, kapd, cfd, x2, iv)
                    q2 = dbesk_sph_scl(k, omg * x2)
                    f21 = p2 * q2
                    q2 = omg * x2 * dbesk_sph_scl(k-1, omg * x2)
                    f22 = p2 * q2
                    scl = exp(omg * (x1 - x2))
                    wgt = weights(n)
                    s21 = s21 + f21 * scl * wgt
                    s22 = s22 + f22 * scl * wgt
                end do
                s21 = d2 * s21
                s22 = d2 * s22
                scl = exp(omg * (x1 - x0))
                h11 = h01 * scl + s21
                h12 = h02 * scl + s22
                wgt = weights(j)
                s11 = s11 + f11 * h11 * wgt
                s12 = s12 + f12 * h12 * wgt
                x0 = x1
                h01 = h11
                h02 = h12
            end do
            s11 = d1 * s11
            s12 = d1 * s12
            i12 = i12 + s11
            i22 = i22 + s12
        end do

        i1 = i11 + i12
        i2 = i21 + i22

        intji_dy_ia_drv = (2 * k + 1) * (i1 - i2) / pihalf
        return
    end function intji_dy_ia_drv

    double precision function dbesj_sph(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x

        double precision, dimension(1) :: ys
        double precision :: xi, alpha
        integer :: sgn, nz

        if (x < 0.0d0) then
            xi = -x
            sgn = (-1)**n
        else
            xi = x
            sgn = 1
        end if
        alpha = 0.5d0 + n
        call dbesj(xi, alpha, 1, ys, nz)
        dbesj_sph = sqrt(pihalf / xi) * ys(1) * sgn
        return
    end function dbesj_sph

    double precision function dbesy_sph(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x

        double precision, dimension(1) :: ys
        double precision :: xi, fnu
        integer :: sgn, nz

        if (x < 0.0d0) then
            xi = -x
            sgn = (-1)**(n+1)
        else
            xi = x
            sgn = 1
        end if
        if (n >= 0) then
            fnu = 0.5d0 + n
            call dbesy(xi, fnu, 1, ys)
        else
            fnu = 0.5d0 - n - 1
            call dbesj(xi, fnu, 1, ys, nz)
            ys(1) = ys(1) * (-1)**(n+1)
        end if
        dbesy_sph = sqrt(pihalf / xi) * ys(1) * sgn
        return
    end function dbesy_sph

    double precision function dbesi1_sph_scl(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x

        double precision, dimension(1) :: ys
        double precision :: alpha
        integer :: nz

        alpha = 0.5d0 + n
        call dbesi(x, alpha, 2, 1, ys, nz)
        dbesi1_sph_scl = sqrt(pihalf / x) * ys(1)
        return
    end function dbesi1_sph_scl

    double precision function dbesk_sph_scl(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x

        double precision, dimension(1) :: ys
        double precision :: fnu
        integer :: nz

        fnu = n + 0.5d0
        if (n < 0) fnu = -fnu
        call dbesk(x, fnu, 2, 1, ys, nz)
        dbesk_sph_scl = sqrt(pihalf / x) * ys(1)
        return
    end function dbesk_sph_scl

end module bspqed_integrals
