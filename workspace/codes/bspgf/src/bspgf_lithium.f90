module bspgf_lithium
    use bspqed_dkb, only: solve, vfunc, wfunc
    use bspqed_function, only: function_type
    use bspqed_integrals, only: intif_re, intji_re, intif_re_drv, intji_re_drv, intic_re, intic_re_drv, intif_ia,&
    & intji_ia, intif_ia_drv, intji_ia_drv, intif_re_drv_2, intic_re_drv_2, intji_re_drv_2
    implicit none
    private

    public :: fbza
    public :: fbza_el
    public :: fbza_ep
    public :: cbza
    public :: cbza_el

contains

    subroutine fbza(knt, ord, nuclpot, magnpot, nodes, weights, irr, red)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: magnpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: irr
        double precision, intent(out) :: red

        double precision, dimension(2*(size(knt)-ord), 2*(size(knt)-ord-2)) :: cfmat
        double precision, dimension((size(knt)-2*(ord-1)-1) * size(nodes)) :: magnvals
        double precision, dimension(2*(size(knt)-ord-2)) :: spc
        double precision, dimension(2*(size(knt)-ord)) :: cf_a, cf_b, cf_vm1_a, cf_vm1_b, cf_vp2_b
        double precision :: epsa, epsb, deab, epsn, v_na, v_nb, i_abab, i_baab, j_baab, j_abab, v_aa, v_bb, tmp
        integer :: ngl, nknots, nbsplines, nspectrum, kap1, kapa, na, kapb, nb, nn, kap2

        ngl = size(nodes)
        nknots = size(knt)
        nbsplines = nknots - ord
        nspectrum = 2 * (nbsplines - 2)

        call storevals(knt, ord, nodes, magnpot, magnvals)

        kap1 = -1
        call solve(knt, ord, kap1, nuclpot, nodes, weights, spc, cfmat)

        kapa = -1
        na = nspectrum / 2 + 2
        epsa = spc(na)
        cf_a = cfmat(:, na)

        kapb = -1
        nb = nspectrum / 2 + 1
        epsb = spc(nb)
        cf_b = cfmat(:, nb)
        deab = epsa - epsb

        cf_vm1_a = 0.0d0
        cf_vm1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsa - epsn) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, magnvals, nodes, weights)
                cf_vm1_a = cf_vm1_a + cfmat(:, nn) * v_na / (epsa - epsn)
            end if
            if (abs(epsb - epsn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vm1_b = cf_vm1_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        kap2 = 2
        call solve(knt, ord, kap2, nuclpot, nodes, weights, spc, cfmat)

        cf_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsb - epsn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vp2_b = cf_vp2_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        irr = 0.0d0

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        irr = irr + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        irr = irr + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        irr = irr + tmp

        red = 0.0d0

        v_aa = intlva(knt, ord, kapa, cf_a, kapa, cf_a, magnvals, nodes, weights)
        v_bb = intlva(knt, ord, kapb, cf_b, kapb, cf_b, magnvals, nodes, weights)

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -2.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        red = red + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 2.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        red = red + tmp

        return
    end subroutine fbza

    subroutine fbza_el(knt, ord, nuclpot, magnpot, elecpot, nodes, weights, a, b, e, g1, g2, g3)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: magnpot
        class(function_type), intent(in) :: elecpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: a
        double precision, intent(out) :: b
        double precision, intent(out) :: e
        double precision, intent(out) :: g1
        double precision, intent(out) :: g2
        double precision, intent(out) :: g3

        double precision, dimension(2*(size(knt)-ord), 2*(size(knt)-ord-2)) :: cfmat
        double precision, dimension((size(knt)-2*(ord-1)-1) * size(nodes)) :: magnvals, elecvals
        double precision, dimension(2*(size(knt)-ord-2)) :: spc
        double precision, dimension(2*(size(knt)-ord)) :: cf_a, cf_b, cf_vm1_a, cf_um1_a, cf2_vm1_a, cf2_um1_a,&
        & cf_vm1_b, cf_um1_b, cf2_vm1_b, cf2_um1_b, cf_vm1_um1_a, cf_um1_vm1_a, cf_vm1_um1_b, cf_um1_vm1_b,&
        & cf_vp2_b, cf2_vp2_b, cf_vp2_um1_b, cf_up2_vp2_b
        double precision :: epsa, epsb, deab, epsn, de_an, v_na, u_na, de_bn, v_nb, u_nb, u_aa, u_bb, v_aa, v_bb,&
        & vu2_aa, vu2_bb, vu_aa, vu_bb, i_abab, i_baab, j_baab, j_abab, tmp
        integer :: ngl, nknots, nbsplines, nspectrum, kap1, kapa, na, kapb, nb, nn, kap2

        a = 0.0d0
        b = 0.0d0
        e = 0.0d0
        g1 = 0.0d0
        g2 = 0.0d0
        g3 = 0.0d0

        ngl = size(nodes)
        nknots = size(knt)
        nbsplines = nknots - ord
        nspectrum = 2 * (nbsplines - 2)

        call storevals(knt, ord, nodes, magnpot, magnvals)
        call storevals(knt, ord, nodes, elecpot, elecvals)

        kap1 = -1
        call solve(knt, ord, kap1, nuclpot, nodes, weights, spc, cfmat)

        kapa = -1
        na = nspectrum / 2 + 2
        epsa = spc(na)
        cf_a = cfmat(:, na)

        kapb = -1
        nb = nspectrum / 2 + 1
        epsb = spc(nb)
        cf_b = cfmat(:, nb)
        deab = epsa - epsb

        cf_vm1_a = 0.0d0
        cf_um1_a = 0.0d0
        cf2_vm1_a = 0.0d0
        cf2_um1_a = 0.0d0
        cf_vm1_b = 0.0d0
        cf_um1_b = 0.0d0
        cf2_vm1_b = 0.0d0
        cf2_um1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_an = epsa - epsn
            if (abs(de_an) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, magnvals, nodes, weights)
                cf_vm1_a = cf_vm1_a + cfmat(:, nn) * v_na / de_an
                cf2_vm1_a = cf2_vm1_a + cfmat(:, nn) * v_na / de_an**2
                u_na = intlua(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, elecvals, nodes, weights)
                cf_um1_a = cf_um1_a + cfmat(:, nn) * u_na / de_an
                cf2_um1_a = cf2_um1_a + cfmat(:, nn) * u_na / de_an**2
            end if
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vm1_b = cf_vm1_b + cfmat(:, nn) * v_nb / de_bn
                cf2_vm1_b = cf2_vm1_b + cfmat(:, nn) * v_nb / de_bn**2
                u_nb = intlua(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, elecvals, nodes, weights)
                cf_um1_b = cf_um1_b + cfmat(:, nn) * u_nb / de_bn
                cf2_um1_b = cf2_um1_b + cfmat(:, nn) * u_nb / de_bn**2
            end if
        end do

        cf_vm1_um1_a = 0.0d0
        cf_um1_vm1_a = 0.0d0
        cf_vm1_um1_b = 0.0d0
        cf_um1_vm1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_an = epsa - epsn
            if (abs(de_an) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kap1, cf_um1_a, magnvals, nodes, weights)
                cf_vm1_um1_a = cf_vm1_um1_a + cfmat(:, nn) * v_na / de_an
                u_na = intlua(knt, ord, kap1, cfmat(:, nn), kap1, cf_vm1_a, elecvals, nodes, weights)
                cf_um1_vm1_a = cf_um1_vm1_a + cfmat(:, nn) * u_na / de_an
            end if
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kap1, cf_um1_b, magnvals, nodes, weights)
                cf_vm1_um1_b = cf_vm1_um1_b + cfmat(:, nn) * v_nb / de_bn
                u_nb = intlua(knt, ord, kap1, cfmat(:, nn), kap1, cf_vm1_b, elecvals, nodes, weights)
                cf_um1_vm1_b = cf_um1_vm1_b + cfmat(:, nn) * u_nb / de_bn
            end if
        end do

        kap2 = 2
        call solve(knt, ord, kap2, nuclpot, nodes, weights, spc, cfmat)

        cf_vp2_b = 0.0d0
        cf2_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vp2_b = cf_vp2_b + cfmat(:, nn) * v_nb / de_bn
                cf2_vp2_b = cf2_vp2_b + cfmat(:, nn) * v_nb / de_bn**2
            end if
        end do

        cf_vp2_um1_b = 0.0d0
        cf_up2_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kap1, cf_um1_b, magnvals, nodes, weights)
                cf_vp2_um1_b = cf_vp2_um1_b + cfmat(:, nn) * v_nb / de_bn
                u_nb = intlua(knt, ord, kap2, cfmat(:, nn), kap2, cf_vp2_b, elecvals, nodes, weights)
                cf_up2_vp2_b = cf_up2_vp2_b + cfmat(:, nn) * u_nb / de_bn
            end if
        end do

        u_aa = intlua(knt, ord, kapa, cf_a, kapa, cf_a, elecvals, nodes, weights)
        u_bb = intlua(knt, ord, kapb, cf_b, kapb, cf_b, elecvals, nodes, weights)

        v_aa = intlva(knt, ord, kapa, cf_a, kapa, cf_a, magnvals, nodes, weights)
        v_bb = intlva(knt, ord, kapb, cf_b, kapb, cf_b, magnvals, nodes, weights)

        vu2_aa = intlva(knt, ord, kapa, cf_a, kap1, cf2_um1_a, magnvals, nodes, weights)
        vu2_bb = intlva(knt, ord, kapb, cf_b, kap1, cf2_um1_b, magnvals, nodes, weights)

        vu_aa = intlva(knt, ord, kapa, cf_a, kap1, cf_um1_a, magnvals, nodes, weights)
        vu_bb = intlva(knt, ord, kapb, cf_b, kap1, cf_um1_b, magnvals, nodes, weights)

        !! --------------------------------------------------------------------------------------------------------
        !! term 01
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * 2.0d0 * i_abab / 3.0d0
        a = a + tmp
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        e = e + tmp
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        e = e + tmp
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0
        a = a + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 02
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * 2.0d0 * i_abab / 3.0d0
        e = e + tmp
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        a = a + tmp
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        a = a + tmp
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0
        e = e + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 03
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_a, kap1, cf_um1_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kap1, cf_um1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kap1, cf_um1_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        e = e + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 04
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_um1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_um1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_um1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_um1_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_um1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        a = a + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 05
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        b = b + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        b = b + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_up2_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_up2_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        b = b + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 06
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf2_vm1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * u_aa
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_b, 1, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_vm1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_vm1_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * u_bb
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf2_vp2_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf2_vp2_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0 * u_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 07
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf2_um1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * v_aa
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_b, 1, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_um1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_um1_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * v_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 08
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 09
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        tmp = -4.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_b, 0, deab, nodes,&
        & weights)
        tmp = -4.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 10
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intif_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intif_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kapb, cf_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * vu2_aa
        g1 = g1 + tmp
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * vu2_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 11
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intif_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -4.0d0 * i_baab * (vu_aa - vu_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (vu_aa + vu_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 12
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intif_re_drv_2(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -2.0d0 * i_baab * (v_aa - v_bb) / 3.0d0 * (u_aa - u_bb)
        g3 = g3 + tmp

        j_baab = intji_re_drv_2(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 2.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0 * (u_aa - u_bb)
        g3 = g3 + tmp

        return
    end subroutine fbza_el

    subroutine fbza_ep(knt, ord, nuclpot, magnpot, nodes, weights, irr, red)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: magnpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: irr
        double precision, intent(out) :: red

        double precision, parameter :: epsabs = 0.0d0
        double precision, parameter :: epsrel = 1.0d-04

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision, dimension(2*(size(knt)-ord), 2*(size(knt)-ord-2)) :: cfmat
        double precision, dimension((size(knt)-2*(ord-1)-1) * size(nodes)) :: magnvals
        double precision, dimension(2*(size(knt)-ord-2)) :: spc
        double precision, dimension(2*(size(knt)-ord)) :: cf_a, cf_b, cf_vm1_a, cf_vm1_b, cf_vp2_b
        double precision :: epsa, epsb, deab, epsn, v_na, v_nb, v_aa, v_bb, abserr
        integer :: ngl, nknots, nbsplines, nspectrum, kap1, kapa, na, kapb, nb, nn, kap2, neval, ier, last

        ngl = size(nodes)
        nknots = size(knt)
        nbsplines = nknots - ord
        nspectrum = 2 * (nbsplines - 2)

        call storevals(knt, ord, nodes, magnpot, magnvals)

        kap1 = -1
        call solve(knt, ord, kap1, nuclpot, nodes, weights, spc, cfmat)

        kapa = -1
        na = nspectrum / 2 + 2
        epsa = spc(na)
        cf_a = cfmat(:, na)

        kapb = -1
        nb = nspectrum / 2 + 1
        epsb = spc(nb)
        cf_b = cfmat(:, nb)
        deab = epsa - epsb

        cf_vm1_a = 0.0d0
        cf_vm1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsa - epsn) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, magnvals, nodes, weights)
                cf_vm1_a = cf_vm1_a + cfmat(:, nn) * v_na / (epsa - epsn)
            end if
            if (abs(epsb - epsn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vm1_b = cf_vm1_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        kap2 = 2
        call solve(knt, ord, kap2, nuclpot, nodes, weights, spc, cfmat)

        cf_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsb - epsn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vp2_b = cf_vp2_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        v_aa = intlva(knt, ord, kapa, cf_a, kapa, cf_a, magnvals, nodes, weights)
        v_bb = intlva(knt, ord, kapb, cf_b, kapb, cf_b, magnvals, nodes, weights)

        call dqagi(firr, 1.0d0, 1, epsabs, epsrel, irr, abserr, neval, ier, limit, lenw, last, iwork, work)
        call dqagi(fred, 1.0d0, 1, epsabs, epsrel, red, abserr, neval, ier, limit, lenw, last, iwork, work)
        return

    contains

        double precision function firr(t)
            double precision, intent(in) :: t

            double precision :: omg0, omg1, i_abab, i_baab, j_baab, j_abab, tmp, f

            omg0 = 2.0d0 * t
            omg1 = sqrt(omg0**2 - deab**2)

            firr = 0.0d0

            i_abab = intif_ia(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_a, kapb, cf_b, 0, omg0, nodes,&
            & weights)
            i_baab = intif_ia(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, omg1, nodes,&
            & weights)
            j_baab = intji_ia(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, omg1, nodes,&
            & weights)
            tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
            firr = firr + tmp

            j_abab = intji_ia(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_b, 1, omg0, nodes,&
            & weights)
            i_baab = intif_ia(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, omg1, nodes,&
            & weights)
            j_baab = intji_ia(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, omg1, nodes,&
            & weights)
            tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
            firr = firr + tmp

            j_abab = intji_ia(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_b, 1, omg0, nodes,&
            & weights)
            j_baab = intji_ia(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, omg1, nodes,&
            & weights)
            tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
            firr = firr + tmp

            f = 2.0d0 * (1.0d0 + 0.5d0 / t**2) * sqrt(t**2 - 1.0d0) / t**2 / 3.0d0
            firr = f * firr

            return
        end function firr

        double precision function fred(t)
            double precision, intent(in) :: t

            double precision :: omg0, omg1, i_baab, j_baab, tmp, f

            omg0 = 2.0d0 * t
            omg1 = sqrt(omg0**2 - deab**2)

            fred = 0.0d0

            i_baab = -intif_ia_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, omg1, nodes,&
            & weights) / omg1 * deab
            tmp = -2.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
            fred = fred + tmp

            j_baab = -intji_ia_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, omg1, nodes,&
            & weights) / omg1 * deab
            tmp = 2.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 &
            & / 3.0d0
            fred = fred + tmp

            f = 2.0d0 * (1.0d0 + 0.5d0 / t**2) * sqrt(t**2 - 1.0d0) / t**2 / 3.0d0
            fred = f * fred

            return
        end function fred
    end subroutine fbza_ep

    subroutine cbza(knt, ord, nuclpot, magnpot, nodes, weights, irr, red)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: magnpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: irr
        double precision, intent(out) :: red

        double precision, dimension(2*(size(knt)-ord), 2*(size(knt)-ord-2)) :: cfmat
        double precision, dimension((size(knt)-2*(ord-1)-1) * size(nodes)) :: magnvals
        double precision, dimension(2*(size(knt)-ord-2)) :: spc
        double precision, dimension(2*(size(knt)-ord)) :: cf_a, cf_b, cf_vm1_a, cf_vm1_b, cf_vp2_b
        double precision :: epsa, epsb, deab, epsn, v_na, v_nb, i_abab, i_baab, j_baab, j_abab, v_aa, v_bb, tmp
        integer :: ngl, nknots, nbsplines, nspectrum, kap1, kapa, na, kapb, nb, nn, kap2

        ngl = size(nodes)
        nknots = size(knt)
        nbsplines = nknots - ord
        nspectrum = 2 * (nbsplines - 2)

        call storevals(knt, ord, nodes, magnpot, magnvals)

        kap1 = -1
        call solve(knt, ord, kap1, nuclpot, nodes, weights, spc, cfmat)

        kapa = -1
        na = nspectrum / 2 + 2
        epsa = spc(na)
        cf_a = cfmat(:, na)

        kapb = -1
        nb = nspectrum / 2 + 1
        epsb = spc(nb)
        cf_b = cfmat(:, nb)
        deab = epsa - epsb

        cf_vm1_a = 0.0d0
        cf_vm1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsa - epsn) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, magnvals, nodes, weights)
                cf_vm1_a = cf_vm1_a + cfmat(:, nn) * v_na / (epsa - epsn)
            end if
            if (abs(epsb - epsn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vm1_b = cf_vm1_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        kap2 = 2
        call solve(knt, ord, kap2, nuclpot, nodes, weights, spc, cfmat)

        cf_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            if (abs(epsb - epsn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vp2_b = cf_vp2_b + cfmat(:, nn) * v_nb / (epsb - epsn)
            end if
        end do

        irr = 0.0d0

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        irr = irr + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        irr = irr + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        irr = irr + tmp

        red = 0.0d0

        v_aa = intlva(knt, ord, kapa, cf_a, kapa, cf_a, magnvals, nodes, weights)
        v_bb = intlva(knt, ord, kapb, cf_b, kapb, cf_b, magnvals, nodes, weights)

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -2.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        red = red + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 2.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        red = red + tmp

        return
    end subroutine cbza

    subroutine cbza_el(knt, ord, nuclpot, magnpot, elecpot, nodes, weights, a, b, e, g1, g2, g3)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        class(function_type), intent(in) :: nuclpot
        class(function_type), intent(in) :: magnpot
        class(function_type), intent(in) :: elecpot
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision, intent(out) :: a
        double precision, intent(out) :: b
        double precision, intent(out) :: e
        double precision, intent(out) :: g1
        double precision, intent(out) :: g2
        double precision, intent(out) :: g3

        double precision, dimension(2*(size(knt)-ord), 2*(size(knt)-ord-2)) :: cfmat
        double precision, dimension((size(knt)-2*(ord-1)-1) * size(nodes)) :: magnvals, elecvals
        double precision, dimension(2*(size(knt)-ord-2)) :: spc
        double precision, dimension(2*(size(knt)-ord)) :: cf_a, cf_b, cf_vm1_a, cf_um1_a, cf2_vm1_a, cf2_um1_a,&
        & cf_vm1_b, cf_um1_b, cf2_vm1_b, cf2_um1_b, cf_vm1_um1_a, cf_um1_vm1_a, cf_vm1_um1_b, cf_um1_vm1_b,&
        & cf_vp2_b, cf2_vp2_b, cf_vp2_um1_b, cf_up2_vp2_b
        double precision :: epsa, epsb, deab, epsn, de_an, v_na, u_na, de_bn, v_nb, u_nb, u_aa, u_bb, v_aa, v_bb,&
        & vu2_aa, vu2_bb, vu_aa, vu_bb, i_abab, i_baab, j_baab, j_abab, tmp
        integer :: ngl, nknots, nbsplines, nspectrum, kap1, kapa, na, kapb, nb, nn, kap2

        a = 0.0d0
        b = 0.0d0
        e = 0.0d0
        g1 = 0.0d0
        g2 = 0.0d0
        g3 = 0.0d0

        ngl = size(nodes)
        nknots = size(knt)
        nbsplines = nknots - ord
        nspectrum = 2 * (nbsplines - 2)

        call storevals(knt, ord, nodes, magnpot, magnvals)
        call storevals(knt, ord, nodes, elecpot, elecvals)

        kap1 = -1
        call solve(knt, ord, kap1, nuclpot, nodes, weights, spc, cfmat)

        kapa = -1
        na = nspectrum / 2 + 2
        epsa = spc(na)
        cf_a = cfmat(:, na)

        kapb = -1
        nb = nspectrum / 2 + 1
        epsb = spc(nb)
        cf_b = cfmat(:, nb)
        deab = epsa - epsb

        cf_vm1_a = 0.0d0
        cf_um1_a = 0.0d0
        cf2_vm1_a = 0.0d0
        cf2_um1_a = 0.0d0
        cf_vm1_b = 0.0d0
        cf_um1_b = 0.0d0
        cf2_vm1_b = 0.0d0
        cf2_um1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_an = epsa - epsn
            if (abs(de_an) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, magnvals, nodes, weights)
                cf_vm1_a = cf_vm1_a + cfmat(:, nn) * v_na / de_an
                cf2_vm1_a = cf2_vm1_a + cfmat(:, nn) * v_na / de_an**2
                u_na = intlua(knt, ord, kap1, cfmat(:, nn), kapa, cf_a, elecvals, nodes, weights)
                cf_um1_a = cf_um1_a + cfmat(:, nn) * u_na / de_an
                cf2_um1_a = cf2_um1_a + cfmat(:, nn) * u_na / de_an**2
            end if
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vm1_b = cf_vm1_b + cfmat(:, nn) * v_nb / de_bn
                cf2_vm1_b = cf2_vm1_b + cfmat(:, nn) * v_nb / de_bn**2
                u_nb = intlua(knt, ord, kap1, cfmat(:, nn), kapb, cf_b, elecvals, nodes, weights)
                cf_um1_b = cf_um1_b + cfmat(:, nn) * u_nb / de_bn
                cf2_um1_b = cf2_um1_b + cfmat(:, nn) * u_nb / de_bn**2
            end if
        end do

        cf_vm1_um1_a = 0.0d0
        cf_um1_vm1_a = 0.0d0
        cf_vm1_um1_b = 0.0d0
        cf_um1_vm1_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_an = epsa - epsn
            if (abs(de_an) /= 0.0d0 .and. nn /= na) then
                v_na = intlva(knt, ord, kap1, cfmat(:, nn), kap1, cf_um1_a, magnvals, nodes, weights)
                cf_vm1_um1_a = cf_vm1_um1_a + cfmat(:, nn) * v_na / de_an
                u_na = intlua(knt, ord, kap1, cfmat(:, nn), kap1, cf_vm1_a, elecvals, nodes, weights)
                cf_um1_vm1_a = cf_um1_vm1_a + cfmat(:, nn) * u_na / de_an
            end if
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0 .and. nn /= nb) then
                v_nb = intlva(knt, ord, kap1, cfmat(:, nn), kap1, cf_um1_b, magnvals, nodes, weights)
                cf_vm1_um1_b = cf_vm1_um1_b + cfmat(:, nn) * v_nb / de_bn
                u_nb = intlua(knt, ord, kap1, cfmat(:, nn), kap1, cf_vm1_b, elecvals, nodes, weights)
                cf_um1_vm1_b = cf_um1_vm1_b + cfmat(:, nn) * u_nb / de_bn
            end if
        end do

        kap2 = 2
        call solve(knt, ord, kap2, nuclpot, nodes, weights, spc, cfmat)

        cf_vp2_b = 0.0d0
        cf2_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kapb, cf_b, magnvals, nodes, weights)
                cf_vp2_b = cf_vp2_b + cfmat(:, nn) * v_nb / de_bn
                cf2_vp2_b = cf2_vp2_b + cfmat(:, nn) * v_nb / de_bn**2
            end if
        end do

        cf_vp2_um1_b = 0.0d0
        cf_up2_vp2_b = 0.0d0
        do nn = 1, nspectrum
            epsn = spc(nn)
            de_bn = epsb - epsn
            if (abs(de_bn) /= 0.0d0) then
                v_nb = intlva(knt, ord, kap2, cfmat(:, nn), kap1, cf_um1_b, magnvals, nodes, weights)
                cf_vp2_um1_b = cf_vp2_um1_b + cfmat(:, nn) * v_nb / de_bn
                u_nb = intlua(knt, ord, kap2, cfmat(:, nn), kap2, cf_vp2_b, elecvals, nodes, weights)
                cf_up2_vp2_b = cf_up2_vp2_b + cfmat(:, nn) * u_nb / de_bn
            end if
        end do

        u_aa = intlua(knt, ord, kapa, cf_a, kapa, cf_a, elecvals, nodes, weights)
        u_bb = intlua(knt, ord, kapb, cf_b, kapb, cf_b, elecvals, nodes, weights)

        v_aa = intlva(knt, ord, kapa, cf_a, kapa, cf_a, magnvals, nodes, weights)
        v_bb = intlva(knt, ord, kapb, cf_b, kapb, cf_b, magnvals, nodes, weights)

        vu2_aa = intlva(knt, ord, kapa, cf_a, kap1, cf2_um1_a, magnvals, nodes, weights)
        vu2_bb = intlva(knt, ord, kapb, cf_b, kap1, cf2_um1_b, magnvals, nodes, weights)

        vu_aa = intlva(knt, ord, kapa, cf_a, kap1, cf_um1_a, magnvals, nodes, weights)
        vu_bb = intlva(knt, ord, kapb, cf_b, kap1, cf_um1_b, magnvals, nodes, weights)

        !! --------------------------------------------------------------------------------------------------------
        !! term 01
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * 2.0d0 * i_abab / 3.0d0
        a = a + tmp
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        e = e + tmp
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kap1, cf_um1_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kap1, cf_um1_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        e = e + tmp
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0
        a = a + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 02
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kap1, cf_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * 2.0d0 * i_abab / 3.0d0
        e = e + tmp
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kapa, cf_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        a = a + tmp
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kap1, cf_um1_b, kapa, cf_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kap1, cf_um1_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * j_abab / 9.0d0 / 3.0d0
        a = a + tmp
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0
        e = e + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 03
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_a, kap1, cf_um1_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kap1, cf_um1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kap1, cf_um1_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_a, kap1, cf_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        e = e + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_a, kap2, cf_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        e = e + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 04
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_vm1_um1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_um1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_um1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_um1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        a = a + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_vp2_um1_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_um1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        a = a + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 05
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf_um1_vm1_a, kapb, cf_b, 0, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0
        b = b + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf_um1_vm1_b, 1, 0.0d0, nodes,&
        & weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0
        b = b + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf_up2_vp2_b, 1, 0.0d0, nodes,&
        & weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_up2_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0
        b = b + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 06
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf2_vm1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * u_aa
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf2_vm1_b, 1, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_vm1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_vm1_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * u_bb
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap2, cf2_vp2_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf2_vp2_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * 2.0d0 * (j_abab - j_baab) / 9.0d0 / 3.0d0 * u_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 07
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kap1, cf2_um1_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_a, kapb, cf_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * v_aa
        g1 = g1 + tmp

        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kap1, cf2_um1_b, 1, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_um1_b, 0, deab, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf2_um1_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * v_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 08
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_vm1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * (-i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 0, deab, nodes,&
        & weights)
        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_vm1_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * (i_baab + 2.0d0 * j_baab / 9.0d0) / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap2, cf_vp2_b, 1, deab, nodes,&
        & weights)
        tmp = -4.0d0 * 2.0d0 * (-j_baab) / 9.0d0 / 3.0d0 * (u_aa - u_bb)
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 09
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kapb, cf_b, 0, deab, nodes,&
        & weights)
        tmp = -4.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kap1, cf_um1_a, kapb, cf_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_b, 0, deab, nodes,&
        & weights)
        tmp = -4.0d0 * i_baab * (v_aa - v_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kap1, cf_um1_b, 1, deab, nodes,&
        & weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 10
        !! --------------------------------------------------------------------------------------------------------

        i_abab = intic_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kapb, cf_b, 0, 0.0d0, nodes, weights)
        i_baab = intic_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        j_abab = intji_re(knt, ord, kapa, cf_a, kapb, cf_b, kapa, cf_a, kapb, cf_b, 1, 0.0d0, nodes, weights)
        j_baab = intji_re(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = -4.0d0 * (2.0d0 * i_abab - i_baab + 2.0d0 * j_baab / 3.0d0) / 3.0d0 * vu2_aa
        g1 = g1 + tmp
        tmp = 4.0d0 * (i_baab + 2.0d0 * (2.0d0 * j_abab + j_baab) / 9.0d0) / 3.0d0 * vu2_bb
        g1 = g1 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 11
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intic_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -4.0d0 * i_baab * (vu_aa - vu_bb) / 3.0d0
        g2 = g2 + tmp

        j_baab = intji_re_drv(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 4.0d0 * j_baab * 2.0d0 * (vu_aa + vu_bb / 3.0d0) / 3.0d0 / 3.0d0
        g2 = g2 + tmp

        !! --------------------------------------------------------------------------------------------------------
        !! term 12
        !! --------------------------------------------------------------------------------------------------------

        i_baab = intic_re_drv_2(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 0, deab, nodes, weights)
        tmp = -2.0d0 * i_baab * (v_aa - v_bb) / 3.0d0 * (u_aa - u_bb)
        g3 = g3 + tmp

        j_baab = intji_re_drv_2(knt, ord, kapb, cf_b, kapa, cf_a, kapa, cf_a, kapb, cf_b, 1, deab, nodes, weights)
        tmp = 2.0d0 * j_baab * 2.0d0 * (v_aa + v_bb / 3.0d0) / 3.0d0 / 3.0d0 * (u_aa - u_bb)
        g3 = g3 + tmp

        return
    end subroutine cbza_el

    double precision function intlva(knt, ord, kapa, cfa, kapb, cfb, potvals, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        double precision, dimension(:), intent(in) :: potvals
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x, s, t, a, b, d, vab, pot
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
                pot = potvals((i-1) * ngl + j)
                t = t + pot * vab * weights(j)
            enddo
            t = d * t
            s = s + t
        enddo
        intlva = s
        return
    end function intlva

    double precision function intlua(knt, ord, kapa, cfa, kapb, cfb, potvals, nodes, weights)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: cfa
        integer, intent(in) :: kapb
        double precision, dimension(:), intent(in) :: cfb
        double precision, dimension(:), intent(in) :: potvals
        double precision, dimension(:), intent(in) :: nodes
        double precision, dimension(:), intent(in) :: weights

        double precision :: x, s, t, a, b, d, wab, pot
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
                wab = wfunc(knt, ord, kapa, cfa, kapb, cfb, x, iv)
                pot = potvals((i-1) * ngl + j)
                t = t + pot * wab * weights(j)
            enddo
            t = d * t
            s = s + t
        enddo
        intlua = s
        return
    end function intlua

    subroutine storevals(knt, ord, nodes, func, vals)
        double precision, dimension(:), intent(in) :: knt
        integer, intent(in) :: ord
        double precision, dimension(:), intent(in) :: nodes
        class(function_type), intent(in) :: func

        double precision, dimension(:), intent(out) :: vals

        double precision :: a, b, d, x
        integer :: ngl, npoints, nintervals, il, i, j

        ngl = size(nodes)
        npoints = size(knt) - 2 * (ord - 1)
        nintervals = npoints - 1
        do i = 1, nintervals
            il = ord - 1 + i
            a = knt(il)
            b = knt(il+1)
            d = 0.5d0 * (b - a)
            do j = 1, ngl
                x = d * (nodes(j) + 1.0d0) + a
                vals((i-1) * ngl + j) = func%eval(x)
            end do
        end do
        return
    end subroutine storevals

end module bspgf_lithium
