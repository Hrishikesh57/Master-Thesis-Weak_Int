module bspse_manypotential
    implicit none
    private

    public :: manypot_low
    public :: manypot_high

    double complex, parameter :: ii = dcmplx(0.0d0, 1.0d0)
    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
    double precision, parameter :: pihalf = 0.5d0 * pi

contains

    double precision function manypot_low(na, kapa, speca, ga1, fa1, kapi, specb, specf, gf1, ff1, vmat, eps1s,&
    & nod1, wgt1, epsabs, epsrel)
        integer, intent(in) :: na
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: speca
        double precision, dimension(:), intent(in) :: ga1
        double precision, dimension(:), intent(in) :: fa1
        integer, intent(in) :: kapi
        double precision, dimension(:), intent(in) :: specb
        double precision, dimension(:), intent(in) :: specf
        double precision, dimension(:, :), intent(in) :: gf1
        double precision, dimension(:, :), intent(in) :: ff1
        double precision, dimension(:, :), intent(in) :: vmat
        double precision, intent(in) :: eps1s
        double precision, dimension(:), intent(in) :: nod1
        double precision, dimension(:), intent(in) :: wgt1
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision :: epsa, hlf0, delta, res, abserr
        integer :: nspca, lmain, nmain, neval, ier, last

        epsa = speca(na)
        hlf0 = 0.5d0 * (epsa - 0.9d0 * eps1s)
        nspca = size(speca)
        lmain = (abs(2 * kapa + 1) - 1) / 2
        nmain = na - nspca / 2 + lmain
        delta = epsa / (5.0d0 * nmain**(1.0d0/3.0d0))

        call dqags(f, 0.0d0, 1.0d0, epsabs, epsrel, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        manypot_low = res
        return

    contains

        double precision function f(t)
            double precision, intent(in) :: t

            double complex, dimension(size(gf1, 1), size(specb)) :: gphi1, fphi1
            double complex :: emit, omega, weight, rint
            double precision :: cost, sint, deno

            cost = cos(pi * t)
            sint = sin(pi * t)
            emit = cost - ii * sint
            deno = sqrt((delta * cost)**2 + sint**2)
            omega = hlf0 * (1.0d0 + delta * emit / deno)
            weight = pi * hlf0 * delta * (sint + ii * delta**2 * cost) / deno**3
            call phiwf(epsa - omega, specf, gf1, ff1, vmat, gphi1, fphi1)
            rint = sumnl1(omega, kapa, epsa, ga1, fa1, kapi, specb, gphi1, fphi1, nod1, wgt1)
            f = dble(ii * rint * weight)
            return
        end function f
    end function manypot_low

    double precision function manypot_high(na, kapa, speca, ga1, fa1, ga2, fa2, kapi, specb, specf, gf1, ff1, gf2,&
    & ff2, vmat, eps1s, nod1, wgt1, nod2, wgt2, epsabs, epsrel)
        integer, intent(in) :: na
        integer, intent(in) :: kapa
        double precision, dimension(:), intent(in) :: speca
        double precision, dimension(:), intent(in) :: ga1
        double precision, dimension(:), intent(in) :: fa1
        double precision, dimension(:), intent(in) :: ga2
        double precision, dimension(:), intent(in) :: fa2
        integer, intent(in) :: kapi
        double precision, dimension(:), intent(in) :: specb
        double precision, dimension(:), intent(in) :: specf
        double precision, dimension(:, :), intent(in) :: gf1
        double precision, dimension(:, :), intent(in) :: ff1
        double precision, dimension(:, :), intent(in) :: gf2
        double precision, dimension(:, :), intent(in) :: ff2
        double precision, dimension(:, :), intent(in) :: vmat
        double precision, intent(in) :: eps1s
        double precision, dimension(:), intent(in) :: nod1
        double precision, dimension(:), intent(in) :: wgt1
        double precision, dimension(:), intent(in) :: nod2
        double precision, dimension(:), intent(in) :: wgt2
        double precision, intent(in) :: epsabs
        double precision, intent(in) :: epsrel

        integer, parameter :: limit = 100
        integer, parameter :: lenw = 4 * limit

        double precision, dimension(lenw) :: work
        integer, dimension(limit) :: iwork

        double precision :: epsa, eps0, res, abserr
        integer :: neval, ier, last

        epsa = speca(na)
        eps0 = epsa - 0.9d0 * eps1s

        call dqagi(f, 0.0d0, 1, epsabs, epsrel, res, abserr, neval, ier, limit, lenw, last, iwork, work)
        manypot_high = res
        return

    contains

        double precision function f(t)
            double precision, intent(in) :: t

            double complex, dimension(size(gf1, 1), size(specb)) :: gphi1, fphi1
            double complex, dimension(size(gf2, 1), size(specb)) :: gphi2, fphi2
            double complex :: omega, weight, rint

            omega = eps0 + ii * t
            weight = ii
            call phiwf(epsa - omega, specf, gf1, ff1, vmat, gphi1, fphi1)
            call phiwf(epsa - omega, specf, gf2, ff2, vmat, gphi2, fphi2)
            rint = sumnl2(omega, kapa, epsa, ga1, fa1, ga2, fa2, kapi, specb, gphi1, fphi1, gphi2, fphi2, nod1,&
            & wgt1, nod2, wgt2)
            f = dble(ii * rint * weight)
            return
        end function f
    end function manypot_high

    subroutine phiwf(omega, specf, gf, ff, vmat, gphi, fphi)
        double complex, intent(in) :: omega
        double precision, dimension(:), intent(in) :: specf
        double precision, dimension(:, :), intent(in) :: gf
        double precision, dimension(:, :), intent(in) :: ff
        double precision, dimension(:, :), intent(in) :: vmat

        double complex, dimension(:, :), intent(out) :: gphi
        double complex, dimension(:, :), intent(out) :: fphi

        double complex :: denom
        double precision :: epsp, vpi
        integer :: nspecb, nspecf, ni, np

        nspecb = size(gphi, 2)
        nspecf = size(specf)
        gphi = 0.0d0
        fphi = 0.0d0
        do ni = 1, nspecb
            do np = 1, nspecf
                epsp = specf(np)
                vpi = vmat(np, ni)
                denom = omega - epsp
                gphi(:, ni) = gphi(:, ni) + vpi * gf(:, np) / denom
                fphi(:, ni) = fphi(:, ni) + vpi * ff(:, np) / denom
            end do
        end do
        return
    end subroutine phiwf

    double complex function sumnl1(omega, kapa, epsa, ga, fa, kapi, speci, gi, fi, nod, wgt)
        double complex, intent(in) :: omega
        integer, intent(in) :: kapa
        double precision, intent(in) :: epsa
        double precision, dimension(:), intent(in) :: ga
        double precision, dimension(:), intent(in) :: fa
        integer, intent(in) :: kapi
        double precision, dimension(:), intent(in) :: speci
        double complex, dimension(:, :), intent(in) :: gi
        double complex, dimension(:, :), intent(in) :: fi
        double precision, dimension(:), intent(in) :: nod
        double precision, dimension(:), intent(in) :: wgt

        double complex, dimension(max(0, abs(abs(kapa)-abs(kapi))-1):(abs(kapa)+abs(kapi)), size(nod)) :: sphj
        double complex, dimension(size(nod)) :: uai, vai, wai, uia, via, wia
        double complex, dimension(max(0, abs(abs(kapa)-abs(kapi))-1):(abs(kapa)+abs(kapi))) :: iuu, iuv, ivu, ivv,&
        & iww
        double precision, dimension(abs(abs(kapa)-abs(kapi)):(abs(kapa)+abs(kapi)-1)) :: c3j
        double complex :: z, sumn, suml, rl
        double precision :: ja, ji, l1min, l1max, epsi
        integer :: ndim, ierr, lmin, lmax, lmin_bes, nnod, nz, dkap, skap, la, li, nspectrum, ni, l, i

        ja = abs(kapa) - 0.5d0
        ji = abs(kapi) - 0.5d0
        ndim = size(c3j)
        call drc3jj(ja, ji, -0.5d0, 0.5d0, l1min, l1max, c3j, ndim, ierr)
        lmin = int(l1min)
        lmax = int(l1max)

        lmin_bes = max(0, lmin-1)
        nnod = size(nod)
        do i = 1, nnod
            z = omega * nod(i)
            call zbesj_sph(z, lmin_bes, 1, sphj(:, i), nz, ierr)
        end do

        dkap = kapa - kapi
        skap = kapa + kapi
        la = (abs(2 * kapa + 1) - 1) / 2
        li = (abs(2 * kapi + 1) - 1) / 2
        nspectrum = size(speci)
        sumn = 0.0d0
        do ni = 1, nspectrum
            epsi = speci(ni)

            uai = -ga * fi(:, ni) + fa * gi(:, ni)
            vai = -ga * fi(:, ni) - fa * gi(:, ni)
            wai =  ga * gi(:, ni) + fa * fi(:, ni)

            uia = -gi(:, ni) * fa + fi(:, ni) * ga
            via = -gi(:, ni) * fa - fi(:, ni) * ga
            wia =  gi(:, ni) * ga + fi(:, ni) * fa

            call intb(sphj, uai, uia, vai, via, wai, wia, nod, wgt, iuu, iuv, ivu, ivv, iww)

            suml = 0.0d0
            do l = lmin, lmax
                if (mod(la + li + l, 2) == 0) then
                    rl = (2 * l + 1) * iww(l)
                    rl = rl + (l + 1) * iuu(l+1) - dkap * (iuv(l+1) - ivu(l+1)) - dkap**2 * ivv(l+1) / (l + 1)
                    if (l /= 0) then
                        rl = rl + l * iuu(l-1) + dkap * (iuv(l-1) - ivu(l-1)) - dkap**2 * ivv(l-1) / l
                    end if
                    rl = ii * omega * rl
                else
                    if (skap /= 0) then
                        rl = -skap**2 * (2 * l + 1) * ivv(l) / (l * (l + 1))
                        rl = ii * omega * rl
                    else
                        rl = 0.0d0
                    end if
                end if
                suml = suml + c3j(l)**2 * rl
            end do
            sumn = sumn + suml / (epsa - epsi - omega)
        end do
        sumnl1 = 2.0d0 * abs(kapi) * sumn
        return
    end function sumnl1

    double complex function sumnl2(omega, kapa, epsa, ga1, fa1, ga2, fa2, kapi, speci, gi1, fi1, gi2, fi2, nod1,&
    & wgt1, nod2, wgt2)
        double complex, intent(in) :: omega
        integer, intent(in) :: kapa
        double precision, intent(in) :: epsa
        double precision, dimension(:), intent(in) :: ga1
        double precision, dimension(:), intent(in) :: fa1
        double precision, dimension(:), intent(in) :: ga2
        double precision, dimension(:), intent(in) :: fa2
        integer, intent(in) :: kapi
        double precision, dimension(:), intent(in) :: speci
        double complex, dimension(:, :), intent(in) :: gi1
        double complex, dimension(:, :), intent(in) :: fi1
        double complex, dimension(:, :), intent(in) :: gi2
        double complex, dimension(:, :), intent(in) :: fi2
        double precision, dimension(:), intent(in) :: nod1
        double precision, dimension(:), intent(in) :: wgt1
        double precision, dimension(:), intent(in) :: nod2
        double precision, dimension(:), intent(in) :: wgt2

        double complex, dimension(max(0, abs(abs(kapa)-abs(kapi))-1):(abs(kapa)+abs(kapi)), size(nod1)) :: sphh1
        double complex, dimension(max(0, abs(abs(kapa)-abs(kapi))-1):(abs(kapa)+abs(kapi)), size(nod2)) :: sphj2
        double complex, dimension(size(nod1)) :: uai, vai, wai
        double complex, dimension(size(nod2)) :: uia, via, wia
        double complex, dimension(max(0, abs(abs(kapa)-abs(kapi))-1):(abs(kapa)+abs(kapi))) :: iuu, iuv, ivu, ivv,&
        & iww
        double precision, dimension(abs(abs(kapa)-abs(kapi)):(abs(kapa)+abs(kapi)-1)) :: c3j
        double complex :: z, sumn, suml, rl
        double precision :: ja, ji, l1min, l1max, epsi
        integer :: ndim, ierr, lmin, lmax, lmin_bes, nnod, nz, dkap, skap, la, li, nspectrum, ni, l, i

        ja = abs(kapa) - 0.5d0
        ji = abs(kapi) - 0.5d0
        ndim = size(c3j)
        call drc3jj(ja, ji, -0.5d0, 0.5d0, l1min, l1max, c3j, ndim, ierr)
        lmin = int(l1min)
        lmax = int(l1max)

        lmin_bes = max(0, lmin-1)
        nnod = size(nod1)
        do i = 1, nnod
            z = omega * nod1(i)
            call zbesh1_sph(z, lmin_bes, 2, sphh1(:, i), nz, ierr)
        end do
        nnod = size(nod2)
        do i = 1, nnod
            z = omega * nod2(i)
            call zbesj_sph(z, lmin_bes, 2, sphj2(:, i), nz, ierr)
        end do

        dkap = kapa - kapi
        skap = kapa + kapi
        la = (abs(2 * kapa + 1) - 1) / 2
        li = (abs(2 * kapi + 1) - 1) / 2
        nspectrum = size(speci)
        sumn = 0.0d0
        do ni = 1, nspectrum
            epsi = speci(ni)

            uai = -ga1 * fi1(:, ni) + fa1 * gi1(:, ni)
            vai = -ga1 * fi1(:, ni) - fa1 * gi1(:, ni)
            wai =  ga1 * gi1(:, ni) + fa1 * fi1(:, ni)

            uia = -gi2(:, ni) * fa2 + fi2(:, ni) * ga2
            via = -gi2(:, ni) * fa2 - fi2(:, ni) * ga2
            wia =  gi2(:, ni) * ga2 + fi2(:, ni) * fa2

            call inta(omega, sphj2, sphh1, uai, uia, vai, via, wai, wia, nod1, wgt1, nod2, wgt2, iuu, iuv, ivu,&
            & ivv, iww)

            suml = 0.0d0
            do l = lmin, lmax
                if (mod(la + li + l, 2) == 0) then
                    rl = (2 * l + 1) * iww(l)
                    rl = rl + (l + 1) * iuu(l+1) - dkap * (iuv(l+1) - ivu(l+1)) - dkap**2 * ivv(l+1) / (l + 1)
                    if (l /= 0) then
                        rl = rl + l * iuu(l-1) + dkap * (iuv(l-1) - ivu(l-1)) - dkap**2 * ivv(l-1) / l
                    end if
                    rl = ii * omega * rl
                else
                    if (skap /= 0) then
                        rl = -skap**2 * (2 * l + 1) * ivv(l) / (l * (l + 1))
                        rl = ii * omega * rl
                    else
                        rl = 0.0d0
                    end if
                end if
                suml = suml + c3j(l)**2 * rl
            end do
            sumn = sumn + suml / (epsa - epsi - omega)
        end do
        sumnl2 = 2.0d0 * abs(kapi) * sumn
        return
    end function sumnl2

    subroutine inta(omega, sphj2, sphh1, uai, uia, vai, via, wai, wia, nod1, wgt1, nod2, wgt2, iuu, iuv, ivu, ivv,&
    & iww)
        double complex, intent(in) :: omega
        double complex, dimension(:, :), intent(in) :: sphj2
        double complex, dimension(:, :), intent(in) :: sphh1
        double complex, dimension(:), intent(in) :: uai
        double complex, dimension(:), intent(in) :: uia
        double complex, dimension(:), intent(in) :: vai
        double complex, dimension(:), intent(in) :: via
        double complex, dimension(:), intent(in) :: wai
        double complex, dimension(:), intent(in) :: wia
        double precision, dimension(:), intent(in) :: nod1
        double precision, dimension(:), intent(in) :: wgt1
        double precision, dimension(:), intent(in) :: nod2
        double precision, dimension(:), intent(in) :: wgt2

        double complex, dimension(:), intent(out) :: iuu
        double complex, dimension(:), intent(out) :: iuv
        double complex, dimension(:), intent(out) :: ivu
        double complex, dimension(:), intent(out) :: ivv
        double complex, dimension(:), intent(out) :: iww

        double complex, dimension(size(iuu)) :: i1uu, i1uv, i1vu, i1vv, i1ww, i2uu, i2uv, i2vu, i2vv, i2ww

        call inta1(omega, sphj2, sphh1, uai, uia, vai, via, wai, wia, nod1, wgt1, nod2, wgt2, i1uu, i1uv, i1vu,&
        & i1vv, i1ww)

        i2uu =  i1uu
        i2uv = -i1vu
        i2vu = -i1uv
        i2vv =  i1vv
        i2ww =  i1ww

        iuu = i1uu + i2uu
        iuv = i1uv + i2uv
        ivu = i1vu + i2vu
        ivv = i1vv + i2vv
        iww = i1ww + i2ww
        return
    end subroutine inta

    subroutine inta1(omega, sphj, sphh, uac, ubd, vac, vbd, wac, wbd, nod1, wgt1, nod2, wgt2, iuu, iuv, ivu, ivv,&
    & iww)
        double complex, intent(in) :: omega
        double complex, dimension(:, :), intent(in) :: sphj
        double complex, dimension(:, :), intent(in) :: sphh
        double complex, dimension(:), intent(in) :: uac
        double complex, dimension(:), intent(in) :: ubd
        double complex, dimension(:), intent(in) :: vac
        double complex, dimension(:), intent(in) :: vbd
        double complex, dimension(:), intent(in) :: wac
        double complex, dimension(:), intent(in) :: wbd
        double precision, dimension(:), intent(in) :: nod1
        double precision, dimension(:), intent(in) :: wgt1
        double precision, dimension(:), intent(in) :: nod2
        double precision, dimension(:), intent(in) :: wgt2

        double complex, dimension(:), intent(out) :: iuu
        double complex, dimension(:), intent(out) :: iuv
        double complex, dimension(:), intent(out) :: ivu
        double complex, dimension(:), intent(out) :: ivv
        double complex, dimension(:), intent(out) :: iww

        double complex, dimension(size(sphj, 1)) :: s1uu, s1uv, s1vu, s1vv, s1ww, h0u, h0v, h0w, s2u, s2v, s2w,&
        & spj, sph, h1u, h1v, h1w
        double complex :: scl
        double precision :: romega, iomega, x0, x1, w1, x2, w2
        integer :: m, n, j2, i, j

        m = size(nod1)
        n = size(nod2) / m
        romega = dble(omega)
        iomega = aimag(omega)
        s1uu = 0.0d0
        s1uv = 0.0d0
        s1vu = 0.0d0
        s1vv = 0.0d0
        s1ww = 0.0d0
        x0 = 0.0d0
        h0u = 0.0d0
        h0v = 0.0d0
        h0w = 0.0d0
        do i = 1, m
            x1 = nod1(i)
            w1 = wgt1(i)
            s2u = 0.0d0
            s2v = 0.0d0
            s2w = 0.0d0
            do j = 1, n
                j2 = (i-1) * n + j
                x2 = nod2(j2)
                w2 = wgt2(j2)
                scl = exp(iomega * (x2 - x1))
                spj = sphj(:, j2)
                s2u = s2u + spj * ubd(j2) * scl * w2
                s2v = s2v + spj * vbd(j2) * scl * w2
                s2w = s2w + spj * wbd(j2) * scl * w2
            end do
            scl = exp(iomega * (x0 - x1))
            h1u = h0u * scl + s2u
            h1v = h0v * scl + s2v
            h1w = h0w * scl + s2w
            sph = sphh(:, i) * exp(ii * romega * x1)
            s1uu = s1uu + sph * uac(i) * h1u * w1
            s1uv = s1uv + sph * uac(i) * h1v * w1
            s1vu = s1vu + sph * vac(i) * h1u * w1
            s1vv = s1vv + sph * vac(i) * h1v * w1
            s1ww = s1ww + sph * wac(i) * h1w * w1
            x0 = x1
            h0u = h1u
            h0v = h1v
            h0w = h1w
        end do
        iuu = s1uu
        iuv = s1uv
        ivu = s1vu
        ivv = s1vv
        iww = s1ww
        return
    end subroutine inta1

    subroutine intb(sphj, uac, ubd, vac, vbd, wac, wbd, nod, wgt, iuu, iuv, ivu, ivv, iww)
        double complex, dimension(:, :), intent(in) :: sphj
        double complex, dimension(:), intent(in) :: uac
        double complex, dimension(:), intent(in) :: ubd
        double complex, dimension(:), intent(in) :: vac
        double complex, dimension(:), intent(in) :: vbd
        double complex, dimension(:), intent(in) :: wac
        double complex, dimension(:), intent(in) :: wbd
        double precision, dimension(:), intent(in) :: nod
        double precision, dimension(:), intent(in) :: wgt

        double complex, dimension(:), intent(out) :: iuu
        double complex, dimension(:), intent(out) :: iuv
        double complex, dimension(:), intent(out) :: ivu
        double complex, dimension(:), intent(out) :: ivv
        double complex, dimension(:), intent(out) :: iww

        double complex, dimension(size(sphj, 1)) :: spj, iu1, iv1, iw1, iu2, iv2, iw2
        double precision :: w
        integer :: m, i

        m = size(nod)
        iu1 = 0.0d0
        iu2 = 0.0d0
        iv1 = 0.0d0
        iv2 = 0.0d0
        iw1 = 0.0d0
        iw2 = 0.0d0
        do i = 1, m
            w = wgt(i)
            spj = sphj(:, i)
            iu1 = iu1 + spj * uac(i) * w
            iv1 = iv1 + spj * vac(i) * w
            iw1 = iw1 + spj * wac(i) * w
            iu2 = iu2 + spj * ubd(i) * w
            iv2 = iv2 + spj * vbd(i) * w
            iw2 = iw2 + spj * wbd(i) * w
        end do
        iuu = iu1 * iu2
        iuv = iu1 * iv2
        ivu = iv1 * iu2
        ivv = iv1 * iv2
        iww = iw1 * iw2
        return
    end subroutine intb

    subroutine zbesj_sph(z, l, kode, sphj, nz, ierr)
        double complex, intent(in) :: z
        integer, intent(in) :: l
        integer, intent(in) :: kode

        double complex, dimension(:), intent(out) :: sphj
        integer, intent(out) :: nz
        integer, intent(out) :: ierr

        double precision, dimension(size(sphj)) :: cyr, cyi
        double precision :: zr, zi, nu
        integer :: n

        zr = dble(z)
        zi = aimag(z)
        nu = l + 0.5d0
        n = size(sphj)
        call zbesj(zr, zi, nu, kode, n, cyr, cyi, nz, ierr)
        sphj = cyr + ii * cyi
        sphj = sqrt(pihalf/z) * sphj
        return
    end subroutine zbesj_sph

    subroutine zbesh1_sph(z, l, kode, sphh1, nz, ierr)
        double complex, intent(in) :: z
        integer, intent(in) :: l
        integer, intent(in) :: kode

        double complex, dimension(:), intent(out) :: sphh1
        integer, intent(out) :: nz
        integer, intent(out) :: ierr

        double precision, dimension(size(sphh1)) :: cyr, cyi
        double precision :: zr, zi, nu
        integer :: n

        zr = dble(z)
        zi = aimag(z)
        nu = l + 0.5d0
        n = size(sphh1)
        call zbesh(zr, zi, nu, kode, 1, n, cyr, cyi, nz, ierr)
        sphh1 = cyr + ii * cyi
        sphh1 = sqrt(pihalf/z) * sphh1
        return
    end subroutine zbesh1_sph

end module bspse_manypotential
