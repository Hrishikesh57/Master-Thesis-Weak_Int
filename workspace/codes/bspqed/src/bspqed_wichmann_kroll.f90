module bspqed_wichmann_kroll
    use iso_fortran_env, only: error_unit
    use bspqed_function, only: function_type
    implicit none
    private

    double precision, parameter :: pi = 4.0d0 * atan(1.0d0)

    public :: wichmann_kroll_potential
    public :: wkfunc

    !+ ========================================================================
    !> @brief Abstract type for the Wichmann-Kroll potential energy function
    !>
    !> @author Halil Cakir
    !> @date 2017/12/15
    !>
    !> #### Description
    !>
    !> Definition of an abstract type for the potential energy function
    !> resulting from the first order vacuum polarization correction and three
    !> and more nuclear potential insertion. The nucleus is considered to be a
    !> point. The potential energy function is given by
    !>
    !>     V(x) = (alpha / pi) * (Z * alpha)**3 * v(x),
    !>
    !> where alpha is the fine-structure constant and Z is the nuclear charge
    !> number. This function returns the value of v for positive values of the
    !> argument x. The function v is evaluated using expressions from [1].
    !>
    !> #### References
    !>
    !> [1] A. G. Fainshtein et al., J. Phys. B 23, 559 (1990)
    !>
    !+ #### Parameters
    !>
    !+ ========================================================================
    type, extends(function_type) :: wichmann_kroll_potential
        double precision :: zalpha
    contains
        procedure, pass(self) :: eval => wichmann_kroll_potential_eval
    end type wichmann_kroll_potential

contains

    !+ ========================================================================
    !+ @brief Evaluation of the Wichmann-Kroll potential energy function
    !+
    !+ @author Halil Cakir
    !+ @date 2017/04/27
    !+ ========================================================================
    double precision function wichmann_kroll_potential_eval(x, self)
        double precision, intent(in) :: x
        class(wichmann_kroll_potential), intent(in) :: self

        double precision :: zalpha

        zalpha = self%zalpha

        if (x <= 0.0d0) then
            write(error_unit, *) "* ERROR in eval of wichmann_kroll_potential"
            write(error_unit, *) "* x must be positive."
            stop
        endif

        wichmann_kroll_potential_eval = wkfunc(x, zalpha) / x
        return
    end function wichmann_kroll_potential_eval

    !+ ========================================================================
    !+ The Wichmann-Kroll potential is given by
    !+     U(x) = (alpha / pi) * (Z * alpha)**3 * wkfunc(x, Z * alpha) / x
    !+ ========================================================================
    double precision function wkfunc(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        double precision :: nominator, denominator

        if (x < 0.0d0) then
            wkfunc = 0.0d0
        else if (x <= 0.15d0) then
            nominator = vacfunc3(x) + charge5p(zalpha)
            denominator = denom1(x, zalpha)
            wkfunc = nominator / denominator
        else if (x <= 1.5d0) then
            nominator = vacfunc3(x)
            denominator = denom2(x, zalpha)
            wkfunc = nominator / denominator
        else if (x <= 4.0d0) then
            nominator = vacfunc3(x)
            denominator = denom3(x, zalpha)
            wkfunc = nominator / denominator
        else if (x <= 10.0d0) then
            nominator = vacfunc3(x)
            denominator = denom4(x, zalpha)
            wkfunc = nominator / denominator
        else
            wkfunc = wkfunc_asy(x, zalpha)
        end if
        return
    end function wkfunc

    double precision function wkfunc_asy(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        double precision, dimension(5, 5:15), parameter :: f = reshape([&
            & 0.3812769d0, 0.0266795d0, 0.0000860d0,    0.0d+00, 0.0d+00, &
            & 0.4946161d0, 0.0458093d0, 0.0002918d0,    6.0d-08, 0.0d+00, &
            & 0.6150374d0, 0.0707397d0, 0.0007216d0,    5.0d-07, 0.0d+00, &
            & 0.7412290d0, 0.1014082d0, 0.0014705d0,   23.0d-07, 1.0d-10, &
            & 0.8721898d0, 0.1376368d0, 0.0026288d0,   69.0d-07, 2.0d-09, &
            & 1.0071432d0, 0.1791838d0, 0.0042776d0,  167.0d-07, 8.0d-09, &
            & 1.1454776d0, 0.2557770d0, 0.0064866d0,  346.0d-07, 3.0d-08, &
            & 1.2867042d0, 0.2771342d0, 0.0093134d0,  643.0d-07, 9.0d-08, &
            & 1.4304276d0, 0.3329753d0, 0.0128043d0, 1097.0d-07, 2.0d-07, &
            & 1.5763244d0, 0.3930296d0, 0.0169953d0, 1750.0d-07, 5.0d-07, &
            & 1.7241274d0, 0.4570401d0, 0.0219132d0, 2644.0d-07, 9.0d-07  &
        &], shape(f))

        double precision :: xarg, zarg, summ, term1, term2
        integer :: m

        xarg = 1.0d0 / x**2
        zarg = zalpha**2
        summ = 0.0d0
        summ = summ + 2.0d0 / 225.0d0
        summ = summ + 59.0d0 * xarg / 1323.0d0
        summ = summ + (1977.0d0 + 20.0d0 * zarg) * xarg**2 / 4725.0d0
        summ = summ + (1258680.0d0 + 34144.0d0 * zarg) * xarg**3 / 190575.0d0
        summ = summ + (1960420032.0d0 + 93618070.0d0 * zarg + 96096.0d0 * zarg**2) * xarg**4 / 12297285.0d0

        term1 = 576.0d0 * xarg**4
        do m = 5, 15
            term1 = term1 * m**2 * xarg
            term2 = polyeval(f(:, m), zarg)
            summ = summ + term1 * term2
        end do
        wkfunc_asy = summ / x**4
        return
    end function wkfunc_asy

    !+ ========================================================================
    !+ The Vacuum polarization potential to 3rd order in Z * alpha is
    !+     U^(3)(x) = (alpha / pi) * (Z * alpha)**3 * vacfunc3(x) / x
    !+ ========================================================================
    double precision function vacfunc3(x)
        double precision, intent(in) :: x

        double precision, dimension(0:35), parameter :: p = [&
            &  2.094001968517539d-02, -6.328682056509706d-02, &
            &  1.656675047469226d-01,  6.218254052794603d-02, &
            &  8.561450880149488d-01, -1.324099530525646d+00, &
            &  1.022894990055742d-01,  1.818708880809046d-01, &
            &  5.002561255243687d-04,  2.167778174657628d-01, &
            &  5.453835857080859d-01,  3.137663571113079d-01, &
            & -1.619774880547481d-01, -6.749903602516955d-02, &
            & -2.935144138820853d-04,  9.755800734714044d-02, &
            & -1.994016822313688d-01, -4.777623414220549d-02, &
            &  7.503472520308676d-03,  3.440742164594281d-05, &
            &  2.8294209807d-03,  1.4195518237d-04, &
            &  1.3489833978d-05,  2.0134008556d-07, &
            &  8.6029184536d-06, -1.7870165724d-05, &
            &  2.3269607856d-05, -1.7798463760d-05, &
            &  8.7274744996d-06, -2.9235606710d-06, &
            &  6.8929360802d-07, -1.1491699126d-07, &
            &  1.3302835807d-08, -1.0190988753d-09, &
            &  4.6508041566d-11, -9.5785924366d-13  &
        &]

        double precision :: term, xlogx, xarg

        vacfunc3 = 0.0d0
        if (x <= 4.25d0) then
            term = polyeval(p(0:8), x)
            vacfunc3 = term
            term = polyeval(p(9:14), x)
            xlogx = xlog(x)
            vacfunc3 = vacfunc3 + term * x**2 * xlogx
            term = polyeval(p(15:19), x)
            vacfunc3 = vacfunc3 + term * x**2 * xlogx**2
        else
            xarg = (10.0d0 / x)**2
            term = polyeval(p(20:35), xarg)
            vacfunc3 = term / x**4
        end if
        vacfunc3 = vacfunc3 * pi
        return
    end function vacfunc3

    !+ ========================================================================
    !+ The function -e * Q^(5+) is given by:
    !+     (alpha / pi) * (Z * alpha)**3 * charge5p(zalpha)
    !+ ========================================================================
    double precision function charge5p(zalpha)
        double precision, intent(in) :: zalpha

        double precision, dimension(15), parameter :: q = [&
            &  1.42425723499d-02,  0.76419630941d-02,  0.36124746165d-02, &
            &  0.15118743043d-02,  0.05247168246d-02,  0.01797758276d-02, &
            & -0.04761316284d-02,  0.19713924639d-02, -0.62421312791d-02, &
            &  1.29259217859d-02, -1.89507086023d-02,  1.89089237003d-02, &
            & -1.24611625668d-02,  0.48793792307d-02, -0.08858438257d-02  &
        &]

        double precision :: lambda

        lambda = lambda_function(zalpha)
        charge5p = polyeval(q, lambda) * lambda
        charge5p = charge5p * pi
        return
    end function charge5p

    double precision function denom1(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        integer, dimension(4), parameter :: s = [1, 1, 0, 0]
        double precision, dimension(0:7, 4), parameter :: a = reshape([&
            &  2.9675d+00,  7.8894d-01, -3.0201d-04,  4.9725d-04, &
            & -4.4288d+00, -2.9946d+01, -9.0441d-01,  2.4826d+00, &
            &  5.9339d+01,  8.4881d+02,  5.8388d+01, -1.0745d+02, &
            & -5.6820d+02, -1.4888d+04, -1.5899d+03,  2.7128d+03, &
            &  4.6131d+03,  1.6003d+05,  2.2880d+04, -3.8249d+04, &
            & -2.6282d+04, -1.0259d+06, -1.8036d+05,  3.0100d+05, &
            &  8.9448d+04,  3.5980d+06,  7.3716d+05, -1.2368d+06, &
            & -1.3412d+05, -5.3067d+06, -1.2228d+06,  2.0678d+06  &
        &], shape(a), order=[2, 1])

        denom1 = denom(x, zalpha, a, s)
        return
    end function denom1

    double precision function denom2(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        integer, dimension(3), parameter :: s = [0, -3, 0]
        double precision, dimension(0:6, 3), parameter :: a = reshape([&
            & -6.64673d-01,  4.68647d-05,  2.56839d-02, &
            &  5.70540d-01, -9.83836d-04, -3.21021d-03, &
            & -9.43955d-01,  8.28746d-03,  1.01163d-02, &
            &  1.18405d+00,  5.94372d-02, -1.62093d-02, &
            & -8.91689d-01, -5.31431d-02,  1.17762d-02, &
            &  3.64707d-01,  7.52653d-03, -4.69814d-03, &
            & -6.17902d-02, -4.54081d-04,  7.76469d-04  &
        &], shape(a), order=[2, 1])

        denom2 = denom(x, zalpha, a, s)
        return
    end function denom2

    double precision function denom3(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        integer, dimension(3), parameter :: s = [-1, 2, 0]
        double precision, dimension(0:6, 3), parameter :: a = reshape([&
            &  1.93743d-02,  1.64381d-01,  2.53250d-02, &
            & -6.43399d-01, -3.05149d-01,  3.56651d-03, &
            &  2.56565d-01,  2.31560d-01, -8.14980d-03, &
            & -1.09235d-01, -9.49948d-02,  6.21168d-03, &
            &  4.14226d-02,  2.22450d-02, -3.29527d-03, &
            & -7.28744d-03, -2.80168d-03,  7.60885d-04, &
            &  4.54914d-04,  1.47382d-04, -6.03993d-05  &
        &], shape(a), order=[2, 1])

        denom3 = denom(x, zalpha, a, s)
        return
    end function denom3

    double precision function denom4(x, zalpha)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha

        integer, dimension(2), parameter :: s = [-13, -6]
        double precision, dimension(0:6, 2), parameter :: a = reshape([&
            & -2.4734758d+09, -1.5115262d+04, &
            &  3.2270784d+09,  1.4945586d+04, &
            & -1.7207915d+09, -5.8361389d+03, &
            &  4.7635740d+08,  1.1476935d+03, &
            & -7.1414759d+07, -1.2098290d+02, &
            &  5.4305859d+06,  6.5439613d+00, &
            & -1.6479928d+05, -1.4296421d-01  &
        &], shape(a), order=[2, 1])

        denom4 = denom(x, zalpha, a, s)
        return
    end function denom4

    double precision function denom(x, zalpha, a, s)
        double precision, intent(in) :: x
        double precision, intent(in) :: zalpha
        double precision, dimension(:, :), intent(in) :: a
        integer, dimension(:), intent(in) :: s
        double precision :: lambda, summ, lmbn, term
        integer :: nn, n
        nn = size(s)
        lambda = lambda_function(zalpha)
        summ = 0.0d0
        lmbn = 1.0d0
        do n = 1, nn
            lmbn = lmbn * lambda
            term = polyeval(a(:, n), x)
            term = term * x**s(n) * lmbn
            summ = summ + term
        end do
        denom = 1.0d0 + summ
        return
    end function denom

    double precision function polyeval(coefs, x)
        double precision, dimension(:), intent(in) :: coefs
        double precision, intent(in) :: x
        double precision :: y
        integer :: n, k
        n = size(coefs)
        y = coefs(n)
        do k = n-1, 1, -1
            y = y * x + coefs(k)
        end do
        polyeval = y
        return
    end function polyeval

    double precision function xlog(x)
        double precision, intent(in) :: x
        if (x == 0.0d0) then
            xlog = 0.0d0
            return
        else
            xlog = x * log(x)
            return
        end if
    end function xlog

    double precision function lambda_function(x)
        double precision, intent(in) :: x
        double precision :: z
        if (abs(x) < 0.005d0) then
            z = x**2
            lambda_function = 0.5d0 + 0.125d0 * z * (1.0d0 + 0.5d0 * z * &
            & (1.0d0 + 0.625d0 * z))
            lambda_function =  lambda_function * z
            return
        else
            lambda_function = 1.0d0 - sqrt((1.0d0 - x) * (1.0d0 + x))
            return
        end if
    end function lambda_function

end module bspqed_wichmann_kroll
