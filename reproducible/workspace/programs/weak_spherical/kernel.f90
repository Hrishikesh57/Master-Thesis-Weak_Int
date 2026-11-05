! Source-complete central-potential DKB interface. All radial coordinates are in
! hbar/(m_e*c), energies in m_e*c^2, with mu=m_e/m_e(reference).
! Spinors: (G Omega_kappa, i F Omega_-kappa)/r; integral (G^2+F^2) dr=1.
subroutine solve_spherical_basis(z, charge_c_fm, charge_a_fm, npoints, order, nq, &
    cavity_bohr, mu, kappa, energies, r, w, g, f, ierr)
    use bspqed_dkb, only: solve, radialgf
    use bspqed_physconsts, only: alpha, electron_reduced_compton
    use bspqed_potentials, only: fermi_potential
    use bspqed_quadrature, only: nodes20, weights20, nodes40, weights40
    implicit none
    integer, intent(in) :: z, npoints, order, nq, kappa
    double precision, intent(in) :: charge_c_fm, charge_a_fm, cavity_bohr, mu
    double precision, intent(out) :: energies(2*(npoints+order-4))
    double precision, intent(out) :: r(nq*(npoints-1)), w(nq*(npoints-1))
    double precision, intent(out) :: g(nq*(npoints-1),2*(npoints+order-4))
    double precision, intent(out) :: f(nq*(npoints-1),2*(npoints+order-4))
    integer, intent(out) :: ierr
    double precision :: knots(npoints+2*(order-1)), points(npoints), nodes(nq), weights(nq)
    double precision :: coefs(2*(npoints+order-2),2*(npoints+order-4))
    double precision :: rmesh, cavity, left, half, phase
    integer :: ninside, i, j, k, idx, iv, ref
    type(fermi_potential) :: potential
    energies=0; r=0; w=0; g=0; f=0; ierr=1
    if (z<=0 .or. npoints<40 .or. order<4 .or. kappa==0) return
    if (nq/=20 .and. nq/=40) return
    if (charge_c_fm<=0 .or. charge_a_fm<=0 .or. mu<=0 .or. cavity_bohr<=0) return
    potential%zalpha=dble(z)*alpha
    potential%c=charge_c_fm*mu/electron_reduced_compton
    potential%a=charge_a_fm*mu/electron_reduced_compton
    ! Fixed physical cavity under mass differentiation, as for all nuclear lengths.
    cavity=cavity_bohr*mu/alpha
    rmesh=potential%c+5*potential%a
    if (cavity<=rmesh) return
    ninside=npoints/10
    ! Same linear/interior + exponential/exterior construction as kntseq1.
    ! Keep ALL repeated origin knots zero. Gauss nodes never touch the endpoint.
    do i=1,ninside
        points(i)=rmesh*dble(i-1)/dble(ninside-1)
    end do
    do i=1,npoints-ninside
        points(ninside+i)=rmesh*(cavity/rmesh)**(dble(i)/dble(npoints-ninside))
    end do
    knots(1:order-1)=0
    knots(order:order+npoints-1)=points
    knots(order+npoints:)=cavity
    if (any(knots(2:)<knots(:size(knots)-1))) return
    if (nq==20) then
        nodes=nodes20; weights=weights20
    else
        nodes=nodes40; weights=weights40
    end if
    call solve(knots,order,kappa,potential,nodes,weights,energies,coefs)
    do i=1,npoints-1
        left=points(i); half=(points(i+1)-left)/2
        do j=1,nq
            idx=(i-1)*nq+j
            r(idx)=left+half*(nodes(j)+1)
            w(idx)=half*weights(j)
            do k=1,size(energies)
                iv=order-1+i
                call radialgf(knots,order,kappa,coefs(:,k),r(idx),g(idx,k),f(idx,k),iv)
            end do
        end do
    end do
    ! Whole-state rephasing only; do not change the relative G/F sign.
    do k=1,size(energies)
        ref=maxloc(abs(g(:,k)),dim=1)
        phase=sign(1.0d0,g(ref,k))
        g(:,k)=phase*g(:,k); f(:,k)=phase*f(:,k)
    end do
    ierr=0
end subroutine solve_spherical_basis

subroutine spherical_constants(alpha_out,lambda_fm,gf_rel,cve_out,cvp_out)
    use bspqed_physconsts, only: alpha,electron_reduced_compton,Gfermi_relativistic,Cv_e,Cv_p
    implicit none
    double precision,intent(out)::alpha_out,lambda_fm,gf_rel,cve_out,cvp_out
    alpha_out=alpha; lambda_fm=electron_reduced_compton
    gf_rel=Gfermi_relativistic; cve_out=Cv_e; cvp_out=Cv_p
end subroutine spherical_constants

subroutine charge_potential(z,c_fm,a_fm,r_fm,v)
    use bspqed_physconsts, only: alpha,electron_reduced_compton
    use bspqed_potentials, only: fermi_potential
    implicit none
    integer,intent(in)::z
    double precision,intent(in)::c_fm,a_fm,r_fm
    double precision,intent(out)::v
    type(fermi_potential)::potential
    potential%zalpha=z*alpha
    potential%c=c_fm/electron_reduced_compton; potential%a=a_fm/electron_reduced_compton
    v=potential%eval(r_fm/electron_reduced_compton)
end subroutine charge_potential
