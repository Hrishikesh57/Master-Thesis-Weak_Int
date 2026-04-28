from __future__ import annotations

from .direct_modules import (
    load_hartree_suite_module,
    load_lithium_suite_module,
    load_rad_kernel_suite_module,
    load_selfen_suite_module,
)


def compute_lithium_dirac_tables(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_dirac_gfactor_tables(znucl, rms_fm, order)


def compute_lithium_delbrueck_tables(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_delbrueck_gfactor_tables(znucl, rms_fm, order)


def compute_lithium_propagator_loop_feynman_table(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_propagator_loop_feynman_table(znucl, rms_fm, order)


def compute_lithium_uehling_tables(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_uehling_gfactor_tables(znucl, rms_fm, order)


def compute_lithium_uehling_all_tables(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_uehling_all_tables(znucl, rms_fm, order)


def compute_lithium_wichmann_kroll_tables(znucl: int, rms_fm: float, order: int = 9):
    module = load_lithium_suite_module()
    return module.gfactor_lithium_suite_f2py.compute_wichmann_kroll_tables(znucl, rms_fm, order)


def compute_selfen_zero(znucl: int, rms_fm: float, kappa: int, nmain: int, npoints: int, epsabs: float, epsrel: float):
    module = load_selfen_suite_module()
    return module.selfen_suite_f2py.compute_zero_potential_self_energy(znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel)


def compute_selfen_one(znucl: int, rms_fm: float, kappa: int, nmain: int, npoints: int, epsabs: float, epsrel: float):
    module = load_selfen_suite_module()
    return module.selfen_suite_f2py.compute_one_potential_self_energy(znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel)


def compute_selfen_many(
    znucl: int,
    rms_fm: float,
    kappa: int,
    nmain: int,
    npoints: int,
    epsabs: float,
    epsrel: float,
    modkapilow: int,
    modkapiupp: int,
    modkapicount: int,
):
    module = load_selfen_suite_module()
    return module.selfen_suite_f2py.compute_many_potential_self_energy(
        znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel, modkapilow, modkapiupp, modkapicount
    )


def compute_hartree_dhfcb(znucl: int, rms_fm: float, core_nmains, core_kappas, valence_nmains, valence_kappas):
    module = load_hartree_suite_module()
    return module.hartree_dhfcb_f2py.compute_dhfcb_energies(
        znucl, rms_fm, core_nmains, core_kappas, valence_nmains, valence_kappas
    )


def compute_hartree_shf(znucl: int, core_nmains, core_orbitals, valence_nmains, valence_orbitals):
    module = load_hartree_suite_module()
    return module.hartree_shf_f2py.compute_shf_energies(
        znucl, core_nmains, core_orbitals, valence_nmains, valence_orbitals
    )


def compute_rad_weak_potential_summary(znucl: int, nneutrons: int, r_fm: float):
    module = load_rad_kernel_suite_module()
    return module.rad_kernel_f2py.weak_potential_summary(znucl, nneutrons, r_fm)


def compute_rad_first_order_weak_shift(znucl: int, nneutrons: int, r_fm: float, r, g, f):
    module = load_rad_kernel_suite_module()
    return module.rad_kernel_f2py.first_order_weak_shift(znucl, nneutrons, r_fm, r, g, f)


def compute_rad_second_order_shift(
    znucl: int,
    nneutrons: int,
    r_fm: float,
    e_initial: float,
    e_intermediate: float,
    r_intermediate,
    g_intermediate,
    f_intermediate,
    r_initial,
    g_initial,
    f_initial,
):
    module = load_rad_kernel_suite_module()
    return module.rad_kernel_f2py.second_order_shift_sampled(
        znucl,
        nneutrons,
        r_fm,
        e_initial,
        e_intermediate,
        r_intermediate,
        g_intermediate,
        f_intermediate,
        r_initial,
        g_initial,
        f_initial,
    )


def compute_rad_discrete_amplitude(
    znucl: int,
    nneutrons: int,
    r_fm: float,
    e_initial: float,
    e_final: float,
    e_intermediate: float,
    e0_field: float,
    state_factor: float,
    weak_enabled: int,
    r_initial,
    g_initial,
    f_initial,
    r_intermediate,
    g_intermediate,
    f_intermediate,
    r_final,
    g_final,
    f_final,
):
    module = load_rad_kernel_suite_module()
    return module.rad_kernel_f2py.discrete_amp_sampled(
        znucl,
        nneutrons,
        r_fm,
        e_initial,
        e_final,
        e_intermediate,
        e0_field,
        state_factor,
        weak_enabled,
        r_initial,
        g_initial,
        f_initial,
        r_intermediate,
        g_intermediate,
        f_intermediate,
        r_final,
        g_final,
        f_final,
    )
