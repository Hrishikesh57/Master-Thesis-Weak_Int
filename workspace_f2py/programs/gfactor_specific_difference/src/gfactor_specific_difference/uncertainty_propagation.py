from __future__ import annotations

from dataclasses import dataclass


@dataclass
class DerivativeStabilitySummary:
    derivative_values: tuple[float, ...]
    relative_spread: float
    stable: bool


def centered_finite_difference(x_minus: float, x_plus: float, f_minus: float, f_plus: float) -> float:
    if x_plus == x_minus:
        raise ValueError("Finite-difference abscissae must differ.")
    return (f_plus - f_minus) / (x_plus - x_minus)


def secant_slope(xs: tuple[float, ...], fs: tuple[float, ...]) -> float:
    if len(xs) < 2 or len(fs) < 2:
        raise ValueError("At least two points are required for a secant slope.")
    return centered_finite_difference(xs[0], xs[-1], fs[0], fs[-1])


def relative_spread(values: tuple[float, ...]) -> float:
    if not values:
        return 0.0
    reference = max(abs(value) for value in values)
    if reference == 0.0:
        return 0.0
    return (max(values) - min(values)) / reference


def assess_derivative_stability(derivative_values: tuple[float, ...], *, tolerance_rel: float = 5.0e-2) -> DerivativeStabilitySummary:
    spread = abs(relative_spread(derivative_values))
    return DerivativeStabilitySummary(
        derivative_values=derivative_values,
        relative_spread=spread,
        stable=spread <= tolerance_rel,
    )
