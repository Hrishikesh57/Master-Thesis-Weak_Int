import numpy as np

def calculate_rho0(R, Cv_p, Z, N):
    # Physical constants in MeV·fm units
    hbar_c = 137  # MeV·fm

    prefactor = 3 / (4 * np.pi)
    ratio = (hbar_c / R) ** 3
    rho0 = prefactor * ratio * (Cv_p * Z - N)
    return rho0

# Example usage:
R = 1.1069e-4        # Nuclear radius in fm
Cv_p = 0.08       # Example value for C_v^p
Z = 92            # Atomic number
N = 146           # Number of neutrons

rho0 = calculate_rho0(R, Cv_p, Z, N)
print(f"rho0 = {rho0:.6e}")