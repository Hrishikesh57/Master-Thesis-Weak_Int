#!/usr/bin/env python3
"""
Check normalization of relativistic wavefunctions.
Calculates ∫₀^∞ (G² + F²) dr for Dirac spinors.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from pathlib import Path

RAD_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = RAD_ROOT / "data" / "grasp" / "rel" / "6_2s_rel.txt"
PLOT_OUTPUT_DIR = RAD_ROOT / "outputs" / "plots"

def check_normalization(filename):
    """
    Check normalization of a wavefunction file.
    """
    print(f"Analyzing: {filename}")
    print("=" * 50)
    
    try:
        # Read data (skip header line)
        data = np.loadtxt(filename, skiprows=1)
        
        r = data[:, 0]  # Radial coordinate
        G = data[:, 1]  # Large component G(r)
        F = data[:, 2]  # Small component F(r)
        
        print(f"Data points: {len(r)}")
        print(f"r range: {r[0]:.6e} to {r[-1]:.6e}")
        print(f"Max |G(r)|: {np.max(np.abs(G)):.6e}")
        print(f"Max |F(r)|: {np.max(np.abs(F)):.6e}")
        
        # Calculate integrand
        integrand = G**2 + F**2
        
        # Calculate normalization using different methods
        norm_simpson = simpson(integrand, x=r)
        norm_trapz = np.trapz(integrand, x=r)
        
        print(f"\nNormalization Results:")
        print(f"Simpson's rule: {norm_simpson:.10e}")
        print(f"Trapezoidal:    {norm_trapz:.10e}")
        print(f"Expected:       1.000000e+00")
        print(f"Error (Simpson): {abs(norm_simpson - 1.0):.10e}")
        print(f"Error (Trapz):   {abs(norm_trapz - 1.0):.10e}")
        
        # Check normalization quality
        if abs(norm_simpson - 1.0) < 1e-6:
            print("✓ Excellent normalization!")
        elif abs(norm_simpson - 1.0) < 1e-3:
            print("⚠ Acceptable normalization")
        else:
            print("✗ Poor normalization - possible issues:")
            print("  - Unit conversion problem")
            print("  - Missing normalization factor")
            print("  - Integration range truncation")
            
            # Check if we need scaling
            suggested_factor = 1.0 / np.sqrt(norm_simpson)
            print(f"  - Try scaling by factor: {suggested_factor:.6e}")
        
        # Create plots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
        
        # Plot G(r) and F(r)
        ax1.plot(r, G, 'b-', label='G(r) - Large component', linewidth=2)
        ax1.plot(r, F, 'r-', label='F(r) - Small component', linewidth=2)
        ax1.set_xlabel('r (relativistic units)')
        ax1.set_ylabel('Wavefunction amplitude')
        ax1.set_title('Wavefunction Components')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(left=0)
        
        # Plot G² and F²
        ax2.plot(r, G**2, 'b-', label='G²(r)', linewidth=2)
        ax2.plot(r, F**2, 'r-', label='F²(r)', linewidth=2)
        ax2.set_xlabel('r (relativistic units)')
        ax2.set_ylabel('Squared components')
        ax2.set_title('Squared Components')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(left=0)
        
        # Plot total probability density
        ax3.plot(r, integrand, 'g-', label='G²+F² - Probability density', 
                linewidth=2)
        ax3.fill_between(r, integrand, alpha=0.3, color='green')
        ax3.set_xlabel('r (relativistic units)')
        ax3.set_ylabel('G² + F²')
        ax3.set_title(f'Total Probability Density (∫ = {norm_simpson:.6e})')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim(left=0)
        
        plt.tight_layout()
        
        # Save plot
        PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        plot_name = PLOT_OUTPUT_DIR / (Path(filename).stem + "_normalization_check.png")
        plt.savefig(plot_name, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {plot_name}")
        
        plt.show()
        
        return norm_simpson
        
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    # Check the 2s wavefunction
    check_normalization(str(DEFAULT_INPUT))
