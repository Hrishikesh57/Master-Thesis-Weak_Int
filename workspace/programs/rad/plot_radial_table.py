import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

NEAR_ORIGIN_RMAX = 2.0
TAIL_RELATIVE_THRESHOLD = 1.0e-8
TAIL_BUFFER_POINTS = 3
RAD_ROOT = Path(__file__).resolve().parent
PLOT_OUTPUT_DIR = RAD_ROOT / "outputs" / "plots"
SEARCH_DIRS = [
    RAD_ROOT / "outputs" / "wavefunctions" / "single_electron",
    RAD_ROOT / "outputs" / "wavefunctions" / "scf",
    RAD_ROOT / "data" / "reference" / "wavefunctions" / "single_electron",
    RAD_ROOT / "data" / "reference" / "wavefunctions" / "single_electron_legacy",
    RAD_ROOT / "data" / "reference" / "wavefunctions" / "scf",
    RAD_ROOT / "data" / "reference" / "wavefunctions" / "scf_legacy",
    RAD_ROOT / "data" / "grasp" / "rel",
    RAD_ROOT / "data" / "reference" / "radial_functions",
    RAD_ROOT,
]

def ensure_plot_output_dir():
    PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def resolve_input_path(filename):
    candidate = Path(filename)
    if candidate.exists():
        return candidate.resolve()

    if not candidate.is_absolute():
        for directory in SEARCH_DIRS:
            test_path = directory / candidate.name
            if test_path.exists():
                return test_path.resolve()

    return candidate

def compute_auto_plot_rmax(r, g, f):
    """
    Choose the linear-plot cutoff from the radial probability envelope.
    Stop once the tail maximum of G^2+F^2 stays below a fixed fraction
    of its peak, so nodes in G or F do not trigger premature cutoffs.
    """
    rho = g**2 + f**2
    peak_rho = float(np.max(rho))
    if peak_rho <= 0.0:
        return float(np.max(r)), False, 0.0

    threshold = TAIL_RELATIVE_THRESHOLD * peak_rho
    tail_max = np.maximum.accumulate(rho[::-1])[::-1]
    below_threshold = np.where(tail_max < threshold)[0]

    if len(below_threshold) == 0:
        return float(np.max(r)), False, threshold

    cutoff_index = min(int(below_threshold[0]) + TAIL_BUFFER_POINTS, len(r) - 1)
    return float(r[cutoff_index]), True, threshold

def analyze_data(filename):
    """
    Analyze the data to check for potential issues
    """
    print(f"\nAnalyzing data from: {filename}")
    print("=" * 50)
    
    try:
        data = []
        with open(filename, 'r') as file:
            for line_num, line in enumerate(file, 1):
                line = line.strip()
                if not line.startswith('#') and line:
                    try:
                        # Convert Fortran scientific notation (D) to Python format (E)
                        line = line.replace('D', 'E')
                        values = [float(x) for x in line.split()]
                        if len(values) == 3:
                            data.append(values)
                    except ValueError:
                        print(f"Warning: Could not parse line {line_num}: {line}")
                        continue
        
        if not data:
            print("Error: No valid data found!")
            return
            
        data = np.array(data)
        r = data[:, 0]
        g = data[:, 1]
        f = data[:, 2]
        
        print(f"Data points: {len(data)}")
        print(f"r range: {r.min():.6e} to {r.max():.6e}")
        print(f"g(r) range: {g.min():.6e} to {g.max():.6e}")
        print(f"f(r) range: {f.min():.6e} to {f.max():.6e}")
        print(f"|g(r)| max: {np.abs(g).max():.6e}")
        print(f"|f(r)| max: {np.abs(f).max():.6e}")
        
        # Check for unrealistic values
        if np.abs(g).max() > 1000 or np.abs(f).max() > 1000:
            print("WARNING: Extremely large wavefunction values detected!")
            print("This suggests the data may be incorrect.")
        
        # Check normalization with the actual radial grid
        prob_density = g**2 + f**2
        norm_estimate = np.trapezoid(prob_density, r)
        print(f"Normalization check: {norm_estimate:.6e}")
        if norm_estimate > 100:
            print("WARNING: Normalization seems too large!")
        
        # Check for NaN or infinite values
        if np.any(np.isnan(g)) or np.any(np.isnan(f)):
            print("WARNING: NaN values detected!")
        if np.any(np.isinf(g)) or np.any(np.isinf(f)):
            print("WARNING: Infinite values detected!")
            
        return data
        
    except Exception as e:
        print(f"Error analyzing data: {e}")
        return None

def plot_wavefunction_with_diagnostics(filename):
    """
    Plot radial wavefunction with diagnostic information (linear scale only)
    """
    
    # First analyze the data
    resolved_path = resolve_input_path(filename)
    data = analyze_data(str(resolved_path))
    if data is None:
        return
    
    r = data[:, 0]
    g = data[:, 1]
    f = data[:, 2]
    
    # Extract state info from filename
    base_name = resolved_path.name
    parts = base_name.split('_')
    if len(parts) >= 4:
        z_str = parts[0]
        state = parts[3].replace('.dat', '').replace('.txt', '')
        title_info = f"Z={z_str}, State={state}"
    else:
        title_info = base_name.replace('.dat', '').replace('.txt', '')
    
    # Create plot with 2 subplots: linear scale and near origin
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    plot_r_max, tail_reached, threshold = compute_auto_plot_rmax(r, g, f)
    if tail_reached:
        print(f"Auto tail cutoff: rho(r) stays below {TAIL_RELATIVE_THRESHOLD:.0e} of peak after r = {plot_r_max:.6e}")
    else:
        print(f"Tail threshold {TAIL_RELATIVE_THRESHOLD:.0e} not reached inside tabulated range; plotting to r = {plot_r_max:.6e}")
    
    # Plot 1: Linear scale over the default physically relevant range
    ax1.plot(r, g, 'b-', linewidth=2, label='g(r) - Large component')
    ax1.plot(r, f, 'r-', linewidth=2, label='f(r) - Small component')
    ax1.set_xlabel('r (relativistic units, m_e = 1)')
    ax1.set_ylabel('Wavefunction amplitude')
    ax1.set_title(f'Linear Scale (auto r ≤ {plot_r_max:.1f}) - {title_info}')
    ax1.set_xlim(0.0, plot_r_max)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Plot 2: Close-up near origin
    r_max_closeup = min(NEAR_ORIGIN_RMAX, plot_r_max)
    mask_closeup = r <= r_max_closeup
    if np.any(mask_closeup):
        ax2.plot(r[mask_closeup], g[mask_closeup], 'b-', linewidth=2, label='g(r)')
        ax2.plot(r[mask_closeup], f[mask_closeup], 'r-', linewidth=2, label='f(r)')
    ax2.set_xlabel('r (relativistic units, m_e = 1)')
    ax2.set_ylabel('Wavefunction amplitude')
    ax2.set_title(f'Near Origin (r ≤ {r_max_closeup}) - {title_info}')
    ax2.set_xlim(0.0, r_max_closeup)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save plot
    ensure_plot_output_dir()
    output_name = PLOT_OUTPUT_DIR / f"{resolved_path.stem}_diagnostic_plot.png"
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"\nDiagnostic plot saved as: {output_name}")
    
    # Show plot
    plt.show()

def compare_single_vs_multi_electron(single_file, multi_file):
    """
    Compare single electron vs multi-electron data (linear scale only)
    """
    print("\nComparing single electron vs multi-electron data:")
    print("=" * 60)
    
    # Load both datasets
    single_path = resolve_input_path(single_file)
    multi_path = resolve_input_path(multi_file)
    single_data = analyze_data(str(single_path))
    multi_data = analyze_data(str(multi_path))
    
    if single_data is None or multi_data is None:
        return
    
    r_single = single_data[:, 0]
    g_single = single_data[:, 1]
    f_single = single_data[:, 2]
    
    r_multi = multi_data[:, 0]
    g_multi = multi_data[:, 1]
    f_multi = multi_data[:, 2]
    plot_r_single, _, _ = compute_auto_plot_rmax(r_single, g_single, f_single)
    plot_r_multi, _, _ = compute_auto_plot_rmax(r_multi, g_multi, f_multi)
    plot_r_max = max(plot_r_single, plot_r_multi)
    
    # Create comparison plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Compare g(r) components
    ax1.plot(r_single, g_single, 'b-', linewidth=2, label='Single electron g(r)')
    ax1.plot(r_multi, g_multi, 'r--', linewidth=2, label='Multi-electron g(r)')
    ax1.set_xlabel('r (relativistic units, m_e = 1)')
    ax1.set_ylabel('g(r)')
    ax1.set_title(f'Large Component Comparison (auto r ≤ {plot_r_max:.1f})')
    ax1.set_xlim(0.0, plot_r_max)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Compare f(r) components
    ax2.plot(r_single, f_single, 'b-', linewidth=2, label='Single electron f(r)')
    ax2.plot(r_multi, f_multi, 'r--', linewidth=2, label='Multi-electron f(r)')
    ax2.set_xlabel('r (relativistic units, m_e = 1)')
    ax2.set_ylabel('f(r)')
    ax2.set_title(f'Small Component Comparison (auto r ≤ {plot_r_max:.1f})')
    ax2.set_xlim(0.0, plot_r_max)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Close-up near origin for g(r)
    r_max_closeup = min(NEAR_ORIGIN_RMAX, plot_r_max)
    mask_single_closeup = r_single <= r_max_closeup
    mask_multi_closeup = r_multi <= r_max_closeup
    
    if np.any(mask_single_closeup):
        ax3.plot(r_single[mask_single_closeup], g_single[mask_single_closeup], 'b-', 
                linewidth=2, label='Single electron g(r)')
    if np.any(mask_multi_closeup):
        ax3.plot(r_multi[mask_multi_closeup], g_multi[mask_multi_closeup], 'r--', 
                linewidth=2, label='Multi-electron g(r)')
    ax3.set_xlabel('r (relativistic units, m_e = 1)')
    ax3.set_ylabel('g(r)')
    ax3.set_title(f'Large Component Near Origin (r ≤ {r_max_closeup})')
    ax3.set_xlim(0.0, r_max_closeup)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Close-up near origin for f(r)
    if np.any(mask_single_closeup):
        ax4.plot(r_single[mask_single_closeup], f_single[mask_single_closeup], 'b-', 
                linewidth=2, label='Single electron f(r)')
    if np.any(mask_multi_closeup):
        ax4.plot(r_multi[mask_multi_closeup], f_multi[mask_multi_closeup], 'r--', 
                linewidth=2, label='Multi-electron f(r)')
    ax4.set_xlabel('r (relativistic units, m_e = 1)')
    ax4.set_ylabel('f(r)')
    ax4.set_title(f'Small Component Near Origin (r ≤ {r_max_closeup})')
    ax4.set_xlim(0.0, r_max_closeup)
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    ensure_plot_output_dir()
    output_name = PLOT_OUTPUT_DIR / f"comparison_{single_path.stem}_vs_{multi_path.stem}.png"
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved as: {output_name}")
    plt.show()

def plot_wavefunction(filename):
    """
    Plot radial wavefunction components g(r) and f(r) (linear scale only)
    Data format: .dat file with columns: r, g(r), f(r)
    """
    
    # Check if file exists
    resolved_path = resolve_input_path(filename)
    if not resolved_path.exists():
        print(f"Error: File '{filename}' not found!")
        return
    
    # Read the data file, skipping comment lines
    try:
        # Read data while skipping lines that start with '#'
        data = []
        with open(resolved_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('#') and line:  # Skip comments and empty lines
                    try:
                        # Convert Fortran scientific notation (D) to Python format (E)
                        line = line.replace('D', 'E')
                        values = [float(x) for x in line.split()]
                        if len(values) == 3:  # Ensure we have r, g, f
                            data.append(values)
                    except ValueError:
                        continue  # Skip lines that can't be parsed as numbers
        
        if not data:
            print("Error: No valid data found in file!")
            return
            
        data = np.array(data)
        r = data[:, 0]  # radial coordinate
        g = data[:, 1]  # large component g(r)
        f = data[:, 2]  # small component f(r)
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Extract state info from filename
    base_name = resolved_path.name
    parts = base_name.split('_')
    if len(parts) >= 4:
        z_str = parts[0]
        state = parts[3].replace('.dat', '').replace('.txt', '')
        title_info = f"Z={z_str}, State={state}"
    else:
        title_info = base_name.replace('.dat', '').replace('.txt', '')
    
    # Create the plot with linear scale
    plt.figure(figsize=(12, 8))
    plot_r_max, tail_reached, threshold = compute_auto_plot_rmax(r, g, f)
    if tail_reached:
        print(f"Auto tail cutoff: rho(r) stays below {TAIL_RELATIVE_THRESHOLD:.0e} of peak after r = {plot_r_max:.6e}")
    else:
        print(f"Tail threshold {TAIL_RELATIVE_THRESHOLD:.0e} not reached inside tabulated range; plotting to r = {plot_r_max:.6e}")
    
    # Plot both components on the same plot
    plt.plot(r, g, 'b-', linewidth=2, label='g(r) - Large component')
    plt.plot(r, f, 'r-', linewidth=2, label='f(r) - Small component')
    
    plt.xlabel('r (relativistic units, m_e = 1)')
    plt.ylabel('Wavefunction amplitude')
    plt.title(f'Radial Wavefunction Components (auto r ≤ {plot_r_max:.1f}) - {title_info}')
    plt.xlim(0.0, plot_r_max)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save plot
    ensure_plot_output_dir()
    output_name = PLOT_OUTPUT_DIR / f"{resolved_path.stem}_linear_plot.png"
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Linear plot saved as: {output_name}")
    
    # Show plot
    plt.show()

def main():
    """
    Main function - can be run interactively or with command line arguments
    Expects .dat files with format: r g(r) f(r)
    """
    
    if len(sys.argv) > 1:
        # Command line usage
        filename = sys.argv[1]
        plot_type = sys.argv[2] if len(sys.argv) > 2 else 'diagnostic'
    else:
        # Interactive usage
        print("Radial Wavefunction Plotter (Linear Scale)")
        print("Data format: .txt/.dat files with columns r, g(r), f(r)")
        print("=" * 50)
        
        filename = input("Enter wavefunction filename (.txt or .dat): ")
        
        print("\nPlot options:")
        print("1. Diagnostic analysis (linear + near origin)")
        print("2. Simple linear plot")
        print("3. Compare with another file")
        
        choice = input("Choose plot type (1-3): ")
        
        if choice == '1':
            plot_type = 'diagnostic'
        elif choice == '2':
            plot_type = 'linear'
        elif choice == '3':
            plot_type = 'compare'
        else:
            plot_type = 'diagnostic'
    
    # Generate plots based on choice
    if plot_type == 'diagnostic':
        plot_wavefunction_with_diagnostics(filename)
    elif plot_type == 'compare':
        if len(sys.argv) > 2:
            other_file = sys.argv[2]
        else:
            other_file = input("Enter second filename for comparison: ")
        compare_single_vs_multi_electron(filename, other_file)
    elif plot_type == 'linear':
        plot_wavefunction(filename)
    else:
        print(f"Unknown plot type: {plot_type}")

if __name__ == "__main__":
    main()
