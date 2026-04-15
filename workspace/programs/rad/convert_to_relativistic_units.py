#!/usr/bin/env python3
"""
Script to convert atomic radial function data from atomic units to
relativistic units. Reads .dat files and converts r, G, F data according
to relativistic scaling.
"""

import numpy as np
import sys
import os
from pathlib import Path

RAD_ROOT = Path(__file__).resolve().parent
SEARCH_DIRS = [
    RAD_ROOT / "data" / "reference" / "radial_functions",
    RAD_ROOT / "data" / "grasp" / "raw",
    RAD_ROOT,
]


def resolve_input_path(input_file):
    candidate = Path(input_file)
    if candidate.exists():
        return candidate

    for directory in SEARCH_DIRS:
        test_path = directory / candidate.name
        if test_path.exists():
            return test_path

    return candidate


def convert_to_relativistic_units(input_file, output_base=None):
    """
    Convert atomic units data to relativistic units.
    
    Parameters:
    -----------
    input_file : str
        Path to input .dat file
    output_base : str, optional
        Base name for output files (without extension)
        If None, uses input filename
    """
    
    # Physical constants
    alpha_inv = 137.035999084  # CODATA 2024 fine structure constant inverse
    
    # Read input file and handle Fortran double precision format
    try:
        # First, read the file as text and convert Fortran 'D' format to 'E' format
        input_path = resolve_input_path(input_file)
        with open(input_path, 'r') as f:
            lines = f.readlines()
        
        # Convert Fortran double precision format (D+xx) to standard format (E+xx)
        converted_lines = []
        for line in lines:
            # Replace 'D+' and 'D-' with 'E+' and 'E-'
            converted_line = line.replace('D+', 'E+').replace('D-', 'E-')
            converted_lines.append(converted_line)
        
        # Try to parse the converted data
        try:
            # First try skipping the first row (in case of header)
            data_lines = converted_lines[1:]
            data = np.array([list(map(float, line.split())) for line in data_lines if line.strip()])
            print(f"Successfully read {len(data)} data points from {input_path} (skipped header)")
        except (ValueError, IndexError):
            # If that fails, try without skipping rows
            data = np.array([list(map(float, line.split())) for line in converted_lines if line.strip()])
            print(f"Successfully read {len(data)} data points from {input_path}")
            
    except Exception as e:
        print(f"Error reading file: {e}")
        raise
    
    if data.shape[1] < 3:
        msg = (f"Input file must have at least 3 columns (r, g, f), "
               f"found {data.shape[1]}")
        raise ValueError(msg)
    
    # Extract columns - assuming r, g_nk(r), f_nk(r)
    r = data[:, 0]  # radial coordinate in atomic units
    G = data[:, 1]  # g_nk(r) function
    F = data[:, 2]  # f_nk(r) function
    
    print("Original data ranges:")
    print(f"  r: {r.min():.6e} to {r.max():.6e} a.u.")
    print(f"  G: {G.min():.6e} to {G.max():.6e}")
    print(f"  F: {F.min():.6e} to {F.max():.6e}")
    
    # Convert to relativistic units
    r_rel = r * alpha_inv
    scale = 1 / np.sqrt(alpha_inv)
    G_rel = G * scale
    F_rel = F * scale
    
    print("\nAfter relativistic conversion:")
    print(f"  r_rel: {r_rel.min():.6e} to {r_rel.max():.6e}")
    print(f"  G_rel: {G_rel.min():.6e} to {G_rel.max():.6e}")
    print(f"  F_rel: {F_rel.min():.6e} to {F_rel.max():.6e}")
    print(f"  Scale factor: {scale:.6e}")
    print(f"  Alpha inverse: {alpha_inv}")
    
    # Determine output filenames
    if output_base is None:
        output_base = input_path.stem + "_relunits"
        output_dir = input_path.parent
    else:
        output_dir = Path(output_base).parent
        output_base = Path(output_base).stem
    
    # Save as npz file
    npz_file = output_dir / f"{output_base}.npz"
    np.savez(npz_file, r=r_rel, G=G_rel, F=F_rel)
    print(f"\nSaved NPZ file: {npz_file}")
    
    # Save as txt file with 3 columns: r_rel, G_rel, F_rel
    txt_file = output_dir / f"{output_base}.txt"
    header = "r_rel (rel. units)    G_rel(r)          F_rel(r)"
    
    # Format data for text output
    with open(txt_file, 'w') as f:
        f.write(header + '\n')
        for i in range(len(r_rel)):
            f.write(f"{r_rel[i]:15.6e} {G_rel[i]:15.6e} {F_rel[i]:15.6e}\n")
    
    print(f"Saved TXT file: {txt_file}")
    
    return r_rel, G_rel, F_rel, npz_file, txt_file


def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) < 2:
        usage_msg = ("Usage: python convert_to_relativistic_units.py "
                     "<input_file.dat> [output_base]")
        print(usage_msg)
        print("\nExample:")
        print("  python convert_to_relativistic_units.py data/reference/radial_functions/radial_function.dat")
        example_msg = ("  python convert_to_relativistic_units.py "
                       "data/reference/radial_functions/radial_function.dat my_output")
        print(example_msg)
        print("\nThis will create:")
        print("  - output_base.npz (numpy array file)")
        print("  - output_base.txt (3-column text file: r_rel, G_rel, F_rel)")
        return
    
    input_file = sys.argv[1]
    output_base = sys.argv[2] if len(sys.argv) > 2 else None

    resolved_input = resolve_input_path(input_file)
    if not resolved_input.exists():
        print(f"Error: Input file '{input_file}' not found!")
        return
    
    try:
        result = convert_to_relativistic_units(str(resolved_input), output_base)
        r_rel, G_rel, F_rel, npz_file, txt_file = result
        print("\nConversion completed successfully!")
        print("Files created:")
        print(f"  {npz_file}")
        print(f"  {txt_file}")
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        return


if __name__ == "__main__":
    main()
