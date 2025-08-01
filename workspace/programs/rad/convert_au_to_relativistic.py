#!/usr/bin/env python3
"""
Script to convert energies from atomic units (Hartree) to relativistic units
where c = 1, m_e = 1, hbar = 1.

Conversion formula: E_relativistic = E_au * alpha^2
where alpha is the fine structure constant ≈ 7.297353e-3
"""

import sys
import os
import argparse


def convert_au_to_relativistic(input_file, alpha=7.297353e-3):
    """Convert atomic units to relativistic units using alpha^2 factor."""
    conversion_factor = alpha**2

    print(f"Fine structure constant α = {alpha}")
    print(f"Conversion factor α² = {conversion_factor:.6e}")
    print("=" * 50)

    # Derive output filenames from input filename
    base, ext = os.path.splitext(os.path.basename(input_file))
    output_file = f"{base}_rel.txt"
    summary_file = f"{base}_conversion_summary.txt"

    # Check if file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found!")
        return

    try:
        with open(input_file, 'r') as infile:
            with open(output_file, 'w') as outfile:
                # Write header to output file - Si³⁺ format
                outfile.write("# =================================================================\n")
                outfile.write("# Silicon 3+ (Si³⁺) Energy Levels in Relativistic Units\n")
                outfile.write("# =================================================================\n")
                outfile.write(f"# Converted from atomic units using α² = {conversion_factor:.6e}\n")
                outfile.write(f"# Fine structure constant α = {alpha:.6e}\n")
                outfile.write("# Units: relativistic (c=1, m_e=1, ℏ=1)\n")
                outfile.write("# =================================================================\n")
                outfile.write("#\n")
                outfile.write("# n    j     Energy (au)           Energy (rel units)        Ratio\n")
                outfile.write("# --   ---   -----------           ------------------        -----\n")
                
                print("Converting Si³⁺ energies:")
                print(f"{'n':<4} {'j':<6} {'Energy (au)':<20} {'Energy (rel)':<20} {'Ratio':<15}")
                print("-" * 75)
                
                line_count = 0
                converted_count = 0
                
                for line_num, line in enumerate(infile, 1):
                    line = line.strip()
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        outfile.write(line + "\n")
                        continue
                    
                    line_count += 1
                    
                    try:
                        # Parse the line: n j energy format
                        parts = line.split()
                        if len(parts) >= 3:
                            n = parts[0]
                            j = parts[1]
                            energy_au = float(parts[2])
                            
                            # Convert to relativistic units
                            energy_rel = energy_au * conversion_factor
                            ratio = conversion_factor
                            
                            # Format based on j value for quantum state naming
                            if j == '-1':
                                state = f"{n}s"
                            elif j == '1':
                                state = f"{n}p1/2"
                            elif j == '-2':
                                state = f"{n}p3/2"
                            elif j == '2':
                                state = f"{n}d3/2"
                            elif j == '-3':
                                state = f"{n}d5/2"
                            elif j == '3':
                                state = f"{n}f5/2"
                            elif j == '-4':
                                state = f"{n}f7/2"
                            else:
                                state = f"{n}(j={j})"
                            
                            # Write to output file in C_li_energies format
                            outfile.write(f"  {state:<8}  {energy_au:<20.6e}  {energy_rel:<20.6e}  {ratio:<15.6e}\n")
                            
                            # Print conversion info
                            print(f"{n:<4} {j:<6} {energy_au:<20.6e} {energy_rel:<20.6e} {ratio:<15.6e}")
                            
                            converted_count += 1
                        else:
                            # Skip lines that don't have the expected format
                            outfile.write(f"# {line}\n")
                            
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not process line {line_num}: "
                              f"{line}")
                        print(f"Error: {e}")
                        outfile.write(f"# Error processing: {line}\n")
                        continue
                
                # Add footer like C_li_energies format
                outfile.write("#\n")
                outfile.write("# =================================================================\n")
                outfile.write("# Notes:\n")
                outfile.write("# - These are binding energies (negative values)\n")
                outfile.write("# - j-j coupling scheme used for relativistic calculations\n")
                outfile.write("# - Si³⁺ has 11 electrons (Z=14, 3 electrons removed)\n")
                outfile.write("# - Fine structure splitting visible in p, d, and f orbitals\n")
                outfile.write("# =================================================================\n")
                
                print("-" * 75)
            print("Conversion completed!")
            print(f"Total lines processed: {line_count}")
            print(f"Energy values converted: {converted_count}")
            print(f"Output written to: {output_file}")

    except IOError as e:
        print(f"Error reading/writing files: {e}")
        return

    # Also create a summary file for quick reference
    with open(summary_file, 'w') as summary:
        summary.write("Energy Conversion Summary\n")
        summary.write("========================\n")
        summary.write(f"Input file: {input_file}\n")
        summary.write(f"Output file: {output_file}\n")
        summary.write(f"Conversion factor (α²): {conversion_factor:.12e}\n")
        summary.write(f"Fine structure constant (α): {alpha:.12e}\n")
        summary.write(f"Total energy values converted: {converted_count}\n")
        summary.write("\nConversion formula: E_relativistic = E_au × α²\n")
        summary.write("Units: atomic units → relativistic units "
                      "(c=1, mₑ=1, ℏ=1)\n")
    
    print(f"Summary written to: {summary_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Convert energies from atomic units to relativistic units')
    parser.add_argument('input_file', help='Input file containing energies in atomic units')
    parser.add_argument('--alpha', type=float, default=7.297353e-3,
                        help='Fine structure constant (default: 7.297353e-3)')
    args = parser.parse_args()
    convert_au_to_relativistic(args.input_file, args.alpha)
