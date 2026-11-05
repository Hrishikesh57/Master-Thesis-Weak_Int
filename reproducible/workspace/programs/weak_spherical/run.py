#!/usr/bin/env python3
"""Run an explicitly specified spherical model, without importing legacy binaries."""
import argparse
import hashlib
import json
from pathlib import Path

from spherical import Calculation, FermiDensity, Grid, NuclearInputs, HERE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True, help='JSON with independent charge, proton and neutron c_fm/a_fm.')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--npoints', type=int, default=200)
    parser.add_argument('--order', type=int, default=9)
    parser.add_argument('--quadrature', type=int, choices=[20, 40], default=20)
    parser.add_argument('--cavity-bohr', type=float, default=2.5)
    args = parser.parse_args()
    raw = json.loads(args.inputs.read_text())
    nucleus = NuclearInputs(raw['z'], raw['n'], *(FermiDensity(**raw[k]) for k in ('charge', 'proton', 'neutron')), raw['provenance'])
    calc = Calculation(nucleus, Grid(args.npoints, args.order, args.quadrature, args.cavity_bohr))
    s = calc.basis(-1); i = s.state(1)
    p = calc.basis(1)
    magnetic = calc.magnetic_response(s, i)
    derivatives = {str(h): calc.mass_derivative(h=h) for h in (1e-3, 5e-4)}
    result = {'status': 'spherical one-electron model; no nuclear deformation or precision uncertainty budget',
              'inputs': calc.inputs(), 'input_sha256': hashlib.sha256(args.inputs.read_bytes()).hexdigest(),
              'build_manifest': json.loads((HERE/'build_manifest.json').read_text()),
              'python_source_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                       for name in ('spherical.py', 'run.py')},
              'density_rms_fm': {k: getattr(nucleus, k).rms_fm for k in ('charge', 'proton', 'neutron')},
              'qweak': calc.cvp*nucleus.z-nucleus.n,
              'energy_1s_mc2': float(s.energies[i]), 'normalization_1s': s.norm(i),
              'delta_energy_pc_mc2': float(calc.pc_column(s, i)[i]),
              'pv_2p12_1s_imag_mc2': float(calc.pv_column(p, s, i)[p.state(2)].imag),
              'magnetic_response': magnetic, 'delta_g_mass_derivative': derivatives,
              'magnetic_relative_disagreement': abs(derivatives['0.0005']/magnetic['delta_g']-1),
              'apv_1s_2s': calc.apv_1s_2s(s, p)}
    # Do not overwrite a previous scientific record silently.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('delta_energy_pc_mc2', 'pv_2p12_1s_imag_mc2', 'magnetic_response',
                                           'magnetic_relative_disagreement', 'apv_1s_2s')}, indent=2))
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
