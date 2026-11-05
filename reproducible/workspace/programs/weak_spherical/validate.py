#!/usr/bin/env python3
"""Record discretization and independent magnetic-response checks, not an error budget."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from spherical import Calculation, FermiDensity, Grid, NuclearInputs, HERE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    # Numerical examples preserve the old rms/diffuseness choices, but make the
    # equality of charge, proton and neutron profiles an explicit approximation.
    a = 2.3/(4*np.log(3))
    cases = [('Pb-208',82,126,5.5012), ('U-238',92,146,5.8571)]
    grids = [Grid(npoints=n) for n in (100, 160, 240, 320)]
    grids += [Grid(npoints=160, order=7), Grid(npoints=160, order=11),
              Grid(npoints=160, nq=40), Grid(npoints=160, cavity_bohr=5)]
    records = []
    for label,z,n,rms in cases:
        profile = FermiDensity.from_rms(rms, a)
        nucleus = NuclearInputs(z,n,profile,profile,profile,
                                'Common-profile numerical benchmark: charge=point-proton=neutron assumed; not a physical nuclear-density fit')
        for grid in grids:
            t = time.monotonic()
            c = Calculation(nucleus,grid)
            s,p = c.basis(-1),c.basis(1); i=s.state(1)
            magnetic=c.magnetic_response(s,i)
            row = {'ion':label, 'grid':asdict(grid), 'energy_1s_mc2':float(s.energies[i]),
                   'norm_1s':s.norm(i), 'delta_energy_pc_mc2':float(c.pc_column(s,i)[i]),
                   'pv_2p12_1s_imag_mc2':float(c.pv_column(p,s,i)[p.state(2)].imag),
                   'magnetic':magnetic, 'apv':c.apv_1s_2s(s,p)}
            if grid == Grid(npoints=160):
                row['mass_derivatives']={str(h):c.mass_derivative(h=h) for h in (1e-3,5e-4,2.5e-4)}
                row['mass_vs_direct_relative_difference']=abs(row['mass_derivatives']['0.00025']/magnetic['delta_g']-1)
            row['seconds']=time.monotonic()-t
            records.append(row)
            print(label, asdict(grid), 'dE=',row['delta_energy_pc_mc2'], 'dg=',magnetic['delta_g'],flush=True)
    checks = []
    for label,_,_,_ in cases:
        rows = [r for r in records if r['ion']==label]
        reference = next(r for r in rows if r['grid']['npoints']==320)
        # Compare the 160-point baseline and all refinements/variations to 320.
        # The 100-point result is retained as an additional coarse diagnostic.
        for row in rows:
            if row['grid']['npoints']==100:
                continue
            for field in ('delta_energy_pc_mc2','pv_2p12_1s_imag_mc2'):
                error=abs(row[field]/reference[field]-1)
                checks.append({'ion':label,'grid':row['grid'],'quantity':field,'relative_difference':error,
                               'tolerance':2e-6,'passed':error<2e-6})
            for field, tolerance, container in [('delta_g',2e-6,'magnetic'),('imag_e_bohr',2e-5,'apv')]:
                error=abs(row[container][field]/reference[container][field]-1)
                checks.append({'ion':label,'grid':row['grid'],'quantity':field,'relative_difference':error,
                               'tolerance':tolerance,'passed':error<tolerance})
            if 'mass_vs_direct_relative_difference' in row:
                error=row['mass_vs_direct_relative_difference']
                checks.append({'ion':label,'quantity':'mass derivative vs direct Zeeman','relative_difference':error,
                               'tolerance':3e-7,'passed':error<3e-7})
    result = {'status':'passed' if all(c['passed'] for c in checks) else 'failed',
              'scope':'Internal numerical validation of explicitly approximate spherical tree-level models; no precision APV benchmark',
              'build_manifest':json.loads((HERE/'build_manifest.json').read_text()),
              'python_source_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                      for name in ('spherical.py','validate.py','test_spherical.py')},
              'records':records,'checks':checks}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(result['status'],len(checks),'checks; saved',args.output)
    if result['status']!='passed':
        raise SystemExit(1)


if __name__=='__main__':
    main()
