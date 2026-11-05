"""Plot weak corrections from the validation record; see docs/spherical-model.md."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    data=json.loads(args.record.read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ion in ('Pb-208','U-238'):
        rows=[r for r in data['records'] if r['ion']==ion and r['grid']['order']==9
              and r['grid']['nq']==20 and r['grid']['cavity_bohr']==2.5]
        rows.sort(key=lambda r:r['grid']['npoints'])
        axes[0].plot([r['grid']['npoints'] for r in rows],
                     [r['delta_energy_pc_mc2']*1e12 for r in rows], 'o-',label=ion)
        reference=rows[-1]['magnetic']['delta_g']
        axes[1].semilogy([r['grid']['npoints'] for r in rows[:-1]],
                         [abs(r['magnetic']['delta_g']/reference-1) for r in rows[:-1]],'o-',label=ion)
    axes[0].set_ylabel('Parity-even energy shift / (10^-12 m_e c^2)')
    axes[1].set_ylabel('Relative delta-g difference from 320 points')
    for ax in axes:
        ax.set_xlabel('Radial breakpoint count')
        ax.legend()
        ax.grid(alpha=.25)
    fig.suptitle('One-electron common-profile model: internal numerical convergence')
    fig.savefig(args.output,dpi=180)
    plt.close(fig)

if __name__=='__main__':main()
