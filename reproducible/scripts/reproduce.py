"""Run the source-built one-electron example; see docs/spherical-model.md."""
import argparse
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
PROGRAM=ROOT/'workspace/programs/weak_spherical'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/demo'))
    parser.add_argument('--build',action='store_true',help='rebuild the extension from source first')
    parser.add_argument('--validate',action='store_true',help='also repeat the 16-case numerical convergence study')
    args=parser.parse_args()
    output=args.output.resolve()
    if output.exists():
        parser.error('Output directory already exists; choose a fresh path.')
    if args.build:
        subprocess.run([sys.executable,str(PROGRAM/'build.py')],cwd=ROOT,check=True)
    output.mkdir(parents=True)
    subprocess.run([sys.executable,str(PROGRAM/'run.py'),'--inputs',str(PROGRAM/'example_pb208.json'),
                    '--output',str(output/'pb208.json')],cwd=ROOT,check=True)
    if args.validate:
        subprocess.run([sys.executable,str(PROGRAM/'validate.py'),'--output',str(output/'validation.json')],cwd=ROOT,check=True)
        subprocess.run([sys.executable,str(ROOT/'scripts/plot_validation.py'),str(output/'validation.json'),
                        str(output/'convergence.png')],cwd=ROOT,check=True)
    print('Saved',output)

if __name__=='__main__':main()
