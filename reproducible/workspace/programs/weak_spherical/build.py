#!/usr/bin/env python3
"""Build only from repository source; never link a retained research object/library."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile

import numpy
import numpy.f2py

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CORE = ROOT / 'workspace/codes/bspqed/src'
SLATEC = ROOT / 'workspace/build/slatec/src'
NAME = 'weak_spherical_kernel'


def build():
    sources = [CORE / f'bspqed_{s}.f90' for s in
               ('function', 'physconsts', 'quadrature', 'uehling', 'potentials', 'dkb')]
    sources += [SLATEC / f'{s}.f' for s in
                ('dbvalu', 'dbspvd', 'dbspvn', 'dintrv', 'xermsg', 'xerprn', 'xersve',
                 'j4save', 'xerhlt', 'xgetua', 'i1mach', 'fdump', 'xercnt')]
    sources += [HERE / 'kernel.f90', ROOT / 'workspace/programs/rad/src/weak_angular_labels.f90']
    for path in sources:
        if not path.is_file():
            raise FileNotFoundError(path)
    commands = []
    with tempfile.TemporaryDirectory(prefix='build-', dir=HERE) as tmp:
        work = Path(tmp)
        def run(args):
            commands.append([str(x) for x in args])
            subprocess.run(args, cwd=work, check=True)
        objects = []
        flags = ['-O2', '-fPIC', '-fcheck=bounds', '-fbacktrace', '-ffree-line-length-none']
        for source in sources:
            obj = work / (source.stem + '.o')
            extra = ['-std=legacy'] if source.suffix == '.f' else []
            run(['gfortran', *flags, *extra, '-J', work, '-I', work, '-c', source, '-o', obj])
            objects.append(obj)
        # Resolve transitive SLATEC dependencies from source (including unused
        # Uehling callbacks), never from the copied libslatec or retained objects.
        while True:
            undefined, defined = set(), set()
            for obj in objects:
                for line in subprocess.check_output(['nm', '-g', obj], text=True).splitlines():
                    fields = line.split()
                    if len(fields) == 2 and fields[0] == 'U':
                        undefined.add(fields[1])
                    elif len(fields) >= 3:
                        defined.add(fields[-1])
            needed = sorted({SLATEC / (name.rstrip('_') + '.f')
                             for name in undefined - defined
                             if (SLATEC / (name.rstrip('_') + '.f')).is_file()} - set(sources))
            if not needed:
                break
            for source in needed:
                obj = work / (source.stem + '.o')
                run(['gfortran', *flags, '-std=legacy', '-c', source, '-o', obj])
                sources.append(source)
                objects.append(obj)
        # Compile the sampled driver against these same freshly built modules.
        # Its remaining angular/optical assumptions are labelled diagnostic.
        diagnostic = work / 'sampled_radial_diagnostic'
        driver = ROOT / 'workspace/programs/rad/src/something.f90'
        run(['gfortran', *flags, '-I', work, driver, *objects,
             '-llapack', '-lblas', '-o', diagnostic])
        sources.append(driver)
        # f2py reparses source arguments; use a relative filename so checkout
        # paths containing spaces remain valid. The original source is hashed.
        shutil.copy2(HERE / 'kernel.f90', work / 'kernel.f90')
        run([sys.executable, '-m', 'numpy.f2py', '-m', NAME, 'kernel.f90'])
        includes = [sysconfig.get_paths()['include'], numpy.get_include(), numpy.f2py.get_include()]
        for source in [work / f'{NAME}module.c', Path(numpy.f2py.get_include()) / 'fortranobject.c']:
            obj = work / (source.stem + '.o')
            run(['gcc', '-O2', '-fPIC', '-c', source, *['-I'+s for s in includes], '-o', obj])
            objects.append(obj)
        for source in sorted(work.glob('*-f2pywrappers*.f*')):
            obj = source.with_suffix('.o')
            run(['gfortran', *flags, '-I', work, '-c', source, '-o', obj])
            objects.append(obj)
        output = work / (NAME + sysconfig.get_config_var('EXT_SUFFIX'))
        # Only Fermi potentials are exposed; unused Uehling callbacks need no trampolines.
        run(['gfortran', '-shared', '-Wl,-z,noexecstack', '-o', output, *objects, '-llapack', '-lblas'])
        # Verify all runtime symbols resolve before publishing the new extension.
        run([sys.executable, '-c', f'import {NAME}; print({NAME}.spherical_constants())'])
        destination = HERE / output.name
        os.replace(output, destination)
        os.replace(diagnostic, HERE / diagnostic.name)
    manifest = {
        'source_complete': True,
        'python': sys.version, 'numpy': numpy.__version__,
        'compiler': subprocess.check_output(['gfortran', '--version'], text=True).splitlines()[0],
        'flags': flags,
        'sources_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'extension_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
        'commands': commands,
    }
    (HERE / 'build_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Built from source: {destination}')


if __name__ == '__main__':
    build()
