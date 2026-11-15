#!/usr/bin/env python3
"""Run all finite verification certificates: sage -python verify_all.py."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent


def main():
    if not __debug__ or sys.flags.optimize:
        raise RuntimeError('Assertions must be enabled; do not use -O.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=['ade', 'folded'])
    parser.add_argument('--details', action='store_true', help='Print computed mathematical certificates.')
    args = parser.parse_args()
    manifest = json.loads((ROOT/'SHA256SUMS.json').read_text(encoding='utf-8'))
    for rel, wanted in manifest.items():
        path = (ROOT/rel).resolve()
        if ROOT not in path.parents:
            raise ValueError('Manifest path outside the package: '+rel)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != wanted:
            raise RuntimeError('Source or expected-result checksum mismatch: '+rel)
    from sage.env import SAGE_VERSION
    import numpy
    import sympy
    print('SageMath', SAGE_VERSION, '| Python', sys.version.split()[0],
          '| NumPy', numpy.__version__, '| SymPy', sympy.__version__, flush=True)
    print('PASS: all supplied file checksums', flush=True)
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='brauer_verification_') as temporary:
        if args.only in (None, 'ade'):
            print('Reconstructing E6, E7, E8 ...', flush=True)
            report_path = Path(temporary)/'ade.json'
            result = subprocess.run([sys.executable, str(ROOT/'ade/run.py'), '--report', str(report_path)],
                                    cwd=temporary, env=env, text=True, capture_output=True)
            if result.returncode:
                print(result.stdout+result.stderr, file=sys.stderr)
                raise RuntimeError('ADE verification failed.')
            report = json.loads(report_path.read_text())
            if report['status'] != 'passed' or set(report['groups']) != {'e6','e7','e8'}:
                raise RuntimeError('Incomplete ADE report.')
            for name in report['groups']:
                print(name.upper()+': PASS', flush=True)
            if args.details:
                print(json.dumps(report['groups'], indent=2, sort_keys=True))
        if args.only in (None, 'folded'):
            print('Reconstructing G2, H3, F4, H4 ...', flush=True)
            cmd = [sys.executable, str(ROOT/'folded/verify_all.py')]
            if args.details:
                cmd.append('--details')
            result = subprocess.run(cmd, cwd=temporary, env=env, text=True, capture_output=True)
            if result.returncode:
                print(result.stdout+result.stderr, file=sys.stderr)
                raise RuntimeError('Folded verification failed.')
            report = json.loads(result.stdout)
            if set(report) != {'G2','H3','F4','H4'}:
                raise RuntimeError('Incomplete folded report.')
            for name in report:
                print(name+': PASS', flush=True)
            if args.details:
                print(json.dumps(report, indent=2, sort_keys=True))
    print('PASS: all selected exact calculations and expected results; %.3f seconds' %
          (time.monotonic()-started), flush=True)


if __name__ == '__main__':
    main()
