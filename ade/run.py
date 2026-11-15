#!/usr/bin/env python3
"""Run ADE certificates in disposable scratch space: sage -python ade/run.py."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent
STEPS={
    'e6':['e6_spectrum.py'],
    'e7':['e7_fusion_check.py','e7_pair_branching.py','e7_higher_branching.py',
          'e7_pair_sign.py','e7_pair_standard_sage.py','e7_four_trivial.py'],
    'e8':['e8_singleton.py','e8_pair_branching.py','e8_four_branching.py',
          'e8_final_screen.py','e8_exact_inputs.py','e8_relations.py','e8_annihilator.py'],
}

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def read(directory,name):
    return json.loads((directory/name).read_text())

def summarize(group,d):
    if group=='e6':
        a=read(d,'e6_spectrum_recomputed.json')
        return {'singleton_bad_set':[r['delta'] for r in a['spectrum']],
                'singleton_factorizations':a['singleton_factors'],
                'singleton_branching_sha256':digest(a['singleton_branching']),
                'symbolic_monomial_relation_checks':a['symbolic_monomial_relation_checks'],
                'higher_determinants':a['higher_determinant_degree_exhaustion'],
                'standard_pair_at_one':a['standard_pair_at_one']}
    if group=='e7':
        a=read(d,'e7_singleton_results.json')
        ranks=[]
        for name in ['pair_sign','pair_standard','four_trivial']:
            r=read(d,'e7_'+name+'_result.json')
            assert r['rank_mod_prime']==r['dimension'],(name,r)
            ranks.append({'module':name,**{k:r[k] for k in ['dimension','delta','prime','rank_mod_prime']}})
        return {'singleton_bad_set':a['exceptional_set'],
                'singleton_full_data_sha256':digest(a),
                'pair_branching_sha256':digest(read(d,'e7_pair_branching.json')),
                'higher_branching_sha256':digest(read(d,'e7_higher_branching.json')),
                'residual_ranks':ranks}
    a=read(d,'e8_singleton_results.json')
    r=read(d,'e8_relations.json')
    c=read(d,'annihilator_verification.json')
    c=c.get('E8_residual_certificates',c.get('E8_complete_criterion'))
    certs=c['six_residual_channel_certificates']
    assert len(certs)==4 and sum(len(q['W_E8_labels']) for q in certs)==6
    assert all(q['value_at_delta_target']!=0 for q in certs)
    return {'singleton_bad_set':a['singleton_bad_set'],
            'singleton_fusion_branching_sha256':digest({k:a[k] for k in ['fusion_identity','fusion_reflection','branching','determinant_factors']}),
            'pair_branching_sha256':digest(read(d,'e8_pair_branching_results.json')),
            'four_branching_sha256':digest(read(d,'e8_four_branching_results.json')),
            'final_screen':read(d,'e8_final_screen_results.json'),
            'formal_relations':{k:v for k,v in r.items() if k!='note'},
            'projector_certificates_sha256':digest(certs),
            'residual_projectors':[{k:q[k] for k in ['delta_target','S6_partition','W_E8_labels','Z1_eigenvalue_c0','a_coeff','b_coeff','value_at_delta_target']} for q in certs],
            'multiplicity_one_quadratic_count':c['total_multiplicity_one_constituents']}

def main():
    if not __debug__ or sys.flags.optimize:
        raise RuntimeError('Assertions must remain enabled; do not use -O or PYTHONOPTIMIZE.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',choices=list(STEPS),help='Run one type; default runs all.')
    parser.add_argument('--report',type=Path,help='Optional actual result JSON (never an expected file).')
    args=parser.parse_args()
    if args.report and (ROOT/'expected').resolve() in [args.report.resolve(),*args.report.resolve().parents]:
        raise ValueError('Refusing to overwrite supplied expected results.')
    provenance=read(ROOT,'PROVENANCE.json')
    for record in provenance['files']:
        assert hashlib.sha256((ROOT/record['file']).read_bytes()).hexdigest()==record['distributed_sha256'],record['file']
    expected=read(ROOT/'expected','results.json')
    groups=[args.only] if args.only else list(STEPS)
    env=os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONOPTIMIZE',None)
    report={'status':'passed','groups':{},'runs':[]}
    start=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='brauer_ade_') as tmp:
        scratch=Path(tmp)/'code'
        shutil.copytree(ROOT/'code',scratch)
        for group in groups:
            for script in STEPS[group]:
                print('RUN',group,script,flush=True)
                t=time.monotonic()
                result=subprocess.run([sys.executable,str(scratch/script)],cwd=scratch,env=env,
                                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                if result.returncode:
                    print(result.stdout,flush=True)
                    raise RuntimeError(f'{script} failed with exit code {result.returncode}')
                elapsed=round(time.monotonic()-t,3)
                report['runs'].append({'script':script,'seconds':elapsed})
                print('PASS',script,elapsed,'seconds',flush=True)
            actual=summarize(group,scratch)
            if actual!=expected[group]:
                print(json.dumps({'actual':actual,'expected':expected[group]},indent=2))
                raise AssertionError(group+' differs from the supplied expected certificate')
            report['groups'][group]=actual
            print('PASS',group,'all exact assertions and expected results',flush=True)
    report['seconds']=round(time.monotonic()-start,3)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS ADE certificate run:',', '.join(groups),report['seconds'],'seconds',flush=True)

if __name__=='__main__':
    main()
