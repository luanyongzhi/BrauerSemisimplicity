"""Reconstruct all 14 F4 forms and certify their polynomial determinants.

Each presentation/contravariance identity has entry degree at most six and
is checked at seven integers. Each determinant is checked over ZZ at N*d+1
integers, where N is its order and d its entry-degree bound.
"""
import json
import sympy as s
import numpy as np
from sage.all import matrix, ZZ
import f4_model as model

def verify():
    D=model.roots_data()
    cells=model.make_cells(D)
    relations={c['label']:model.verify_relations(c) for c in cells}
    z=s.Symbol('z')
    expected={
     '2:triv': z**27*(z-2)**15*(z+1)**9*(z+2)**2*(z*z+5*z-8)**2*(z**3+7*z*z+14*z+32),
     '2:sign': z**27*(z-4)**2*(z-3)**9*(z+1)**2*(z+2)**18*(z*z-3*z+8),
     '2:std': 3**30*z**54*(z-3)**9*(z-1)**22*(z+1)**9*(z+4)**12*(z*z-z-8)**4*(z**3-2*z*z+5*z-16)**2,
     '4:triv': z**72*(z-2)**15*(z-1)**38*(z+2)**28*(z+4)**3,
    }
    f0=z**9*(z-4)**2*(z+8)
    f1=(z-2)**8*(z+4)**4
    f2=z**6*(z-4)**3*(z+4)**3
    f3=(z-4)**4*(z+2)**8
    for name,k,kind in model.SINGLE:
        f=[f0,f1,f2,f3][k]
        if k in [1,2]: value=f**3
        elif kind=='std': value=3**12*f**2
        else: value=f
        expected['1:'+name]=value
    certificates=[]
    for c in cells:
        n=c['dim']; degree=n*c['layer']; coeff=c['G']
        assert s.degree(expected[c['label']],z)==degree
        poly=s.Poly(expected[c['label']],z)
        # Use Python/Sage integers at these large evaluation points, avoiding int64.
        for a in range(degree+1):
            entries=[sum(int(coeff[k,i,j])*a**k for k in range(len(coeff)))
                     for i in range(n) for j in range(n)]
            actual=matrix(ZZ,n,n,entries).det()
            assert actual==int(poly.eval(a)),(c['label'],a)
        certificates.append({'label':c['label'],'order':n,'entry_degree_bound':c['layer'],
                             'determinant_degree_bound':degree,'integer_points':degree+1,
                             'evaluation_range':[0,degree],'factorization':str(s.factor(expected[c['label']]))})
    # At the small relation-check points |z|<=4, all contraction entries have
    # absolute value <=32, group entries <=2, Gram entries <=512 and N<=60.
    # The largest relation products use three contractions, or two contractions
    # and two group matrices. The explicit upper bounds are far below int64 max.
    max_e=max(abs(int(v)) for c in cells for ec in c['Es'] for a in [-2,-1,0,1,2,3,4]
              for v in sum((a**k*A for k,A in enumerate(ec)),np.zeros_like(ec[0])).flat)
    max_r=max(abs(int(v)) for c in cells for R in c['Rs'] for v in R.flat)
    max_g=max(abs(int(v)) for c in cells for a in [-2,-1,0,1,2,3,4]
              for v in sum((a**k*A for k,A in enumerate(c['G'])),np.zeros_like(c['G'][0])).flat)
    N=max(c['dim'] for c in cells)
    bound=max(N**2*max_e**3,N**3*max_e**2*max_r**2,N*max_e*max_g,N**3*max_r**4,16*max_e)
    assert bound<2**63
    return {'status':'PASS','positive_E6_roots':len(D['R']),
            'nongroup_cells':len(cells),'polynomial_identities':sum(relations.values())//7,
            'relation_evaluations':sum(relations.values()),
            'determinant_evaluations':sum(c['integer_points'] for c in certificates),
            'int64_bound':bound,'determinants':certificates}

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
