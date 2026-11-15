"""Exact G2 models: nu_1=1, nu_2=3; triple pair (i,j)=(2,1)."""
import json
import sympy as s

def verify():
    z=s.Symbol('delta')
    G=s.Matrix([[z,1,1],[1,z,1],[1,1,z]])
    models=[]
    for eps in [1,-1]:
        models.append((f'middle_{eps}',
          [s.Matrix([[1,0,0],[0,0,1],[0,1,0]]),eps*s.Matrix([[0,1,0],[1,0,0],[0,0,1]])],
          [s.Matrix([[z,1,1],[0,0,0],[0,0,0]]),s.zeros(3)],G))
    models.append(('minimal',
        [s.Matrix([[0,1,0],[1,0,0],[0,0,1]]),s.Matrix([[1,0,0],[0,0,1],[0,1,0]])],
        [s.Matrix([[0,0,0],[0,0,0],[1,1,z]]),s.Matrix([[z**3,z**2,z**2],[0,0,0],[0,0,0]])],z**2*G))
    results=[]
    for name,R,E,H in models:
        checks=0
        def eq(A,B,label):
            nonlocal checks
            assert all(s.expand(v)==0 for v in A-B),(name,label,A-B)
            checks+=1
        for a,nu in [(0,1),(1,3)]:
            eq(R[a]**2,s.eye(3),'r^2')
            eq(R[a]*E[a],E[a],'re')
            eq(E[a]*R[a],E[a],'er')
            eq(E[a]**2,z**nu*E[a],'e^2')
        i,j=1,0
        eq((R[j]*R[i])**6,s.eye(3),'triple Coxeter')
        eq(R[i]*E[j]*E[i],R[j]*E[i],'triple 38')
        eq(E[i]*E[j]*R[i],E[i]*R[j],'triple 39')
        eq(E[j]*R[i]*E[j]*R[i]*E[j],E[j],'triple 40')
        eq(E[j]*R[i]*E[j]*R[i]*R[j],E[j]*R[i]*R[j]*R[i],'triple 41')
        eq(E[i]*R[j]*E[i],z**2*E[i],'triple 42')
        eq(R[j]*R[i]*E[j]*R[i]*E[j],R[i]*R[j]*R[i]*E[j],'triple 43')
        relation_checks=checks
        for A in R+E:
            eq(A.T*H,H*A,'contravariance')
        det=s.factor(H.det())
        expected=z**6*(z-1)**2*(z+2) if name=='minimal' else (z-1)**2*(z+2)
        assert s.expand(det-expected)==0
        results.append({'model':name,'relation_checks':relation_checks,
                        'contravariance_checks':checks-relation_checks,'determinant':str(det)})
    return {'status':'PASS','models':results,
            'symbolic_relations':sum(v['relation_checks'] for v in results),
            'contravariance_identities':sum(v['contravariance_checks'] for v in results)}

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
