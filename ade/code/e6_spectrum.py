#!/usr/bin/env python3
"""Reconstruct E6 singleton branching and all higher-layer ranks exactly.

Only CHEVIE's pinned weyle6.jl is mathematical input.  Higher matrices are
rebuilt from the CFW root-set action.  Sage rational rank is used directly,
independently of the old modular-rank/kernel-reconstruction implementation.
No historical result file is used as mathematical input.
Run with sage -python, without -O.
"""
from pathlib import Path
import ast
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
import math
import platform
import re
import time
import numpy as np
from sage.all import matrix, QQ
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parent
S6 = [-20,-14,-12,-8,-6,-5,-4,-3,-2,1,2,4,5,6,7,8,10]

def partitions(n, ceiling=None):
    if n == 0:
        yield ()
        return
    for a in range(min(n, n if ceiling is None else ceiling), 0, -1):
        for rest in partitions(n-a, a):
            yield (a,)+rest

def cells(lam):
    return {(i,j) for i,row in enumerate(lam) for j in range(row)}

@lru_cache(None)
def schar(lam, cycles):
    if not cycles:
        return int(not lam)
    k, *rest = cycles
    full = cells(lam)
    ans = 0
    for mu in partitions(sum(lam)-k):
        inner = cells(mu)
        if not inner <= full:
            continue
        skew = full-inner
        if any({(i,j),(i+1,j),(i,j+1),(i+1,j+1)} <= skew for i,j in skew):
            continue
        reached = {next(iter(skew))}
        frontier = list(reached)
        while frontier:
            i,j = frontier.pop()
            for q in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
                if q in skew and q not in reached:
                    reached.add(q); frontier.append(q)
        if reached == skew:
            ans += (-1)**(len({i for i,j in skew})-1)*schar(mu,tuple(rest))
    return ans

def hook_dimension(lam):
    h = math.prod(lam[i]-j+sum(row>j for row in lam[i+1:])
                  for i,row in enumerate(lam) for j in range(row))
    return math.factorial(sum(lam))//h

def singleton_spectrum():
    src=(ROOT/'weyle6.jl').read_text()
    def extract(marker):
        tail=src.split(marker,1)[1]
        start=tail.index('['); depth=0
        for j in range(start,len(tail)):
            depth += (tail[j]=='[')-(tail[j]==']')
            if depth==0:return tail[start:j+1]
        raise ValueError(marker)
    chars=extract('Vector{Pol{Int64}}')
    chars=re.sub(r'Pol\(\s*(\[[^\]]*\])(?:\s*,\s*-?\d+)?\s*\)',
                 lambda m:str(sum(ast.literal_eval(m[1]))),chars)
    ch=np.array(ast.literal_eval(chars.replace('Pol()','1')),dtype=np.int64)
    sizes=np.array(ast.literal_eval(extract(':classes=>')),dtype=np.int64)
    words=ast.literal_eval(extract(':classwords =>').replace('Int[]','[]'))
    labels=ast.literal_eval(extract(':charparams,'))
    assert np.array_equal((ch*sizes)@ch.T,51840*np.eye(25,dtype=np.int64))
    assert list(ch[:,0])==[p[0] for p in labels]
    simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],
                     [-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],
                     [0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0]],dtype=np.int64)
    roots=set(map(tuple,simple)); todo=list(roots)
    while todo:
        r=np.array(todo.pop())
        for a in simple:
            q=tuple(r-(r@a)*a//4)
            if q not in roots:roots.add(q);todo.append(q)
    roots=np.array(sorted(roots)); lookup={tuple(r):i for i,r in enumerate(roots)}
    assert len(roots)==72
    def reflection(a):return np.array([lookup[tuple(r-(r@a)*a//4)] for r in roots])
    gens=list(map(reflection,simple)); identity=np.arange(72)
    def word(w,gg=gens):
        p=identity.copy()
        for i in w:p=gg[i-1][p]
        return p
    def cyc(p):
        seen=set(); c=Counter()
        for i in range(len(p)):
            if i in seen:continue
            j=i;k=0
            while j not in seen:seen.add(j);j=int(p[j]);k+=1
            if k>1:c[k]+=1
        return tuple(sorted(c.items()))
    classcycles=[cyc(word(w)) for w in words]
    assert len(set(classcycles))==25
    fusion={x:i for i,x in enumerate(classcycles)}
    a0=simple[0];s0=reflection(a0)
    pos={tuple(r) for r in roots if r@a0==0 and next(x for x in r if x)<0}
    assert len(pos)==15
    ss=[np.array(r) for r in pos if not any(tuple(np.array(r)-np.array(q)) in pos for q in pos)]
    assert len(ss)==5
    adj={i:[j for j in range(5) if ss[i]@ss[j]==-4] for i in range(5)}
    order=[next(i for i,v in adj.items() if len(v)==1)]
    while len(order)<5:order.append(next(i for i in adj[order[-1]] if i not in order))
    gg=list(map(reflection,[ss[i] for i in order]))
    assert all(np.array_equal(s0[g],g[s0]) for g in gg)
    parts=list(partitions(6)); fusion0=[];fusion1=[]
    for p in parts:
        w=[];offset=0
        for k in p:w+=list(range(offset+1,offset+k));offset+=k
        h=word(w,gg);fusion0.append(fusion[cyc(h)]);fusion1.append(fusion[cyc(s0[h])])
    symchars=np.array([[schar(l,p) for p in parts] for l in parts],dtype=np.int64)
    symsizes=np.array([720//math.prod(k**v*math.factorial(v) for k,v in Counter(p).items()) for p in parts])
    assert np.array_equal((symchars*symsizes)@symchars.T,720*np.eye(11,dtype=np.int64))
    numerator=((ch[:,fusion0]+ch[:,fusion1])*symsizes)@symchars.T
    assert np.all(numerator%1440==0)
    branch=numerator//1440
    assert branch.min()>=0
    dimensions=np.array([36*hook_dimension(p) for p in parts],dtype=np.int64)
    assert np.array_equal(ch[:,0]@branch,dimensions)
    reflection_col=fusion[cyc(s0)]
    c=36*ch[:,reflection_col]//ch[:,0]
    assert np.array_equal(c*ch[:,0],36*ch[:,reflection_col])
    content=[sum(j-i for i,row in enumerate(p) for j in range(row)) for p in parts]
    factors=[]
    for j in range(11):
        f=Counter()
        for i in range(25):
            if branch[i,j]:f[int(1+content[j]-c[i])]+=int(branch[i,j]*ch[i,0])
        assert sum(f.values())==dimensions[j]
        factors.append(dict(sorted(f.items())))
    assert sorted(set().union(*(set(f) for f in factors))-{0})==S6
    return parts,dimensions,factors,branch.tolist(),labels

def higher_matrices():
    C=2*np.eye(6,dtype=np.int64)
    for a,b in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[a,b]=C[b,a]=-1
    roots=set(map(tuple,np.eye(6,dtype=np.int64))); todo=list(roots)
    while todo:
        r=np.array(todo.pop())
        for i in range(6):
            q=r.copy();q[i]-=(C@r)[i];q=tuple(q)
            if q not in roots:roots.add(q);todo.append(q)
    roots=sorted(r for r in roots if all(x>=0 for x in r))
    R=np.array(roots);G=R@C@R.T;lookup={r:i for i,r in enumerate(roots)}
    def line(r):
        if next(x for x in r if x)<0:r=-r
        return lookup[tuple(r)]
    reflections=np.array([[line(R[b]-G[a,b]*R[a]) for b in range(36)] for a in range(36)])
    def reflect(a,X):return tuple(sorted(int(reflections[a,x]) for x in X))
    pairs=[x for x in itertools.combinations(range(36),2) if G[x]==0]
    fours=[x for x in itertools.combinations(range(36),4) if all(G[a,b]==0 for a,b in itertools.combinations(x,2))]
    assert (len(roots),len(pairs),len(fours))==(36,270,135)
    idx={X:i for i,X in enumerate(pairs)}
    decor={X:tuple(Y for Y in pairs if all(G[a,b]==0 for a in X for b in Y)) for X in pairs}
    assert {len(v) for v in decor.values()}=={3}
    def check_relations(sets, decorations):
        """Compare polynomial monomial operators, before any specialization."""
        index={X:i for i,X in enumerate(sets)}
        width=len(decorations[sets[0]])
        n=len(sets)*width
        identity=(np.arange(n),np.zeros(n,dtype=np.int64))
        def operator(a, is_e):
            dest=[];powers=[]
            for X in sets:
                for D in decorations[X]:
                    if not is_e:
                        Y=reflect(a,X);F=reflect(a,D);k=0
                    elif a in X:
                        Y=X;F=D;k=1
                    else:
                        bs=[b for b in X if G[a,b]]
                        if not bs:dest.append(-1);powers.append(0);continue
                        b=bs[0];Y=reflect(b,reflect(a,X));F=reflect(b,reflect(a,D));k=0
                    dest.append(width*index[Y]+decorations[Y].index(F));powers.append(k)
            return np.array(dest),np.array(powers)
        def product(*args):
            result=identity
            for A in reversed(args):
                d=np.where(result[0]<0,-1,A[0][result[0]])
                powers=np.where(d<0,0,result[1]+A[1][result[0]])
                result=(d,powers)
            return result
        count=0
        def equal(A,B):
            nonlocal count
            assert np.array_equal(A[0],B[0]) and np.array_equal(A[1],B[1])
            count+=1
        simple=[lookup[tuple(x)] for x in np.eye(6,dtype=np.int64)]
        R=[operator(a,False) for a in simple]
        E=[operator(a,True) for a in simple]
        for i in range(6):
            equal(product(R[i],R[i]),identity)
            equal(product(R[i],E[i]),E[i]);equal(product(E[i],R[i]),E[i])
            equal(product(E[i],E[i]),(E[i][0],np.where(E[i][0]<0,0,E[i][1]+1)))
            for j in range(6):
                if i==j:continue
                if C[i,j]==0:
                    for A,B in [(R[i],R[j]),(R[i],E[j]),(E[i],E[j])]:equal(product(A,B),product(B,A))
                else:
                    equal(product(R[i],R[j],R[i]),product(R[j],R[i],R[j]))
                    equal(product(R[j],E[i],R[j]),product(R[i],E[j],R[i]))
                    equal(product(R[j],R[i],E[j]),product(E[i],E[j]))
                    equal(product(E[i],E[j],E[i]),E[i])
            for a in range(36):
                equal(product(R[i],operator(a,True),R[i]),operator(int(reflections[simple[i],a]),True))
        return count
    relation_counts={
        'pair_decoration_dimension_810':check_relations(pairs,decor),
        'four_dimension_135':check_relations(fours,{X:((),) for X in fours}),
    }
    matrices={key:[np.zeros((n,n),dtype=np.int64) for _ in range(3)] for key,n in [('pair_trivial',270),('pair_sign',270),('pair_standard',540)]}
    for j,X in enumerate(pairs):
        for i,Y in enumerate(pairs):
            dest=X;D=decor[X];k=0
            for a in reversed(Y):
                if a in dest:k+=1;continue
                bs=[b for b in dest if G[a,b]]
                if not bs:dest=None;break
                b=bs[0]
                dest=reflect(b,reflect(a,dest))
                D=tuple(reflect(b,reflect(a,Z)) for Z in D)
            if dest is None:continue
            assert dest==Y
            perm=[decor[Y].index(z) for z in D]
            sign=(-1)**sum(perm[a]>perm[b] for a in range(3) for b in range(a+1,3))
            P=np.eye(3,dtype=np.int64)[:,perm]
            matrices['pair_trivial'][k][i,j]=1
            matrices['pair_sign'][k][i,j]=sign
            matrices['pair_standard'][k][2*i:2*i+2,2*j:2*j+2]=(P[:,:2]-P[:,2:3])[:2,:]
    F=[np.zeros((135,135),dtype=np.int64) for _ in range(5)]
    for j,X in enumerate(fours):
        for i,Y in enumerate(fours):
            dest=X;k=0
            for a in reversed(Y):
                if a in dest:k+=1;continue
                bs=[b for b in dest if G[a,b]]
                if not bs:dest=None;break
                dest=reflect(bs[0],reflect(a,dest))
            if dest is not None:
                assert dest==Y
                F[k][i,j]=1
    matrices['four_trivial']=F
    for key,A in matrices.items():
        assert np.array_equal(A[-1],np.eye(A[0].shape[0],dtype=np.int64))
    return matrices,relation_counts

def exact_rank(A):
    return int(matrix(QQ,A.tolist()).rank())

def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    start=time.monotonic()
    parts,dims,factors,branch,labels=singleton_spectrum()
    print('PASS reconstructed singleton characters, branching and 11 spectra',flush=True)
    matrices,relation_counts=higher_matrices()
    print('Built four higher CFW polynomial matrices',flush=True)
    result={'method':'fresh root-set CFW matrices and direct Sage rational ranks; singleton CHEVIE character reconstruction',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__,'sage':SAGE_VERSION},
            'character_source_sha256':hashlib.sha256((ROOT/'weyle6.jl').read_bytes()).hexdigest(),
            'symbolic_monomial_relation_checks':relation_counts,
            'singleton_partitions':[list(p) for p in parts],
            'singleton_generic_dimensions':list(map(int,dims)),
            'singleton_factors':factors,'singleton_branching':branch,'weyl_labels':labels,
            'spectrum':[]}
    for d in S6:
        sr=[int(n)-f.get(d,0) for n,f in zip(dims,factors)]
        ranks={key:exact_rank(sum(A*d**k for k,A in enumerate(coeffs))) for key,coeffs in matrices.items()}
        quotient=51840+sum(x*x for x in sr)+sum(x*x for x in ranks.values())
        degenerate=sum(x!=n for x,n in zip(sr,dims))+sum(ranks[key]!=As[0].shape[0] for key,As in matrices.items())
        row={'delta':d,'singleton_ranks':sr,'singleton_square_sum':sum(x*x for x in sr),'higher_ranks':ranks,
             'semisimple_quotient_dimension':quotient,'radical_dimension':1440585-quotient,
             'degenerate_nonempty_cells':int(degenerate),'number_of_simple_modules':25+sum(x>0 for x in sr)+sum(x>0 for x in ranks.values())}
        result['spectrum'].append(row)
        print('PASS',d,row['radical_dimension'],ranks,flush=True)
        (ROOT/'e6_spectrum_recomputed.json').write_text(json.dumps(result,indent=2)+'\n')
    A=matrices['pair_standard'];T=sum(A);Tp=A[1]+2*A[2]
    jet=np.block([[T,np.zeros_like(T)],[Tp,T]])
    result['standard_pair_at_one']={'matrix_nullity':540-exact_rank(T),'order_two_jet_nullity':1080-exact_rank(jet)}
    assert result['standard_pair_at_one']=={'matrix_nullity':406,'order_two_jet_nullity':466}
    # For a monic polynomial matrix, a nullity n at a gives at least n
    # powers of (delta-a) in its determinant.  Jet nullity gives the stronger
    # sum min(2,nu_i).  Exhaustion of the total degree makes all bounds exact.
    determinant_checks={}
    for key,coeffs in matrices.items():
        order=coeffs[0].shape[0]
        exponents={row['delta']:order-row['higher_ranks'][key] for row in result['spectrum']
                   if order-row['higher_ranks'][key]}
        at_zero=order-exact_rank(coeffs[0])
        if at_zero:exponents[0]=at_zero
        if key=='pair_standard':exponents[1]=result['standard_pair_at_one']['order_two_jet_nullity']
        degree=order*(len(coeffs)-1)
        assert sum(exponents.values())==degree,(key,exponents,degree)
        determinant_checks[key]={'order':order,'degree':degree,
                                 'exact_linear_factor_exponents':dict(sorted(exponents.items())),
                                 'degree_exhausted':True}
    result['higher_determinant_degree_exhaustion']=determinant_checks
    result['elapsed_seconds']=round(time.monotonic()-start,3)
    assert all(row['number_of_simple_modules']==40 for row in result['spectrum'])
    (ROOT/'e6_spectrum_recomputed.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS complete exact E6 spectrum and jet; seconds',result['elapsed_seconds'],flush=True)

if __name__=='__main__':main()
