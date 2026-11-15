"""Exact W(E7) restriction of pair layers induced from N(X).
Uses N=((B2 x B4) cap D6) x A1 and quotient S4 x C2.
"""
import e7_singleton as s
import numpy as np,itertools,json
from collections import Counter
parts=list(s.parts(4)); labels=[(p,e) for p in parts for e in [1,-1]]
acc=np.zeros((60,10),dtype=np.int64)
for a in itertools.permutations(range(2)):
 for b in itertools.permutations(range(2,6)):
  perm=list(a)+list(b)+[6,7]
  # cycle partition of underlying permutation on last4
  seen=set();cy=[]
  for i in range(2,6):
   if i in seen:continue
   j=i;k=0
   while j not in seen:seen.add(j);k+=1;j=perm[j]
   cy.append(k)
  cy=tuple(sorted(cy,reverse=True))
  for signs in itertools.product([-1,1],repeat=6):
   if np.prod(signs)!=1:continue
   sg=list(signs)+[1,1]
   roots=np.array(s.roots);out=np.empty_like(roots)
   for i in range(8):out[:,perm[i]]=sg[i]*roots[:,i]
   for flip in [0,1]:
    if flip:out[:,[6,7]]=out[:,[7,6]]
    p=np.array([s.lookup[tuple(r)] for r in out]);f=s.fusion(p)
    theta=np.array([s.symchar(l,cy)*(eps if flip else 1) for l,eps in labels])
    acc+=np.outer(s.chars[:,f],theta)
assert np.all(acc%3072==0)
b=acc//3072
assert np.all(b>=0)
dim=np.array([s.symchar(l,(1,)*4) for l,e in labels])
assert np.array_equal(s.chars[:,0]@b,945*dim)
cs=[2*sum(j-i for i,x in enumerate(l) for j in range(x))+e for l,e in labels]
extras=set(); candidates=set(); rows=[]
for j,(l,e) in enumerate(labels):
 vals={int(1+(cs[j]-s.cw[i])//2) for i in range(60) if b[i,j] and (cs[j]-s.cw[i])%2==0}-{0}
 candidates|=vals;extras|=vals-s.obs
 rows.append(dict(partition=l,epsilon=e,c=cs[j],s0_candidates=sorted(vals),extra=sorted(vals-s.obs),branching=b[:,j].tolist()))
json.dump(dict(labels=rows,all_s0_candidates=sorted(candidates),extra=sorted(extras)),open('e7_pair_branching.json','w'),indent=2)
print('Pair s0 extras:',sorted(extras));print('s0 candidate union:',sorted(candidates))
