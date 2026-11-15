#!/usr/bin/env python3
"""Exact independent H4 spectra from relation-certified group-ring data.
Requirements: Python3 and NumPy.  No SymPy, rank guesses, floating point,
or Round12 coefficient files are used.  All int64 operations are preceded
by a submultiplicative infinity-norm bound below 2**63-1.
"""
from itertools import permutations,combinations
from collections import Counter
import json
import numpy as np
from verify_h4_relations import geom,even,pmul,pinv
LIM=2**63-1

def norm(A):return max((sum(abs(int(x)) for x in row) for row in A),default=0)
def mul(A,B):
 assert norm(A)*norm(B)<=LIM
 return A@B

def certify(M,rs):
 n=len(M);I=np.eye(n,dtype=np.int64);acc=I.copy();bounds=[]
 for r in rs:
  factor=M-r*I;bounds.append(norm(acc)*norm(factor));acc=mul(acc,factor)
 assert not np.any(acc)
 out={}
 projectors={}
 for r in rs:
  A=I.copy();den=1
  for s in rs:
   if s!=r:A=mul(A,M-s*I);den*=r-s
  assert n*norm(A)<=LIM
  tr=sum(int(A[i,i]) for i in range(n))
  assert tr%den==0
  out[r]=tr//den;projectors[r]=(A,den)
 assert all(m>0 for m in out.values()) and sum(out.values())==n
 return out,projectors,max(bounds)

def representations():
 A5=[p for p in permutations(range(5)) if even(p)]
 id5=tuple(range(5));I5=np.eye(5,dtype=np.int64)
 R4={}
 for p in A5:
  P=I5[:,p];R4[p]=(P[:,:4]-P[:,4:5])[:4]
 wedge_basis=list(combinations(range(4),2))
 def wedge(M):return np.array([[M[a,c]*M[b,d]-M[a,d]*M[b,c] for c,d in wedge_basis] for a,b in wedge_basis],dtype=np.int64)
 R6={p:wedge(R4[p]) for p in A5}
 # Six Sylow5 subgroups, with conjugation permutation representation.
 c=(1,2,3,4,0);cls5=sorted({pmul(p,pmul(c,pinv(p))) for p in A5});assert len(cls5)==12
 subs=set()
 for a in cls5:
  powers=[];x=id5
  for _ in range(4):x=pmul(x,a);powers.append(x)
  subs.add(tuple(sorted(powers)))
 subs=sorted(subs);assert len(subs)==6
 R5={}
 for p in A5:
  induced=tuple(subs.index(tuple(sorted(pmul(p,pmul(a,pinv(p))) for a in S))) for S in subs)
  P=np.eye(6,dtype=np.int64)[:,induced]
  R5[p]=(P[:,:5]-P[:,5:6])[:5]
 J4=np.eye(4,dtype=np.int64)+np.ones((4,4),dtype=np.int64);J6=wedge(J4)
 J5=np.eye(5,dtype=np.int64)+np.ones((5,5),dtype=np.int64)
 for R,J in [(R4,J4),(R6,J6),(R5,J5)]:
  for p in A5:
   assert np.array_equal(mul(mul(R[p].T,J),R[p]),J)
   for q in A5:assert np.array_equal(mul(R[p],R[q]),R[pmul(p,q)])
 classsum=sum(R6[p] for p in cls5);I6=np.eye(6,dtype=np.int64)
 assert not np.any((classsum-2*I6)%2)
 C=(classsum-2*I6)//2;assert np.array_equal(mul(C,C),5*I6)
 assert np.trace(C)==0
 return A5,{'chi1':{p:np.ones((1,1),dtype=np.int64) for p in A5},'chi3sum':R6,'chi4':R4,'chi5':R5}, {'chi1':np.ones((1,1),dtype=np.int64),'chi3sum':J6,'chi4':J4,'chi5':J5},C

def verify(d=None):
 d=geom() if d is None else d;A5,reps,metrics,C=representations();pairs=d['pairs'];dec=d['dec'];frames=d['frames8'];fl=d['fl'];act=d['act'];K=d['K'];pl=d['pl'];rr=d['RR']
 refs=[rr[a][rr[b]] for a,b in pairs]
 # Group-ring class sum: off diagonal = K(0), diagonal =1+class(2A).
 off=0
 for j,X in enumerate(pairs):
  observed={};diagonal=[]
  for g in refs:
   Y=act(g,X);i=pl[Y];p=d['transport'](g,X,Y)
   if i==j:diagonal.append(p)
   else:
    assert i not in observed;observed[i]=p;assert K[i][j]==(0,p);off+=1
  assert len(observed)==44 and len(diagonal)==16
  assert diagonal.count(tuple(range(5)))==1
  invols=[p for p in A5 if p!=tuple(range(5)) and pmul(p,p)==tuple(range(5))]
  assert len(invols)==15 and sorted(diagonal)==sorted(invols+[tuple(range(5))])
 result={'reflection_class_sum_offdiagonal_matches':off,'constant_terms':{},'spectra':{},'int64_product_norm_bounds':{},'rank_method':'squarefree annihilator + exact Lagrange projector traces','pair_representations':'trivial; exterior square of A5 4dim standard; standard4; Sylow5 standard5'}
 roots={'chi1':[-44,-4,1,4],'chi3sum':[-4,8,11,16],'chi4':[-14,-11,1,13],'chi5':[-16,-11,-8,4]}
 for name,rep in reps.items():
  dim=len(rep[tuple(range(5))]);M=np.zeros((60*dim,60*dim),dtype=np.int64)
  for i in range(60):
   for j in range(60):
    if i!=j and K[i][j] is not None:M[dim*i:dim*(i+1),dim*j:dim*(j+1)]=rep[K[i][j][1]]
  # Verify the diagonal coefficient class sum is scalar in each representation.
  S=sum(rep[p] for p in invols);c=int(S[0,0]);assert np.array_equal(S,c*np.eye(dim,dtype=np.int64))
  result['constant_terms'][name]=c
  # Build Z independently from the full 60 reflections and compare.
  Z=np.zeros_like(M)
  for g in refs:
   for j,X in enumerate(pairs):
    Y=act(g,X);i=pl[Y];p=d['transport'](g,X,Y)
    Z[dim*i:dim*(i+1),dim*j:dim*(j+1)]+=rep[p]
  assert np.array_equal(M,Z-(1+c)*np.eye(len(M),dtype=np.int64))
  ms,P,bound=certify(-M,roots[name]);result['spectra'][name]=ms;result['int64_product_norm_bounds'][name]=bound
  if name=='chi3sum':
   CC=np.kron(np.eye(60,dtype=np.int64),C)
   assert np.array_equal(mul(CC,M),mul(M,CC))
   for r,(numer,den) in P.items():
    T=mul(numer,CC);assert len(T)*norm(T)<=LIM
    assert sum(int(T[i,i]) for i in range(len(T)))==0
    assert ms[r]%2==0
   result['spectra']['chi3_each']={r:v//2 for r,v in ms.items()}
 # Maximal graph, polynomial, eigenspaces, exact graph definition.
 A=np.array([[int(token[0]==1) for token in row] for row in d['Gmax']],dtype=np.int64)
 assert np.array_equal(A,A.T) and np.array_equal(A@np.ones(75,dtype=np.int64),18*np.ones(75,dtype=np.int64))
 ms,P,bound=certify(A,[-5,-2,3,18]);result['spectra']['maximal_adjacency']=ms;result['int64_product_norm_bounds']['maximal_adjacency']=bound
 I=np.eye(75,dtype=np.int64);J=np.ones((75,75),dtype=np.int64)
 assert np.array_equal(np.array([[int(x[0]==2) for x in row] for row in d['Gmax']]),I)
 assert np.array_equal(np.array([[int(x[0]==0) for x in row] for row in d['Gmax']]),J-I-A)
 # The all-ones eigenspace is the unique 18-eigenspace; J acts as 75 there.
 # On its orthogonal complement J acts as zero. Thus the exact graph spectrum
 # gives the full maximal Gram determinant without evaluating a large polynomial.
 from sage.all import PolynomialRing, QQ
 t=PolynomialRing(QQ,'t').gen()
 maximal=(t*t+18*t+56)**ms[18]
 for r in [-5,-2,3]: maximal *= (t*t+r*t-1-r)**ms[r]
 assert maximal==(t-1)**82*(t-4)**25*(t+4)**42*(t+14)
 result['maximal_determinant']='(t - 1)^82 * (t - 4)^25 * (t + 4)^42 * (t + 14)'
 result['pair_determinants_up_to_nonzero_constant']={
     name:' * '.join('(t - (%s))^%s'%(r,m) for r,m in result['spectra'][name].items())
     for name in ['chi1','chi3_each','chi4','chi5']}
 result['status']='PASS'
 return result

if __name__=='__main__':
 print(json.dumps(verify(),indent=2,sort_keys=True))
