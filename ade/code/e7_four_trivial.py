"""Exact pair sign module T(-7) in type E7, with modular rank certificate."""
import e7_singleton as s
import numpy as np,itertools,json
roots=np.array([r for r in s.roots if next(x for x in r if x)>0]); n=len(roots)
look={tuple(r):i for i,r in enumerate(roots)}
def canon(r):return tuple(r if next(x for x in r if x)>0 else -r)
gram=roots@roots.T
start=tuple(sorted(look[canon(v)] for v in [np.array([2,2,0,0,0,0,0,0]),np.array([2,-2,0,0,0,0,0,0]),np.array([0,0,2,2,0,0,0,0]),np.array([0,0,2,-2,0,0,0,0])]))
# all root reflection permutations
rp=[]
for a in roots:
 im=roots-np.outer(roots@a,a)//4
 rp.append(np.array([look[canon(v)] for v in im]))
rp=np.array(rp)
from collections import deque
seen={start};q=deque([start])
while q:
 X=q.popleft()
 for p in rp:
  Y=tuple(sorted(p[list(X)]))
  if Y not in seen:seen.add(Y);q.append(Y)
pairs=sorted(seen);assert len(pairs)==945
plook={p:i for i,p in enumerate(pairs)}
N=len(pairs);out=np.full((63,N),-1,dtype=int);sg=np.zeros((63,N),dtype=int);ex=np.zeros((63,N),dtype=int)
for a in range(n):
 for k,X in enumerate(pairs):
  if a in X:out[a,k]=k;sg[a,k]=1;ex[a,k]=1
  elif all(gram[a,b]==0 for b in X):continue
  else:
   b=next(b for b in X if gram[a,b]!=0)
   dest=plook[tuple(sorted(rp[b,rp[a,list(X)]]))]
   out[a,k]=dest;sg[a,k]=1
T=np.zeros((N,N),dtype=np.int64)
for k,X in enumerate(pairs):
 for col in range(N):
  dest=col;power=0
  for a in X:
   old=dest;dest=out[a,old]
   if dest<0:break
   power+=int(ex[a,old])
  if dest<0:continue
  assert dest==k
  T[k,col]=(-7)**power
# rank over prime field. Full rank certifies determinant nonzero over Z.
def rankmod(A,p=1009):
 A=A.copy()%p;r=0;piv=[]
 for j in range(A.shape[1]):
  choices=np.flatnonzero(A[r:,j])
  if not len(choices):continue
  k=r+choices[0];A[[r,k]]=A[[k,r]]
  inv=pow(int(A[r,j]),-1,p)
  A[r,j:]=A[r,j:]*inv%p
  inds=np.flatnonzero(A[r+1:,j])+r+1
  A[inds,j:]=(A[inds,j:]-A[inds,j,None]*A[r,j:])%p
  piv.append(j);r+=1
  if r==A.shape[0]:break
 return r
rank=rankmod(T)
json.dump(dict(dimension=N,delta=-7,prime=1009,rank_mod_prime=rank,nullity_mod_prime=N-rank),open('e7_four_trivial_result.json','w'),indent=2)
print('E7 four-root trivial rank at -7 mod1009:',rank,'nullity',N-rank)
# Verify CFW defining relations as exact Laurent monomials, not sampled values.
zero=np.zeros(N,dtype=int);one=np.ones(N,dtype=int);identity=(np.arange(N),one,zero)
def mul(A,B):
 ai,ac,ae=A;bi,bc,be=B
 oi=np.full(N,-1,dtype=int);oc=np.zeros(N,dtype=int);oe=np.zeros(N,dtype=int)
 good=np.flatnonzero(bi>=0);good=good[ai[bi[good]]>=0]
 oi[good]=ai[bi[good]];oc[good]=bc[good]*ac[bi[good]];oe[good]=be[good]+ae[bi[good]]
 return oi,oc,oe
def eq(A,B):return all(np.array_equal(x,y) for x,y in zip(A,B))
simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0]])
rsim=[];esim=[]
for root in simple:
 a=look[canon(root)]
 ri=np.array([plook[tuple(sorted(rp[a,list(X)]))] for X in pairs])
 rc=one;rsim.append((ri,rc,zero))
 esim.append((out[a],sg[a],ex[a]))
for i in range(7):
 r,e=rsim[i],esim[i]
 assert eq(mul(r,r),identity)
 assert eq(mul(r,e),e) and eq(mul(e,r),e)
 scaled=(e[0],e[1],np.where(e[0]>=0,e[2]+1,0))
 assert eq(mul(e,e),scaled)
 for j in range(i):
  rj,ej=rsim[j],esim[j]
  if simple[i]@simple[j]==0:
   for a,b in [(r,rj),(e,ej),(r,ej),(e,rj)]:assert eq(mul(a,b),mul(b,a))
  else:
   assert eq(mul(r,mul(rj,r)),mul(rj,mul(r,rj)))
   assert eq(mul(rj,mul(e,rj)),mul(r,mul(ej,r)))
   assert eq(mul(rj,mul(r,ej)),mul(e,ej))
   assert eq(mul(r,mul(rj,e)),mul(ej,e))
print('All E7 CFW defining relations passed as exact monomial identities.')
for root in simple:
 a=look[canon(root)]
 ri=np.array([plook[tuple(sorted(rp[a,list(X)]))] for X in pairs])
 r=(ri,one,zero)
 for b in range(n):
  lhs=mul(r,mul((out[b],sg[b],ex[b]),r))
  c=rp[a,b]
  assert eq(lhs,(out[c],sg[c],ex[c]))
print('All root-operator conjugation identities passed exactly.')
result=json.load(open('e7_four_trivial_result.json'))
result.update(cfw_relations_exact=True,root_conjugations_exact=True,coefficient='trivial',cardinality=4)
json.dump(result,open('e7_four_trivial_result.json','w'),indent=2)
