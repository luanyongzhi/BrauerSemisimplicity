"""Exact pair sign module T(-7) in type E7, with modular rank certificate."""
import e7_singleton as s
import numpy as np,itertools,json
roots=np.array([r for r in s.roots if next(x for x in r if x)>0]); n=len(roots)
look={tuple(r):i for i,r in enumerate(roots)}
def canon(r):return tuple(r if next(x for x in r if x)>0 else -r)
gram=roots@roots.T
pairs=[(i,j) for i in range(n) for j in range(i+1,n) if gram[i,j]==0]
assert len(pairs)==945
plook={p:i for i,p in enumerate(pairs)}
eta=[]
for i,j in pairs:
 orth=[k for k in range(n) if gram[i,k]==gram[j,k]==0]
 iso=[k for k in orth if all(k==l or gram[k,l]==0 for l in orth)]
 assert len(orth)==13 and len(iso)==1
 eta.append(iso[0])
eta=np.array(eta)
# reflection permutations of signed/canonical root lines
rp=[];rs=[]
for a in roots:
 im=roots-np.outer(roots@a,a)//4
 rp.append(np.array([look[canon(v)] for v in im]))
 rs.append(np.array([1 if next(x for x in v if x)>0 else -1 for v in im]))
rp=np.array(rp);rs=np.array(rs)
N=len(pairs);out=np.full((63,N),-1,dtype=int);sg=np.zeros((63,N),dtype=int);ex=np.zeros((63,N),dtype=int)
for a in range(n):
 for k,(i,j) in enumerate(pairs):
  if a in (i,j):out[a,k]=k;sg[a,k]=1;ex[a,k]=1
  elif gram[a,i]==gram[a,j]==0:continue
  else:
   b=i if gram[a,i]!=0 else j
   wi,wj=rp[b,rp[a,i]],rp[b,rp[a,j]]
   dest=plook[tuple(sorted((wi,wj)))]
   out[a,k]=dest;sg[a,k]=rs[a,eta[k]]*rs[b,rp[a,eta[k]]]
   assert rp[b,rp[a,eta[k]]]==eta[dest]
T=np.zeros((N,N),dtype=np.int64)
for k,(a,b) in enumerate(pairs):
 for col in range(N):
  mid=out[b,col]
  if mid<0:continue
  dest=out[a,mid]
  if dest<0:continue
  assert dest==k
  T[k,col]=int(sg[b,col]*sg[a,mid])*(-7)**int(ex[b,col]+ex[a,mid])
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
json.dump(dict(dimension=N,delta=-7,prime=1009,rank_mod_prime=rank,nullity_mod_prime=N-rank),open('e7_pair_sign_result.json','w'),indent=2)
print('E7 pair sign rank at -7 mod1009:',rank,'nullity',N-rank)
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
 ri=np.array([plook[tuple(sorted((rp[a,i],rp[a,j])))] for i,j in pairs])
 rc=rs[a,eta];rsim.append((ri,rc,zero))
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
