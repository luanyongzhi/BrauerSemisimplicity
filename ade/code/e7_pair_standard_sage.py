"""E7 pair standard module, epsilon=+, at delta=-7. Exact modular certificate."""
exec(open('e7_pair_sign.py').read().split('T=np.zeros')[0])
from sage.all import matrix as sage_matrix, GF
def nmod_mat(data,p):
 return sage_matrix(GF(p),data)
simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0]])
sids=[look[canon(a)] for a in simple]
def frame_ref(a,F):
 return tuple(sorted(canon(np.array(v)-(np.array(v)@roots[a])*roots[a]//4) for v in F))
def ray_ref(a,v):
 v=np.array(v);return canon(v-(v@roots[a])*roots[a]//4)
x0=tuple(sorted([look[(2,-2,0,0,0,0,0,0)],look[(2,2,0,0,0,0,0,0)]]));i0=plook[x0]
F0=tuple(sorted(tuple(4*int(j==i) for j in range(8)) for i in range(2,6)))
frames={i0:F0};todo=[i0]
while todo:
 k=todo.pop();i,j=pairs[k]
 for a in sids:
  dest=plook[tuple(sorted((rp[a,i],rp[a,j])))];F=frame_ref(a,frames[k])
  if dest in frames:assert frames[dest]==F
  else:frames[dest]=F;todo.append(dest)
assert len(frames)==945
print('945 four-ray frames consistently transported.',flush=True)
# Fiber maps for all contractions and simple Weyl reflections.
perm=np.zeros((63,945,4),dtype=np.int64)
for a in range(63):
 for k,(i,j) in enumerate(pairs):
  dest=out[a,k]
  if dest<0:continue
  if a in (i,j):perm[a,k]=np.arange(4);continue
  b=i if gram[a,i] else j
  im=[ray_ref(b,ray_ref(a,v)) for v in frames[k]]
  perm[a,k]=[frames[dest].index(v) for v in im]
# Check full defining relations on dimension3780 decorated permutation module.
NN=3780;zero=np.zeros(NN,dtype=np.int64);I=(np.arange(NN),zero)
def mul(a,b):
 ai,ae=a;bi,be=b;di=np.full(NN,-1,dtype=np.int64);de=zero.copy();good=np.flatnonzero(bi>=0);good=good[ai[bi[good]]>=0];di[good]=ai[bi[good]];de[good]=be[good]+ae[bi[good]];return di,de
def prod(*args):
 z=I
 for a in args:z=mul(z,a)
 return z
def eq(a,b):assert all(np.array_equal(u,v) for u,v in zip(a,b))
def E(a):
 d=np.where(out[a,:,None]<0,-1,4*out[a,:,None]+perm[a]).reshape(-1);e=np.repeat(ex[a],4);return d,e
rsim=[];esim=[]
for a in sids:
 d=[]
 for k,(i,j) in enumerate(pairs):
  dest=plook[tuple(sorted((rp[a,i],rp[a,j])))];d.extend(4*dest+frames[dest].index(ray_ref(a,v)) for v in frames[k])
 rsim.append((np.array(d),zero));esim.append(E(a))
for i in range(7):
 r,e=rsim[i],esim[i];eq(prod(r,r),I);eq(prod(r,e),e);eq(prod(e,r),e);eq(prod(e,e),(e[0],np.where(e[0]<0,0,e[1]+1)))
 for j in range(i):
  rr,ee=rsim[j],esim[j]
  if simple[i]@simple[j]==0:
   for a,b in [(r,rr),(e,ee),(r,ee),(e,rr)]:eq(prod(a,b),prod(b,a))
  else:
   eq(prod(r,rr,r),prod(rr,r,rr));eq(prod(rr,e,rr),prod(r,ee,r));eq(prod(rr,r,ee),prod(e,ee));eq(prod(r,rr,e),prod(ee,e))
for i,a in enumerate(sids):
 for b in range(63):eq(prod(rsim[i],E(b),rsim[i]),E(int(rp[a,b])))
print('All CFW relations and 441 root-conjugation identities verified symbolically.',flush=True)
T=np.zeros((2835,2835),dtype=np.int64);eye=np.eye(4,dtype=np.int64)
for row,(a,b) in enumerate(pairs):
 for col in range(945):
  mid=out[b,col]
  if mid<0:continue
  dest=out[a,mid]
  if dest<0:continue
  assert dest==row
  pi=perm[a,mid][perm[b,col]];P=eye[:,pi];M=(P[:,:3]-P[:,3:4])[:3,:]
  T[3*row:3*row+3,3*col:3*col+3]=M*(-7)**int(ex[b,col]+ex[a,mid])
print('Integer 2835-square matrix constructed.',flush=True)
p=1009;rank=int(nmod_mat(T.tolist(),p).rank());print('rank mod',p,rank,flush=True)
json.dump(dict(dimension=2835,delta=-7,prime=p,rank_mod_prime=rank,nullity_mod_prime=2835-rank,CFW_relations_verified=True,root_conjugation_identities=441),open('e7_pair_standard_result.json','w'),indent=2)
