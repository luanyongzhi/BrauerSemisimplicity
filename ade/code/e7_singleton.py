"""Exact E7 singleton calculation. Source character data: CHEVIE weyle7.jl.
D6 degenerate pairs are aggregated, which preserves obstruction unions.
"""
import ast,re,itertools,math,json,csv
from functools import lru_cache
from collections import Counter
import numpy as np
src=open('weyle7.jl').read()
def extract_after(marker):
 s=src.split(marker,1)[1]; start=s.index('['); depth=0
 for i in range(start,len(s)):
  depth+=(s[i]=='[')-(s[i]==']')
  if depth==0:return s[start:i+1]
vp=ast.literal_eval(extract_after('Union{Nothing, Vector{Tuple{Vector{Int64}, Int64}}}').replace('nothing','None'))
chars=[]
for row in vp:
 if row is None: chars.append([x*(1 if j<30 else -1) for j,x in enumerate(chars[-1])])
 else: chars.append([sum(c) for c,v in row])
chars=np.array(chars,dtype=np.int64)
classes=ast.literal_eval(extract_after(':classes => '))
assert chars.shape==(60,60)
assert np.array_equal((chars*np.array(classes))@chars.T,2903040*np.eye(60,dtype=np.int64))
params=ast.literal_eval(extract_after(':charparams,'))
assert list(chars[:,0])==[x[0] for x in params]
cw=63*chars[:,31]//chars[:,0]
assert np.array_equal(cw*chars[:,0],63*chars[:,31])
cycles=extract_after(':cyclestructure,').replace('Pair{Int,Int}[]','[]')
cycles=ast.literal_eval(re.sub(r'(\d+)\s*=>\s*(\d+)',r'(\1,\2)',cycles))
cycles=[tuple(x) for x in cycles]
# Roots scaled by 2 in R8, perpendicular to e7+e8.
roots=[]
for i,j in itertools.combinations(range(6),2):
 for a,b in itertools.product([-2,2],repeat=2):
  r=[0]*8;r[i]=a;r[j]=b;roots.append(tuple(r))
roots += [(0,)*6+(2,-2),(0,)*6+(-2,2)]
for r in itertools.product([-1,1],repeat=8):
 if r[6]==-r[7] and sum(x==-1 for x in r)%2==0:roots.append(r)
assert len(roots)==126
lookup={r:i for i,r in enumerate(roots)}
neg=np.array([lookup[tuple(-x for x in r)] for r in roots])
def cyc(p):
 seen=set();count=Counter()
 for i in range(len(p)):
  if i in seen:continue
  j=i;l=0
  while j not in seen:seen.add(j);l+=1;j=p[j]
  if l>1:count[l]+=1
 return tuple(sorted(count.items()))
def fusion(p):
 c=cyc(p); opts=[i for i,v in enumerate(cycles) if v==c]
 if len(opts)==1:return opts[0]
 c2=cyc(neg[p]); ds=[((2,51),),((2,59),),((3,6),(6,17)),((2,2),(3,2),(6,19)),((2,7),(4,25)),((2,11),(4,25))]
 result=[2,4,19,20,39,42][ds.index(c2)]
 assert result in opts
 return result
@lru_cache(None)
def parts(n,cap=99):
 if n==0:return ((),)
 return tuple((k,)+p for k in range(min(cap,n),1-1,-1) for p in parts(n-k,k))
def z(p):
 return math.prod(k**v*math.factorial(v) for k,v in Counter(p).items())
def rootperm(mu,nu,split=False,reflect=False):
 p=list(range(8));sg=[1]*8;pos=0
 for length,negative in [(x,False) for x in mu]+[(x,True) for x in nu]:
  for i in range(pos,pos+length):p[i]=pos+(i-pos+1)%length
  if negative:sg[pos]=-1
  pos+=length
 if split:
  sg=[sg[i]*(-1 if i==0 else 1)*(-1 if p[i]==0 else 1) for i in range(8)]
 if reflect:p[6],p[7]=7,6
 ans=[]
 for r in roots:
  out=[0]*8
  for i in range(8):out[p[i]]=sg[i]*r[i]
  ans.append(lookup[tuple(out)])
 return np.array(ans)
dclasses=[]
for k in range(7):
 for mu in parts(k):
  for nu in parts(6-k):
   if len(nu)%2:continue
   size=2**6*math.factorial(6)//(2**(len(mu)+len(nu))*z(mu)*z(nu))
   split=not nu and all(x%2==0 for x in mu)
   for which in range(2 if split else 1):
    dclasses.append((mu,nu,size//(2 if split else 1),fusion(rootperm(mu,nu,bool(which))),fusion(rootperm(mu,nu,bool(which),True))))
assert sum(c[2] for c in dclasses)==23040
@lru_cache(None)
def symchar(lam,mu):
 if not mu:return int(not lam)
 k=mu[0];total=0
 cells={(i,j) for i,x in enumerate(lam) for j in range(x)}
 for rho in parts(sum(lam)-k):
  if len(rho)>len(lam) or any(x>lam[i] for i,x in enumerate(rho)):continue
  skew=cells-{(i,j) for i,x in enumerate(rho) for j in range(x)}
  if any({(i,j),(i+1,j),(i,j+1),(i+1,j+1)}<=skew for i,j in skew):continue
  seen={next(iter(skew))}
  while True:
   more=seen|{(i+a,j+b) for i,j in seen for a,b in [(1,0),(-1,0),(0,1),(0,-1)] if (i+a,j+b) in skew}
   if more==seen:break
   seen=more
  if seen!=skew:continue
  total+=(-1)**(len({i for i,j in skew})-1)*symchar(rho,mu[1:])
 return total
@lru_cache(None)
def wreath(a,b,mu,nu):
 cycles=[(x,1) for x in mu]+[(x,-1) for x in nu];v=0
 for mask in range(1<<len(cycles)):
  ac=tuple(sorted((x for i,(x,s) in enumerate(cycles) if mask>>i&1),reverse=True))
  if sum(ac)!=sum(a):continue
  bc=tuple(sorted((x for i,(x,s) in enumerate(cycles) if not mask>>i&1),reverse=True))
  sign=math.prod(s for i,(x,s) in enumerate(cycles) if not mask>>i&1)
  v+=sign*symchar(a,ac)*symchar(b,bc)
 return v
bip=[]
for k in range(7):
 for a in parts(k):
  for b in parts(6-k):
   if a>=b:bip.append((a,b))
branch=[];obs=set();data=[]
for a,b in bip:
 ch=np.array([wreath(a,b,mu,nu) for mu,nu,sz,f,g in dclasses],dtype=np.int64)
 dim=wreath(a,b,(1,)*6,())
 ctheta=2*sum(j-i for part in [a,b] for i,x in enumerate(part) for j in range(x))
 assert 30*wreath(a,b,(2,1,1,1,1),())==ctheta*dim
 values=sum((sz*ch[i]*(chars[:,f]+chars[:,g]) for i,(mu,nu,sz,f,g) in enumerate(dclasses)),np.zeros(60,dtype=np.int64))
 assert all(values%46080==0)
 values=values//46080
 assert all(values>=0)
 assert values@chars[:,0]==63*dim
 branch.append(values)
 roots0=sorted({int(1+ctheta-cw[i]) for i,v in enumerate(values) if v})
 obs.update(roots0)
 factors=Counter()
 for i,v in enumerate(values):
  if v:factors[int(1+ctheta-cw[i])]+=int(v*chars[i,0])
 assert sum(factors.values())==63*dim
 assert sum(m*e for m,e in factors.items())==0
 data.append(dict(a=a,b=b,aggregated_split=a==b,dimension=int(dim),c=int(ctheta),roots=roots0,multiplicities=values.tolist(),determinant_factors=dict(sorted(factors.items()))))
B=np.array(branch).T
# D6 character orthogonality; aggregated split pairs have norm two.
D=np.array([[wreath(a,b,mu,nu) for mu,nu,sz,f,g in dclasses] for a,b in bip],dtype=np.int64)
assert np.array_equal((D*np.array([c[2] for c in dclasses]))@D.T,23040*np.diag([2 if a==b else 1 for a,b in bip]))
obs.discard(0)
print('D6 classes:',len(dclasses),'aggregated labels:',len(bip),'max aggregated multiplicity:',B.max())
print('Singleton exceptional set:',sorted(obs))
print('central E7 scalars:',sorted(set(cw)))
json.dump(dict(exceptional_set=sorted(obs),e7_parameters=params,e7_c=cw.tolist(),d6_classes=dclasses,columns=data),open('e7_singleton_results.json','w'),indent=2)
with open('e7_singleton_branching_aggregated.csv','w') as f:
 w=csv.writer(f);w.writerow(['degree','b','c']+[str((a,b)) for a,b in bip]);w.writerows([[*params[i],int(cw[i]),*B[i].tolist()] for i in range(60)])
