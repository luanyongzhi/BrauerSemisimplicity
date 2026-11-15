"""Exact E8 pair-layer W-branching and central-collision screen.
N=(W(B6)xW(B2)) intersect W(D8); A_X=kernel(N -> S6).
Uses E8 character/root parser in e8_singleton and symmetric characters in e7_singleton.
"""
import itertools,math,json,csv
from collections import Counter
from fractions import Fraction
import numpy as np
import e8_singleton as e8
from e7_singleton import parts,z,symchar

def signed_types(n):
 for k in range(n+1):
  for mu in parts(k):
   for nu in parts(n-k):yield mu,nu

def classsize(n,mu,nu):return 2**n*math.factorial(n)//(2**(len(mu)+len(nu))*z(mu)*z(nu))

def perm(mu,nu,lam,kap,split=False):
 p=list(range(8));sg=[1]*8
 for pos,blocks in [(0,[(x,False) for x in mu]+[(x,True) for x in nu]),(6,[(x,False) for x in lam]+[(x,True) for x in kap])]:
  for length,negative in blocks:
   for i in range(pos,pos+length):p[i]=pos+(i-pos+1)%length
   if negative:sg[pos]=-1
   pos+=length
 if split:sg=[sg[i]*(-1 if i==0 else 1)*(-1 if p[i]==0 else 1) for i in range(8)]
 ans=[]
 for r in e8.roots:
  out=[0]*8
  for i in range(8):out[p[i]]=sg[i]*r[i]
  ans.append(e8.lookup[tuple(out)])
 return np.array(ans)
classes=[]
for mu,nu in signed_types(6):
 for lam,kap in signed_types(2):
  if (len(nu)+len(kap))%2:continue
  sz=classsize(6,mu,nu)*classsize(2,lam,kap)
  split=not nu and not kap and all(x%2==0 for x in mu+lam)
  for which in range(2 if split else 1):
   fus=e8.fusdict[e8.cycles(perm(mu,nu,lam,kap,bool(which)))]
   classes.append(dict(mu=mu,nu=nu,lam=lam,kap=kap,size=sz//(2 if split else 1),fusion=fus,split=which))
assert sum(c['size'] for c in classes)==184320
labels=parts(6);cols=[];c2=[];dims=[]
for th in labels:
 chi=np.array([symchar(th,tuple(sorted(c['mu']+c['nu'],reverse=True))) for c in classes])
 raw=(e8.ch8[:,[c['fusion'] for c in classes]]*np.array([c['size'] for c in classes]))@chi
 assert np.all(raw%184320==0)
 col=raw//184320
 assert np.all(col>=0)
 dim=symchar(th,(1,)*6)
 assert col@e8.ch8[:,0]==3780*dim
 cols.append(col);dims.append(dim);c2.append(2*sum(j-i for i,n in enumerate(th) for j in range(n)))
B=np.array(cols).T
# First C2xE7 induction and collision filter data.
r7=json.load(open('e7_singleton_results.json'));sigma7=set(r7['exceptional_set']);sigma8=set(e8.result['singleton_bad_set'])
cs={0:set(e8.c8),1:set(e8.c7),2:set(c2),4:{-12,0,12},8:{0}}
by_theta=[];s0all=set()
for j,th in enumerate(labels):
 vals=set()
 for i,mult in enumerate(B[:,j]):
  if not mult:continue
  d=1+Fraction(int(c2[j]-e8.c8[i]),2)
  if d.denominator==1 and d:vals.add(int(d))
 s0all|=vals
 by_theta.append(dict(partition=th,dimension=dims[j],central_scalar=c2[j],s0_candidates=sorted(vals),s0_outside_singleton=sorted(vals-sigma8)))
remaining={};allhigher=set()
for t in [2,4,8]:
 for s in [i for i in cs if i<t]:
  if (t,s)==(2,0):vals=s0all
  else:
   vals=set()
   for ct in cs[t]:
    for cx in cs[s]:
     if s==0:
      if cx in [-120,-90]:continue
      if cx==120 and ct!=max(cs[t]):continue
      if t==4 and cx==-72:continue # sign-twisted harmonic quadratic exclusion
      if t==4 and cx==72 and ct!=12:continue # harmonic quadratic invariants only trivial S3
     d=1+Fraction(int(ct-cx),t-s)
     if d.denominator==1 and d:
      if (t,s)==(2,1) and int(d) not in sigma7:continue
      vals.add(int(d))
  allhigher|=vals;remaining[f'{t}->{s}']=sorted(vals-sigma8)
out=dict(normalizer_order=184320,normalizer_class_count=len(classes),normalizer_classes=classes,coefficient_labels=labels,W_labels=e8.pa8,branching=B.tolist(),pair_s0_candidates=sorted(s0all),pair_s0_outside_singleton=sorted(s0all-sigma8),pair_s0_by_coefficient=by_theta,remaining_by_collision=remaining,remaining_total=sorted(allhigher-sigma8))
json.dump(out,open('e8_pair_branching_results.json','w'),indent=2)
with open('e8_pair_branching.csv','w') as f:
 w=csv.writer(f);w.writerow(['degree','b','c']+[str(th) for th in labels]);w.writerows([[*e8.pa8[i],int(e8.c8[i]),*B[i].tolist()] for i in range(112)])
print('E8 pair s0 candidates',sorted(s0all))
print('E8 remaining',remaining,'total',sorted(allhigher-sigma8))
