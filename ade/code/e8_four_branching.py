"""Exact normalizer enumeration and W(E8) branching for the four-root orbit."""
import numpy as np,itertools,collections,json,csv,time
import e8_singleton as e
root=e.roots;idx=e.lookup
neg=np.array([idx[tuple(-r)] for r in root]);line=np.minimum(np.arange(240),neg)
def rvec(coords):return np.array(coords)
Xroots=[rvec([0,0,0,0,2,2,0,0]),rvec([0,0,0,0,2,-2,0,0]),rvec([0,0,0,0,0,0,2,2]),rvec([0,0,0,0,0,0,2,-2])]
X=tuple(sorted(line[idx[tuple(r)]] for r in Xroots));XR=[e.refl(r) for r in Xroots]
perp=[i for i,r in enumerate(root) if i==line[i] and all(r@x==0 for x in Xroots)]
frames=[tuple(q) for q in itertools.combinations(perp,4) if all(root[i]@root[j]==0 for i,j in itertools.combinations(q,2))]
assert len(frames)==3
reflection_frame_actions=collections.Counter()
for i in perp:
 p=e.refl(root[i]);im=tuple(frames.index(tuple(sorted(line[p[list(f)]]))) for f in frames)
 reflection_frame_actions[im]+=1
assert sorted(reflection_frame_actions.values())==[4,4,4]
assert all(sum(i==j for i,j in enumerate(pi))==1 for pi in reflection_frame_actions)

def actionframes(p):return tuple(frames.index(tuple(sorted(line[p[list(f)]]))) for f in frames)
Agen=XR.copy()
for a in [r for i,r in enumerate(root) if i==line[i]]:
 ra=e.refl(a)
 for i,j in itertools.combinations(range(4),2):
  if abs(a@Xroots[i])==4 and abs(a@Xroots[j])==4:
   q=ra[XR[i][XR[j][ra]]]
   assert tuple(sorted(line[q[list(X)]]))==X
   Agen.append(q)
Agen=list({q.tobytes():q for q in Agen}.values())
print('A generators',len(Agen),flush=True)
assert all(actionframes(q)==(0,1,2) for q in Agen)
# Schreier orbit transversal, root permutations in uint8.
trans=[np.arange(240,dtype=np.uint8)];sets=[X];lookup={X:0};pos=0
while pos<len(sets):
 x=sets[pos];p=trans[pos]
 for g in e.gens:
  y=tuple(sorted(line[g[list(x)]]))
  if y not in lookup:
   lookup[y]=len(sets);sets.append(y);trans.append(g[p].astype(np.uint8))
 pos+=1
assert len(sets)==9450
print('orbit',len(sets),flush=True)
inv=[np.argsort(p).astype(np.uint8) for p in trans]
rng=np.random.default_rng(187);ngens=[]
for j in rng.choice(len(sets),size=30,replace=False):
 for g in e.gens:
  y=tuple(sorted(line[g[list(sets[j])]]));q=inv[lookup[y]][g[trans[j]]]
  if not np.array_equal(q,np.arange(240)):ngens.append(q.astype(np.uint8))
# A_X together with W(perpD4) generates N by CFW splitting.
# Use a modest Schreier generating sample and exact closure cardinality check.
ngens=ngens[:12]
assert all(tuple(sorted(line[g[list(X)]]))==X for g in ngens)
seen={np.arange(240,dtype=np.uint8).tobytes():0};elts=[np.arange(240,dtype=np.uint8)];pos=0
while pos<len(elts):
 p=elts[pos]
 for g in ngens:
  q=g[p];key=q.tobytes()
  if key not in seen:seen[key]=len(elts);elts.append(q)
 pos+=1
 if pos%10000==0:print('closure',pos,len(elts),flush=True)
assert len(elts)==73728
print('N order',len(elts),flush=True)
# All CFW annihilator generators lie in the action kernel; since action onto S3
# and |N/A|=6, the action kernel equals A_X.
fp=[actionframes(p) for p in elts];assert len(set(fp))==6
assert all(g.astype(np.uint8).tobytes() in seen for g in Agen)
counts=collections.Counter()
for p,pi in zip(elts,fp):
 fus=e.fusdict[e.cycles(p)]
 fixed=sum(i==j for i,j in enumerate(pi));sgn=(-1)**sum(pi[i]>pi[j] for i in range(3) for j in range(i+1,3))
 counts[(fus,fixed,sgn)]+=1
B=[]
for kind in ['trivial','standard','sign']:
 raw=sum((n*(1 if kind=='trivial' else fix-1 if kind=='standard' else sign)*e.ch8[:,f] for (f,fix,sign),n in counts.items()),np.zeros(112,dtype=np.int64))
 assert np.all(raw%73728==0)
 col=raw//73728;assert np.all(col>=0)
 assert col@e.ch8[:,0]==9450*(2 if kind=='standard' else 1)
 B.append(col)
B=np.array(B).T
cs=[12,0,-12];sigma=set(e.result['singleton_bad_set']);candidates=[]
from fractions import Fraction
for j,ct in enumerate(cs):
 vals=set()
 for i,v in enumerate(B[:,j]):
  if v:
   d=1+Fraction(int(ct-e.c8[i]),4)
   if d.denominator==1 and d:vals.add(int(d))
 candidates.append(sorted(vals))
print('t4 s0',candidates,'outside',sorted(set(sum(candidates,[]))-sigma),flush=True)
result=dict(orbit_size=len(sets),normalizer_order=len(elts),annihilator_generators=len(Agen),quotient_action='three perpendicular D4 orthogonal frames',all_A_generators_fix_frames=True,fusion_counts=[dict(E8_class=f,fixed_frames=fix,sign=sign,count=n) for (f,fix,sign),n in sorted(counts.items())],branching=B.tolist(),coefficient_labels=['trivial','standard','sign'],coefficient_scalars=cs,s0_candidates=candidates,s0_outside_singleton=sorted(set(sum(candidates,[]))-sigma))
json.dump(result,open('e8_four_branching_results.json','w'),indent=2)
with open('e8_four_branching.csv','w') as f:
 w=csv.writer(f);w.writerow(['degree','b','c','trivial','standard','sign']);w.writerows([[*e.pa8[i],int(e.c8[i]),*B[i].tolist()] for i in range(112)])
