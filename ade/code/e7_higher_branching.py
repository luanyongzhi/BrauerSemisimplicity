"""Exact stabilizer-induced branching of E7 t=3,4 layers."""
import e7_singleton as s
import numpy as np,itertools,json
from collections import deque
R=np.array(s.roots,dtype=np.int16); ident=np.arange(126,dtype=np.int16)
def refl(a):
 a=np.array(a); out=R-(R@a)[:,None]*a[None,:]//4
 return np.array([s.lookup[tuple(r)] for r in out],dtype=np.int16)
e=np.eye(8,dtype=np.int16)*2
simple=[e[i]-e[i+1] for i in range(5)]+[e[4]+e[5],np.array([1,1,1,1,1,-1,1,-1])]
gens=[refl(a) for a in simple]
canon=np.minimum(ident,s.neg)
def lines(vs):return tuple(sorted(int(canon[s.lookup[tuple(v)]]) for v in vs))
def act(p,X):return tuple(sorted(int(canon[p[i]]) for i in X))
def inverse(p):
 q=np.empty_like(p);q[p]=ident;return q

def stabilizer(X,expected):
 trans={X:ident}; queue=deque([X]); stgens={}
 while queue:
  Y=queue.popleft();g=trans[Y]
  for a in gens:
   ag=a[g];Z=act(a,Y)
   if Z not in trans:trans[Z]=ag;queue.append(Z)
   else:
    h=inverse(trans[Z])[ag]
    if not np.array_equal(h,ident):stgens[h.tobytes()]=h
 print('orbit',len(trans),'Schreier generators',len(stgens),flush=True)
 # Incrementally add only generators outside current subgroup.
 group={ident.tobytes():ident}; kept=[]
 for g in stgens.values():
  if g.tobytes() in group:continue
  kept.append(g);queue=deque(list(group.values()))
  while queue:
   h=queue.popleft()
   for a in kept:
    ah=a[h];k=ah.tobytes()
    if k not in group:group[k]=ah;queue.append(ah)
  if len(group)==expected:break
 assert len(group)==expected,(len(group),expected)
 return list(group.values())

result={}
for t,vs,order in [(3,[e[0]+e[1],e[0]-e[1],e[6]-e[7]],9216),(4,[e[0]+e[1],e[0]-e[1],e[2]+e[3],e[2]-e[3]],3072)]:
 X=lines(vs);G=stabilizer(X,order)
 if t==3:
  frames=[lines([e[a]+e[b],e[a]-e[b],e[c]+e[d],e[c]-e[d]]) for a,b,c,d in [(2,3,4,5),(2,4,3,5),(2,5,3,4)]]
  def quotient(g):
   p=[frames.index(act(g,F)) for F in frames]
   fix=sum(i==p[i] for i in range(3)); sign=1 if fix in [0,3] else -1
   return np.array([1,sign,fix-1])
  cs=[12,-12,0];dims=[1,1,2]
 else:
  perp=[e[4]+e[5],e[4]-e[5],e[6]-e[7]]; pos={s.lookup[tuple(a)] for a in perp};neg={int(s.neg[i]) for i in pos}
  def quotient(g):
   images=[int(g[i]) for i in pos];assert set(images)<=pos|neg
   sg=(-1)**sum(i in neg for i in images)
   return np.array([1,sg])
  cs=[3,-3];dims=[1,1]
 # Certify CFW A_X lies in quotient kernel, and has the expected order.
 xrefs=[refl(v) for v in vs]; agens=list(xrefs)
 for i,j in itertools.combinations(range(t),2):
  for a in R:
   if abs(int(a@vs[i]))==4 and abs(int(a@vs[j]))==4:
    ra=refl(a); agens.append(ra[xrefs[i][xrefs[j][ra]]])
 adict={ident.tobytes():ident}; kept=[]
 for a in agens:
  assert act(a,X)==X
  assert np.array_equal(quotient(a),dims)
  if a.tobytes() in adict:continue
  kept.append(a);queue=deque(list(adict.values()))
  while queue:
   h=queue.popleft()
   for g in kept:
    gh=g[h];key=gh.tobytes()
    if key not in adict:adict[key]=gh;queue.append(gh)
 assert len(adict)==1536,(t,len(adict))
 print('certified A_X order',len(adict),flush=True)
 acc=np.zeros((60,len(cs)),dtype=np.int64);qgram=np.zeros((len(cs),len(cs)),dtype=np.int64)
 for g in G:
  q=quotient(g);acc+=np.outer(s.chars[:,s.fusion(g)],q);qgram+=np.outer(q,q)
 assert np.array_equal(qgram,order*np.eye(len(cs),dtype=int)),qgram
 assert np.all(acc%order==0);B=acc//order;assert np.all(B>=0)
 assert np.array_equal(s.chars[:,0]@B,(2903040//order)*np.array(dims))
 cols=[];allvals=set()
 for j,c in enumerate(cs):
  vals={int(1+(c-s.cw[i])//t) for i in range(60) if B[i,j] and (c-s.cw[i])%t==0}-{0}
  allvals|=vals;cols.append(dict(c=c,dimension=dims[j],branching=B[:,j].tolist(),s0_candidates=sorted(vals),extra=sorted(vals-s.obs)))
 result[t]=dict(columns=cols,s0_candidates=sorted(allvals),extra=sorted(allvals-s.obs));print(t,result[t]['extra'],flush=True)
json.dump(result,open('e7_higher_branching.json','w'),indent=2)
