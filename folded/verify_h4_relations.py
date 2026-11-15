#!/usr/bin/env python3
"""Independent H4 relation and group-ring Gram audit.

Uses doubled E8 roots, but constructs maximal frames by orthogonal quadruples,
not the Round12 admissible-closure/BFS routine.  All presentation identities
are checked as monomial maps with a formal scalar t and A5 permutation labels.
Thus the pair-layer verification covers every A5 representation at once.
Run: sage -python verify_h4_relations.py
"""
from collections import deque,Counter
from itertools import combinations,product,permutations
import json
import numpy as np

ID5=tuple(range(5))
def pmul(p,q):return tuple(p[q[i]] for i in range(len(p)))
def pinv(p):return tuple(p.index(i) for i in range(len(p)))
def even(p):return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2==0

def geom():
 roots=[]
 for i,j in combinations(range(8),2):
  for a,b in product((-2,2),repeat=2):
   v=[0]*8;v[i]=a;v[j]=b;roots.append(tuple(v))
 roots += [s for s in product((-1,1),repeat=8) if sum(x<0 for x in s)%2==0]
 roots=np.array(roots,dtype=np.int64);lookup={tuple(r):i for i,r in enumerate(roots)}
 neg=np.array([lookup[tuple(-r)] for r in roots]);canonical=np.minimum(np.arange(240),neg)
 representatives=sorted(set(canonical));r_to_l={r:i for i,r in enumerate(representatives)}
 lines=roots[representatives];dots=lines@lines.T
 RR=[];rfull=[]
 for v in lines:
  a=lines-np.outer(lines@v,v)//4
  RR.append(np.array([r_to_l[canonical[lookup[tuple(w)]]] for w in a]))
  b=roots-np.outer(roots@v,v)//4
  rfull.append(np.array([lookup[tuple(w)] for w in b],dtype=np.uint8))
 simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0],[0,0,0,0,0,-2,2,0]])
 sl=[r_to_l[canonical[lookup[tuple(v)]]] for v in simple]
 P=[tuple(sorted((sl[a],sl[b]))) for a,b in [(1,4),(3,5),(2,6),(0,7)]]
 wg=[RR[a][RR[b]] for a,b in P];wgfull=[rfull[a][rfull[b]] for a,b in P]
 signed_identity=np.arange(240,dtype=np.uint8)
 for g in wgfull:assert np.array_equal(g[g],signed_identity)
 for i,j,m in [(0,1,5),(1,2,3),(2,3,3),(0,2,2),(0,3,2),(1,3,2)]:
  h=wgfull[i][wgfull[j]];power=signed_identity.copy()
  for k in range(1,m+1):
   power=h[power]
   assert np.array_equal(power,signed_identity)==(k==m)
 def act(g,X):return tuple(sorted(int(g[x]) for x in X))
 seen=set(P);q=deque(P)
 while q:
  X=q.popleft()
  for g in wg:
   Y=act(g,X)
   if Y not in seen:seen.add(Y);q.append(Y)
 pairs=sorted(seen);pl={X:i for i,X in enumerate(pairs)}
 assert len(pairs)==60
 assert sorted(x for X in pairs for x in X)==list(range(120))
 perp=lambda X,Y:all(dots[a,b]==0 for a in X for b in Y)
 orth=np.array([[perp(X,Y) for Y in pairs] for X in pairs],dtype=bool)
 assert all(sum(row)==15 for row in orth)
 # Direct incidence construction: enumerate orthogonal root-pair quadruples.
 frames4=[C for C in combinations(range(60),4) if all(orth[a,b] for a,b in combinations(C,2))]
 assert len(frames4)==75
 frames8=[tuple(sorted(x for i in C for x in pairs[i])) for C in frames4]
 frame_order=sorted(range(75),key=lambda i:frames8[i]);frames4=[frames4[i] for i in frame_order];frames8=[frames8[i] for i in frame_order]
 fl={F:i for i,F in enumerate(frames8)}
 D={X:tuple(i for i,F in enumerate(frames8) if set(X)<=set(F)) for X in pairs}
 assert all(len(v)==5 for v in D.values())
 # Choose basepair the first folded simple pair, making E1's corner explicit.
 base=P[0];dec={base:D[base]};q=deque([base])
 while q:
  X=q.popleft()
  for g in wg:
   Y=act(g,X)
   if Y not in dec:
    dec[Y]=tuple(fl[act(g,frames8[f])] for f in dec[X]);q.append(Y)
 def transport(g,X,Y):return tuple(dec[Y].index(fl[act(g,frames8[f])]) for f in dec[X])
 # ALL choices in CFW's r_b r_a rule; they must give identical A5 labels.
 K=[[None]*60 for _ in range(60)];choices_checked=0
 for j,X in enumerate(pairs):
  for i,Y in enumerate(pairs):
   if i==j:K[i][j]=(1,ID5);continue
   if orth[i,j]:continue
   choices=[]
   for b1 in X:
    if not dots[Y[1],b1]:continue
    p1=RR[b1][RR[Y[1]]];X1=act(p1,X)
    for b0 in X1:
     if not dots[Y[0],b0]:continue
     p0=RR[b0][RR[Y[0]]];g=p0[p1]
     assert act(g,X)==Y
     pi=transport(g,X,Y);assert even(pi)
     choices.append(pi);choices_checked+=1
   assert choices and len(set(choices))==1
   K[i][j]=(0,choices[0])
 # Induced W and pair E action. Labels encode destination, degree in t, A5.
 R=[];E=[]
 for g,U in zip(wg,P):
  rr=[];ee=[];i=pl[U]
  for j,X in enumerate(pairs):
   Y=act(g,X);pi=transport(g,X,Y);assert even(pi)
   rr.append((pl[Y],0,pi))
   ee.append(None if K[i][j] is None else (i,*K[i][j]))
  R.append(rr);E.append(ee)
 # Maximal E8 frame actions never require admissible closure.
 def e_line(a,F):
  if a in F:return F,1
  bs=[b for b in F if dots[a,b]];assert bs
  outputs={act(RR[b][RR[a]],F) for b in bs};assert len(outputs)==1
  return outputs.pop(),0
 def e_pair(U,F):
  F,k=e_line(U[1],F);F,h=e_line(U[0],F);return F,k+h
 Rmax=[];Emax=[]
 for g,U in zip(wg,P):
  Rmax.append([(fl[act(g,F)],0,ID5) for F in frames8]);row=[]
  for F in frames8:
   F1,k=e_pair(U,F);assert k%2==0 and k in (0,2)
   row.append((fl[F1],k//2,ID5))
  Emax.append(row)
 # Normalized maximal corner e_F=e_U e_V for ANY two root-pairs in F.
 Gmax=[[None]*75 for _ in range(75)];pair_choice_checks=0
 for i,C in enumerate(frames4):
  for j,F in enumerate(frames8):
   outputs=[]
   for a,b in combinations(C,2):
    F1,k=e_pair(pairs[b],F);F2,h=e_pair(pairs[a],F1)
    assert F2==frames8[i] and (k+h) in (0,2,4)
    outputs.append((k+h)//2);pair_choice_checks+=1
   assert len(set(outputs))==1
   Gmax[i][j]=(outputs[0],ID5)
 # Recover full signed root group to certify order and normalizer/kernel.
 eye=np.arange(240,dtype=np.uint8);group={eye.tobytes():eye};q=deque([eye]);normalizer=[]
 # full action on a canonical positive representative recovers root-line action.
 def full_line_action(g):return np.array([r_to_l[canonical[int(g[r])]] for r in representatives])
 while q:
  g=q.popleft();gl=full_line_action(g)
  if act(gl,base)==base:normalizer.append(transport(gl,base,base))
  for s in wgfull:
   h=s[g];key=h.tobytes()
   if key not in group:group[key]=h;q.append(h)
 assert len(group)==14400
 remaining=set(group);class_sizes=[]
 while remaining:
  key=min(remaining);orbit={key};cq=deque([group[key]])
  while cq:
   g=cq.popleft()
   for s in wgfull:
    h=s[g[s]];hk=h.tobytes()
    if hk not in orbit:orbit.add(hk);cq.append(h)
  class_sizes.append(len(orbit));remaining.difference_update(orbit)
 assert len(class_sizes)==34 and sum(class_sizes)==14400
 assert len(normalizer)==240 and len(set(normalizer))==60
 assert set(Counter(normalizer).values())=={4}
 h4_full_reflections=[rfull[a][rfull[b]] for a,b in pairs]
 rx=h4_full_reflections[pl[base]]
 containing_frame=frames4[D[base][0]]
 ch3=eye.copy()
 for a in containing_frame:
  if a!=pl[base]:ch3=h4_full_reflections[a][ch3]
 assert np.array_equal(rx[rx],eye) and np.array_equal(ch3[ch3],eye)
 assert np.array_equal(rx[ch3],ch3[rx]) and np.array_equal(rx[ch3],neg)
 assert len({eye.tobytes(),rx.tobytes(),ch3.tobytes(),rx[ch3].tobytes()})==4
 for g in [eye,rx,ch3,rx[ch3]]:
  gl=full_line_action(g);assert act(gl,base)==base
  assert transport(gl,base,base)==ID5
 return locals()

# One generalized monomial with scalar t^degree and coefficient pi.
def token_product(a,b):
 if a is None or b is None:return None
 return (a[0]+b[0],pmul(a[1],b[1]))
def token_star(a):return None if a is None else (a[0],pinv(a[1]))
def act_word(word,j,gen):
 state=(j,0,ID5)
 for name in reversed(word):
  if state is None:return None
  dest,deg,pi=state;step=gen[name][dest]
  if step is None:return None
  state=(step[0],deg+step[1],pmul(step[2],pi))
 return state

def relations(R,E):
 n=len(R[0]);gen={f'r{i+1}':R[i] for i in range(4)}|{f'e{i+1}':E[i] for i in range(4)}
 rel=[]
 def add(name,left,right=(),shift=0):
  for j in range(n):
   l=act_word(left,j,gen);r=act_word(right,j,gen)
   if r is not None:r=(r[0],r[1]+shift,r[2])
   assert l==r,(name,j,l,r)
  rel.append(name)
 for i in range(1,5):
  r=f'r{i}';e=f'e{i}'
  add(f'r{i}^2',[r,r]);add(f'r{i}e{i}',[r,e],[e]);add(f'e{i}r{i}',[e,r],[e]);add(f'e{i}^2=t e{i}',[e,e],[e],1)
 for i,j in [(1,3),(1,4),(2,4)]:
  r,s=f'r{i}',f'r{j}';e,f=f'e{i}',f'e{j}'
  add(f'{r}{s}={s}{r}',[r,s],[s,r]);add(f'{e}{f}={f}{e}',[e,f],[f,e])
  add(f'{e}{s}={s}{e}',[e,s],[s,e]);add(f'{f}{r}={r}{f}',[f,r],[r,f])
 for i,j in [(2,3),(3,4)]:
  r,s=f'r{i}',f'r{j}';e,f=f'e{i}',f'e{j}'
  add(f'braid{i}{j}',[r,s,r],[s,r,s]);add(f'{r}{f}{r}={s}{e}{s}',[r,f,r],[s,e,s])
  add(f'{s}{r}{f}={e}{f}',[s,r,f],[e,f]);add(f'{r}{s}{e}={f}{e}',[r,s,e],[f,e])
 add('five_braid',['r1','r2','r1','r2','r1'],['r2','r1','r2','r1','r2'])
 add('five_edge_covariance',['r1','r2','e1','r2','r1'],['r2','r1','e2','r1','r2'])
 add('e1e2e1=e1',['e1','e2','e1'],['e1'])
 add('e1r2e1=e1',['e1','r2','e1'],['e1'])
 add('e1r2r1r2e1=e1',['e1','r2','r1','r2','e1'],['e1'])
 w=['r2','r1','r2','r1','r3','r2','r1','r2','r3','r1','r2','r1','r2','r3','r4']
 add('extra_H4_relation',['e1']+w*5,['e1'])
 return rel

def contravariance(K,R,E):
 n=len(K);checks=0
 for M in R+E:
  # K M: one summand per column. M* K: inverse label per row.
  for i in range(n):
   for j in range(n):
    a=M[j];b=M[i]
    left=None if a is None else token_product(K[i][a[0]],a[1:])
    right=None if b is None else token_product(token_star(b[1:]),K[b[0]][j])
    assert left==right,('contravariance',i,j,left,right)
    checks+=1
 return checks

def verify(d=None):
 d=geom() if d is None else d
 pr=relations(d['R'],d['E']);mr=relations(d['Rmax'],d['Emax'])
 pc=contravariance(d['K'],d['R'],d['E']);mc=contravariance(d['Gmax'],d['Rmax'],d['Emax'])
 return {'status':'PASS','pairs':len(d['pairs']),'frames':len(d['frames8']),
         'pair_relation_identities':len(pr),'pair_relation_columns':len(pr)*60,
         'maximal_relation_identities':len(mr),'maximal_relation_columns':len(mr)*75,
         'pair_group_ring_contravariance_entries':pc,'maximal_contravariance_entries':mc,
         'cfw_choice_paths_checked':d['choices_checked'],
         'maximal_two_pair_choices_checked':d['pair_choice_checks'],
         'full_signed_root_group_order':len(d['group']),
         'full_group_conjugacy_classes':len(d['class_sizes']),
         'basepair_normalizer_order':len(d['normalizer']),
         'normalizer_A5_image_order':len(set(d['normalizer'])),
         'normalizer_A5_kernel_order':4}

if __name__=='__main__':
 print(json.dumps(verify(),indent=2,sort_keys=True))
