"""Exact E8 singleton branching, using CHEVIE q=1 tables and root permutations.
Source tables: https://github.com/jmichel7/Chevie.jl/tree/master/src/tbl
"""
import ast,re,json,itertools,collections,hashlib,csv
import numpy as np
from pathlib import Path

def bracket(s,start):
 start=s.index('[',start); depth=0
 for i in range(start,len(s)):
  if s[i]=='[': depth+=1
  elif s[i]==']':
   depth-=1
   if depth==0:return s[start:i+1]
 raise ValueError

def parse(n):
 s=Path(f'weyle{n}.jl').read_text()
 pos=s.index(':classwords,')
 pos=s.index('\n',pos)
 raw=ast.literal_eval(bracket(s,pos).replace('Int[]','[]'))
 cox=list(range(1,9)) if n==8 else list(range(7,0,-1))
 words=[[j for i in w for j in (cox if i==0 else [i])] for w in raw]
 vp=ast.literal_eval(bracket(s,s.index('}}}',s.index(':vpolheckeirreducibles'))+3).replace('nothing','None'))
 char=[]
 for row in vp:
  if row is None:char.append([v*(-1)**len(w) for v,w in zip(char[-1],words)])
  else:char.append([sum(coeff) for coeff,val in row])
 classes=ast.literal_eval(bracket(s,s.index(':classes =>')))
 params=ast.literal_eval(bracket(s,s.index(':charparams,')))
 chars=np.array(char,dtype=np.int64);cls=np.array(classes,dtype=np.int64)
 size=sum(classes)
 assert chars.shape==(len(classes),len(classes))
 assert np.array_equal((chars*cls)@chars.T,size*np.eye(len(classes),dtype=np.int64))
 assert list(chars[:,0])==[p[0] for p in params]
 return words,chars,cls,params

w7,ch7,cl7,pa7=parse(7)
w8,ch8,cl8,pa8=parse(8)
# Doubled Euclidean coordinates avoid floating point.
roots=[]
for i,j in itertools.combinations(range(8),2):
 for a,b in itertools.product([-2,2],repeat=2):
  v=[0]*8;v[i]=a;v[j]=b;roots.append(tuple(v))
for signs in itertools.product([-1,1],repeat=8):
 if sum(x<0 for x in signs)%2==0:roots.append(signs)
roots=np.array(roots,dtype=np.int64); lookup={tuple(v):i for i,v in enumerate(roots)}
simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0],[0,0,0,0,0,-2,2,0]])
def refl(v):
 nums=np.outer(roots@v,v)
 assert np.all(nums%4==0)
 return np.array([lookup[tuple(r)] for r in roots-nums//4])
gens=[refl(v) for v in simple]
ident=np.arange(240)
def word(w):
 p=ident.copy()
 for i in w:p=gens[i-1][p]
 return p

def cycles(p):
 seen=set();lens=[]
 for i in range(len(p)):
  if i in seen:continue
  j=i;k=0
  while j not in seen:seen.add(j);k+=1;j=int(p[j])
  if k>1:lens.append(k)
 return tuple(sorted(collections.Counter(lens).items()))
cp=[cycles(word(w)) for w in w8]
assert len(set(cp))==112, 'E8 root cycle structures must distinguish classes'
fusdict={x:i for i,x in enumerate(cp)}
alpha=np.array([0,0,0,0,0,0,2,2]); ralpha=refl(alpha)
assert np.all(simple[:7]@alpha==0)
assert all(np.array_equal(ralpha[g],g[ralpha]) for g in gens[:7])
fus0=[fusdict[cycles(word(w))] for w in w7]
fus1=[fusdict[cycles(ralpha[word(w)])] for w in w7]
# Character induction from C2 x E7 with trivial C2 factor.
num=((ch8[:,fus0]+ch8[:,fus1])*cl7)@ch7.T
den=2*sum(cl7)
assert np.all(num%den==0)
b=num//den
assert np.all(b>=0)
assert np.array_equal(ch8[:,0]@b,120*ch7[:,0])
assert np.array_equal(b@ch7[:,0],(ch8[:,0]+ch8[:,next(i for i,w in enumerate(w8) if len(w)==1)])//2)
# Reflections classes are represented by length-one words.
r7=next(i for i,w in enumerate(w7) if len(w)==1)
r8=next(i for i,w in enumerate(w8) if len(w)==1)
assert cl7[r7]==63 and cl8[r8]==120
assert np.all((63*ch7[:,r7])%ch7[:,0]==0)
assert np.all((120*ch8[:,r8])%ch8[:,0]==0)
c7=63*ch7[:,r7]//ch7[:,0];c8=120*ch8[:,r8]//ch8[:,0]
polys=[];sigma=set()
for j in range(60):
 factors=collections.Counter()
 for i in range(112):
  if b[i,j]:factors[int(1+c7[j]-c8[i])]+=int(b[i,j]*ch8[i,0])
 assert sum(factors.values())==120*ch7[j,0]
 assert sum(m*e for m,e in factors.items())==0
 assert sum(m*m*e for m,e in factors.items())==56*120*ch7[j,0]
 sigma.update(m for m in factors if m)
 polys.append(dict(sorted(factors.items())))
assert {-56,-44,-32,-20,-8,4,16,28}<=sigma
# Identify the 35-dimensional harmonic quadratics and sign twist exactly.
square_fusion=[fusdict[cycles(word(w)[word(w)])] for w in w8]
ref=ch8[pa8.index([8,1])]
harm=(ref*ref+ref[square_fusion])//2-1
assert np.array_equal(harm,ch8[pa8.index([35,2])])
sgn=ch8[pa8.index([1,120])]
assert np.array_equal(harm*sgn,ch8[pa8.index([35,74])])
result={'source':'CHEVIE weyle7.jl and weyle8.jl, q=1','sha256':{f'weyle{n}.jl':hashlib.sha256(Path(f'weyle{n}.jl').read_bytes()).hexdigest() for n in [7,8]},'checks':['both row orthogonality identities','dimensions match charparams','112 distinct E8 root cycle types','all 120 subgroup class fusions unique','all branching entries nonnegative integers','all 60 induced dimensions correct','all 60 determinant traces zero','exterior-power subset present','reflection-fixed dimension of every restricted E8 character','all 60 quadratic trace identities','35-dimensional harmonic quadratic and sign-twist character identities'],'E7_labels':pa7,'E8_labels':pa8,'fusion_identity':fus0,'fusion_reflection':fus1,'reflection_scalars_E7':c7.tolist(),'reflection_scalars_E8':c8.tolist(),'branching':b.tolist(),'determinant_factors':polys,'singleton_bad_set':sorted(sigma),'max_branching_multiplicity':int(b.max())}
Path('e8_singleton_results.json').write_text(json.dumps(result,indent=2))
print('E8 singleton set',sorted(sigma),'size',len(sigma),'maxmult',b.max())
print('within [-56,56] remaining candidates',sorted(set(range(-56,57))-sigma-{0}))

with open('e8_singleton_branching.csv','w') as f:
 w=csv.writer(f);w.writerow(['degree','b','c']+[str(p) for p in pa7]);w.writerows([[*pa8[i],int(c8[i]),*b[i].tolist()] for i in range(112)])
