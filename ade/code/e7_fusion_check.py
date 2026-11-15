"""Independent check of all 60 E7 class identifications against classwords."""
import e7_singleton as s
import ast,numpy as np
src=open('weyle7.jl').read();raw=src.split(':classwords,map',1)[1].split('\n',1)[1]
i=raw.index('[');lev=0
for j in range(i,len(raw)):
 lev+=(raw[j]=='[')-(raw[j]==']')
 if lev==0:break
words=ast.literal_eval(raw[i:j+1].replace('Int[]','[]'))
words=[[k for x in w for k in (list(range(7,0,-1)) if x==0 else [x])] for w in words]
simple=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0]])
rr=np.array(s.roots)
gens=[np.array([s.lookup[tuple(r)] for r in rr-np.outer(rr@a,a)//4]) for a in simple]
for i,w in enumerate(words):
 p=np.arange(126)
 for x in w:p=gens[x-1][p]
 assert s.fusion(p)==i,(i,s.fusion(p))
print('All 60 E7 class fusions agree with independent coordinate classword evaluation.')
