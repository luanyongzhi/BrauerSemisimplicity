"""Formal E8 pair relations and root conjugations, in scratch.
The geometry prefix is shared with e8_annihilator.py.
"""
from pathlib import Path
import json
source=Path("e8_annihilator.py").read_text()
exec(compile(source.split("B1_cols =",1)[0],"e8_annihilator.py (geometry prefix)","exec"),globals())
# Compare BOTH choices of bridge when a meets both roots in the pair.
bridge_cases = bridge_ray_checks = 0
for a in range(120):
    ks = np.flatnonzero((gram_120[a,pairs_i] != 0) & (gram_120[a,pairs_j] != 0) & (pairs_i != a) & (pairs_j != a))
    rd = r_dest[a,ks]
    b, c = pairs_i[ks], pairs_j[ks]
    d1, d2 = r_dest[b,rd], r_dest[c,rd]
    p1 = s6_mul_id[r_perm[b,rd], r_perm[a,ks]]
    p2 = s6_mul_id[r_perm[c,rd], r_perm[a,ks]]
    assert np.array_equal(d1,d2), ('bridge destinations', a)
    assert np.array_equal(p1,p2), ('bridge permutations', a)
    assert np.array_equal(d1,out[a,ks]) and np.array_equal(p1,e_perm[a,ks])
    bridge_cases += len(ks)
    bridge_ray_checks += 6*len(ks)

N=3780*6
z=np.zeros(N,dtype=np.int64)
identity=(np.arange(N,dtype=np.int64),z)
fiber=np.array(s6_list,dtype=np.int64)
def operator(dest,perm,power):
    d=np.where(dest[:,None]<0,-1,6*dest[:,None]+fiber[perm]).reshape(-1)
    p=np.repeat(power,6)
    p[d<0]=0
    return d,p
def mul(A,B):
    ad,ap=A; bd,bp=B
    d=np.full(N,-1,dtype=np.int64); p=z.copy()
    v=np.flatnonzero(bd>=0); v=v[ad[bd[v]]>=0]
    d[v]=ad[bd[v]]; p[v]=bp[v]+ap[bd[v]]
    return d,p
def prod(*ops):
    x=identity
    for op in ops: x=mul(x,op)
    return x
def check(A,B,label):
    assert np.array_equal(A[0],B[0]) and np.array_equal(A[1],B[1]), label
R=[operator(r_dest[a],r_perm[a],np.zeros(3780,dtype=np.int64)) for a in simple_indices]
E=[operator(out[a],e_perm[a],ex[a]) for a in range(120)]
counts={'single_node':0,'nonadjacent':0,'adjacent':0}
def ck(a,b,kind,label):
    check(a,b,label); counts[kind]+=1
for i,a in enumerate(simple_indices):
    r,e=R[i],E[a]
    ck(prod(r,r),identity,'single_node',(i,'r2'))
    ck(prod(r,e),e,'single_node',(i,'re'))
    ck(prod(e,r),e,'single_node',(i,'er'))
    ck(prod(e,e),(e[0],np.where(e[0]<0,0,e[1]+1)),'single_node',(i,'e2'))
    for j in range(i):
        rr,ee=R[j],E[simple_indices[j]]
        if gram_120[a,simple_indices[j]]==0:
            for name,x,y in [('rr',r,rr),('ee',e,ee),('re',r,ee),('er',e,rr)]:
                ck(prod(x,y),prod(y,x),'nonadjacent',(i,j,name))
        else:
            ck(prod(r,rr,r),prod(rr,r,rr),'adjacent',(i,j,'braid'))
            ck(prod(rr,e,rr),prod(r,ee,r),'adjacent',(i,j,'conjugation'))
            ck(prod(rr,r,ee),prod(e,ee),'adjacent',(i,j,'mixed1'))
            ck(prod(r,rr,e),prod(ee,e),'adjacent',(i,j,'mixed2'))
            ck(prod(ee,r,rr),prod(ee,e),'adjacent',(i,j,'mixed1 reversed'))
            ck(prod(e,rr,r),prod(e,ee),'adjacent',(j,i,'mixed2 reversed'))
conjugations=0
for i,a in enumerate(simple_indices):
    for b in range(120):
        check(prod(R[i],E[b],R[i]),E[int(refl_120[a,b])],('root conjugation',i,b))
        conjugations+=1
result={'basis_dimension':N,'pairs':3780,'rays_per_pair':6,'bridge_cases':bridge_cases,'bridge_ray_agreements':bridge_ray_checks,'cfw_relation_counts':counts,'cfw_relations_total':sum(counts.values()),'root_conjugation_identities':conjugations,'all_passed':True,'note':'Operators compared as partial monomial maps with formal delta exponents; the six-ray permutation action need not detect all elements of S6 in each individual comparison, but is faithful to S6 as a group.'}
path=Path('e8_relations.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
