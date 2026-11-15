"""F4 cells reconstructed from E6 roots; no stored matrices are read.

Extracted from the project Round15 integer audit. See PROVENANCE.txt.
"""
import itertools
import numpy as np

def roots_data():
    C=2*np.eye(6,dtype=np.int64)
    for a,b in [(0,2),(2,3),(3,4),(4,5),(1,3)]: C[a,b]=C[b,a]=-1
    rr=set(map(tuple,np.eye(6,dtype=np.int64))); todo=list(rr)
    while todo:
        r=np.array(todo.pop())
        for a in range(6):
            q=r.copy();q[a]-=(C@r)[a];q=tuple(q)
            if q not in rr:rr.add(q);todo.append(q)
    roots=sorted(r for r in rr if all(x>=0 for x in r));R=np.array(roots);G=R@C@R.T
    lookup={r:i for i,r in enumerate(roots)}
    def line(v):
        v=np.array(v);sg=1 if next(x for x in v if x)>0 else -1
        return lookup[tuple(sg*v)],sg
    refl=np.array([[line(R[b]-G[a,b]*R[a])[0] for b in range(36)]for a in range(36)])
    def rset(a,X):return tuple(sorted(int(refl[a,b]) for b in X))
    def rvec(a,v):return v-int(R[a]@C@v)*R[a]
    idx=[lookup[tuple(v)]for v in np.eye(6,dtype=np.int64)]
    sigma=lambda a:lookup[tuple(R[a][[5,1,4,3,2,0]])]
    fixed=tuple(a for a in range(36) if sigma(a)==a)
    pairs=sorted({tuple(sorted((a,sigma(a))))for a in range(36) if sigma(a)!=a})
    pairs+=sorted(X for X in itertools.combinations(fixed,2) if G[X]==0)
    def closure(B):
        B=set(B);changed=True
        while changed:
            changed=False
            for triple in itertools.combinations(sorted(B),3):
                for a in range(36):
                    if all(abs(G[a,b])==1 for b in triple):
                        z=a
                        for b in triple:z=int(refl[b,z])
                        z=int(refl[a,z])
                        if z not in B:B.add(z);changed=True
        return tuple(sorted(B))
    def actE(a,B):
        if a in B:return B
        bs=[b for b in B if G[a,b]]
        return rset(bs[0],rset(a,B)) if bs else closure(B+(a,))
    words=[[idx[1]],[idx[3]],[idx[2],idx[4]],[idx[0],idx[5]]]
    seen={()};todo=[()]
    for B in todo:
        for word in words:
            for typ in ['r','e']:
                X=B
                for a in reversed(word):X=rset(a,X) if typ=='r' else actE(a,X)
                if X not in seen:seen.add(X);todo.append(X)
    fours=sorted(B for B in seen if len(B)==4)
    allpairs=[X for X in itertools.combinations(range(36),2) if G[X]==0]
    B=pairs[0];D0=tuple(Y for Y in allpairs if all(G[a,b]==0 for a in B for b in Y))
    def action(a,X,D,typ,oriented=False):
        power=0
        if typ=='e':
            if a in X:return X,D,1
            bs=[b for b in X if G[a,b]]
            if not bs:return None,None,0
            reflword=[a,bs[0]]
        else:reflword=[a]
        for b in reflword:
            X=rset(b,X)
            D=tuple(rvec(b,Z) for Z in D) if oriented else tuple(rset(b,Z) for Z in D)
        return X,D,power
    decor={B:D0};todo=[B]
    for X in todo:
        for word in words:
            for typ in ['r','e']:
                Y,D=X,decor[X]
                for a in reversed(word):
                    Y,D,_=action(a,Y,D,typ)
                    if Y is None:break
                if Y in pairs and Y not in decor:decor[Y]=D;todo.append(Y)
    assert len(roots)==36 and len(fixed)==12 and len(pairs)==30 and len(fours)==39
    return dict(C=C,R=R,G=G,refl=refl,rset=rset,rvec=rvec,line=line,idx=idx,
                fixed=fixed,pairs=pairs,fours=fours,words=words,action=action,decor=decor)

def pmat(p,kind):
    if kind=='triv':return np.ones((1,1),dtype=np.int64)
    if kind=='sign':return np.array([[(-1)**sum(p[a]>p[b]for a in range(3)for b in range(a+1,3))]],dtype=np.int64)
    P=np.eye(3,dtype=np.int64)[:,p]
    return (P[:,:2]-P[:,2:3])[:2,:]

SINGLE=[('3|0',0,'triv'),('111|0',0,'sign'),('21|0',0,'std'),
        ('0|3',3,'triv'),('0|111',3,'sign'),('0|21',3,'std'),
        ('2|1',1,'triv'),('11|1',1,'sign'),('1|2',2,'triv'),('1|11',2,'sign')]

def signed_rep(p,sg,k,kind):
    if k in [0,3]:return pmat(p,kind)*(int(np.prod(sg))if k==3 else 1)
    subsets=list(itertools.combinations(range(3),k));M=np.zeros((3,3),dtype=np.int64)
    for j,B in enumerate(subsets):
        image=tuple(sorted(p[b]for b in B));i=subsets.index(image)
        coeff=int(np.prod([sg[b]for b in B]))
        if kind=='sign':
            H=[b for b in range(3)if (b in B)==(k==2)]
            coeff*=(-1)**sum(p[a]>p[b]for a,b in itertools.combinations(H,2))
        M[i,j]=coeff
    return M

def make_cells(D):
    R,G,words=D['R'],D['G'],D['words'];cells=[]
    for layer,name,k,kind in [(1,*spec)for spec in SINGLE]+[(2,'triv',0,'triv'),(2,'sign',0,'sign'),(2,'std',0,'std'),(4,'triv',0,'triv')]:
        sets=[(a,)for a in D['fixed']]if layer==1 else D['pairs']if layer==2 else D['fours']
        if layer==1:
            dec={X:tuple(R[b].copy()for b in D['fixed']if G[X[0],b]==0)for X in sets}
            dim=signed_rep([0,1,2],[1,1,1],k,kind).shape[0]
        elif layer==2:dec=D['decor'];dim=2 if kind=='std'else 1
        else:dec={X:()for X in sets};dim=1
        n=len(sets)*dim;si={X:i for i,X in enumerate(sets)}
        def block(Y,deco):
            if layer==1:
                lines=[D['line'](Z)for Z in deco];dest=[D['line'](Z)[0]for Z in dec[Y]]
                return signed_rep([dest.index(a)for a,sg in lines],[sg for a,sg in lines],k,kind)
            if layer==2:return pmat([dec[Y].index(Z)for Z in deco],kind)
            return np.ones((1,1),dtype=np.int64)
        def apply(X,deco,word,typ):
            power=0
            for a in reversed(word):
                X,deco,l=D['action'](a,X,deco,typ,layer==1);power+=l
                if X is None:break
            return X,deco,power
        Rs=[];Es=[]
        for word in words:
            r=np.zeros((n,n),dtype=np.int64);ec=np.zeros((3,n,n),dtype=np.int64)
            for j,X in enumerate(sets):
                Y,dd,_=apply(X,dec[X],word,'r');i=si[Y]
                r[i*dim:(i+1)*dim,j*dim:(j+1)*dim]=block(Y,dd)
                Y,dd,power=apply(X,dec[X],word,'e')
                if Y is not None:
                    i=si[Y];ec[power,i*dim:(i+1)*dim,j*dim:(j+1)*dim]=block(Y,dd)
            Rs.append(r);Es.append(ec)
        H=np.kron(np.eye(len(sets),dtype=np.int64),np.array([[2,1],[1,2]],dtype=np.int64))if kind=='std'and dim==2 else np.eye(n,dtype=np.int64)
        gram=np.zeros((layer+1,n,n),dtype=np.int64)
        for i,Y in enumerate(sets):
            for j,X in enumerate(sets):
                Z,dd,power=apply(X,dec[X],list(Y),'e')
                if Z is not None:
                    assert Z==Y
                    gram[power,i*dim:(i+1)*dim,j*dim:(j+1)*dim]=block(Y,dd)
        gram=np.array([H@A for A in gram])
        assert all(np.array_equal(A,A.T)for A in gram),(layer,name,'Gram symmetry')
        cells.append(dict(label=f'{layer}:{name}',layer=layer,name=name,sets=sets,dim=n,
                          Rs=Rs,Es=Es,G=gram,H=H))
    return cells

def verify_relations(cell):
    # Liu order is reverse of our manuscript long-long-short-short order.
    R=cell['Rs'][::-1];EC=cell['Es'][::-1];I=np.eye(cell['dim'],dtype=np.int64);checks=0
    def eq(a,b,tag):
        nonlocal checks
        assert np.array_equal(a,b),(cell['label'],z,tag);checks+=1
    for z in [-2,-1,0,1,2,3,4]:
        E=[sum((z**k*A for k,A in enumerate(es)),np.zeros_like(es[0]))for es in EC]
        B=sum((z**k*A for k,A in enumerate(cell['G'])),np.zeros_like(cell['G'][0]))
        for i,(r,e)in enumerate(zip(R,E)):
            eq(r@r,I,('r2',i));eq(r@e,e,('re',i));eq(e@r,e,('er',i));eq(e@e,z**(2 if i<2 else 1)*e,('e2',i))
            eq(r.T@B,B@r,('rform',i));eq(e.T@B,B@e,('eform',i))
            for j,(rr,ee)in enumerate(zip(R,E)):
                if i==j:continue
                if abs(i-j)>1:
                    eq(r@rr,rr@r,('commrr',i,j));eq(e@rr,rr@e,('commer',i,j));eq(e@ee,ee@e,('commee',i,j))
                elif {i,j}!={1,2}:
                    eq(r@rr@r,rr@r@rr,('braid',i,j));eq(rr@r@ee,e@ee,('rree',i,j));eq(r@ee@r,rr@e@rr,('rer',i,j))
        r,q=R[1:3];e,f=E[1:3]
        for lhs,rhs,tag in [(r@q@r@q,q@r@q@r,'braid4'),(r@q@e,q@e,'rse'),(r@f@r@f,f@e@f,'rfrf'),(r@q@r@f,f@r@q@r,'rsrf'),(e@q@e,z*e,'ese'),(e@f@e,z*e,'efe'),(e@q@r,e@q,'esr'),(e@f@r,e@f,'efr')]:eq(lhs,rhs,tag)
    return checks

