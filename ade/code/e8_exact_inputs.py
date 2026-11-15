"""Recheck reconstructed E8 character products with Python unbounded integers."""
import e8_singleton as e
import e8_pair_branching as p
import numpy as np
for n in (7,8):
 c=getattr(e,'ch'+str(n)); cl=getattr(e,'cl'+str(n));
 x=(c.astype(object)*cl.astype(object)) @ c.T.astype(object)
 assert np.array_equal(x,(c*cl)@c.T)
 print('W', n, 'orthogonality exact max', max(abs(int(t)) for t in x.flat))
x=((e.ch8[:,e.fus0].astype(object)+e.ch8[:,e.fus1].astype(object))*e.cl7.astype(object))@e.ch7.T.astype(object)
assert np.array_equal(x,e.num.astype(object))
print('singleton induction exact max',max(abs(int(t)) for t in x.flat))
for th in p.labels:
 chi=np.array([p.symchar(th,tuple(sorted(c['mu']+c['nu'],reverse=True))) for c in p.classes],dtype=object)
 x=(e.ch8[:,[c['fusion'] for c in p.classes]].astype(object)*np.array([c['size'] for c in p.classes],dtype=object))@chi
 assert np.array_equal(x,184320*p.B[:,p.labels.index(th)].astype(object))
print('pair branching all Python integer products match')
