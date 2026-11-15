"""Final exact E8 candidate screen after certified branching computations."""
import json
from fractions import Fraction
s=json.load(open('e8_singleton_results.json'));p=json.load(open('e8_pair_branching_results.json'));q=json.load(open('e8_four_branching_results.json'));e7=json.load(open('e7_singleton_results.json'))
sigma=set(s['singleton_bad_set']);sigma7=set(e7['exceptional_set'])
c={0:set(s['reflection_scalars_E8']),1:set(s['reflection_scalars_E7']),2:{-30,-18,-10,-6,0,6,10,18,30},4:{-12,0,12},8:{0}}
rows=[]
for t in [2,4,8]:
 for u in sorted(x for x in c if x<t):
  if (t,u)==(2,0):v=set(p['pair_s0_candidates'])
  elif (t,u)==(4,0):v=set(sum(q['s0_candidates'],[]))
  else:
   v=set()
   for a in c[t]:
    for b in c[u]:
     d=1+Fraction(a-b,t-u)
     if d and d.denominator==1:
      if (t,u)==(2,1) and int(d) not in sigma7:continue
      v.add(int(d))
  rows.append({'t':t,'s':u,'all_candidates':sorted(v),'outside_singleton':sorted(v-sigma)})
extra=set().union(*(set(r['outside_singleton']) for r in rows))
assert extra=={-3,5}
out={'singleton_bad_set':sorted(sigma),'higher_collision_rows':rows,'unresolved_values':sorted(extra),'proved_nonsemisimple_count':len(sigma),'containing_set':sorted(sigma|extra)}
json.dump(out,open('e8_final_screen_results.json','w'),indent=2)
print(json.dumps(rows,indent=2));print('unresolved',sorted(extra))
