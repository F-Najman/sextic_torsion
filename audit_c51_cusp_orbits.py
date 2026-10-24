"""Primitive-column check of the mixed-sign Y51 quotient cusp signature."""
from collections import Counter
from math import gcd,comb
import json
N=51
def canon(x,y):
    y%=N;d=gcd(y,N)
    return min((x%d,y),((-x)%d,(-y)%N))
cs={canon(x,y) for y in range(N) for x in range(gcd(y,N))
    if gcd(x,gcd(y,N))==1}
def f(c):return canon(2*c[0],c[1])
def h(c):return canon(16*c[0],pow(16,-1,N)*c[1])
def orbit(x,fn):
    out=set()
    while x not in out:out.add(x);x=fn(x)
    return frozenset(out)
source={orbit(c,f) for c in cs}
hs={orbit(c,h) for c in cs}
assert len(cs)==64 and len(hs)==32 and {len(o) for o in hs}=={2}
def fy(o):return frozenset(f(c) for c in o)
target={orbit(o,fy) for o in hs}
assert Counter(map(len,source))=={1:16,2:8,8:4}
assert Counter(map(len,target))=={1:16,4:4}
total=comb(21,6)+4*comb(17,2)
assert total==54808
print(json.dumps(dict(source_signature=sorted(Counter(map(len,source)).items()),
                      quotient_signature=sorted(Counter(map(len,target)).items()),
                      all_degree6_cusp_divisors=total,quotient_etale=True,
                      exact=True),indent=2))
