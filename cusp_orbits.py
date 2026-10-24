from math import gcd
from collections import Counter

# Tate cusp labels: Q=(a,b) modulo15, P=(5c,5d), det(P/5,Q)=1 mod3.
# Changing q^(1/15) by a root of unity gives the upper unipotent action;
# simultaneous negation identifies isomorphic level structures.
labels={(a,b,c,d) for a in range(15) for b in range(15)
        if gcd(gcd(a,b),15)==1 for c in range(3) for d in range(3)
        if (c*b-d*a)%3==1}
def orbit(x):
    a,b,c,d=x
    return frozenset(((s*(a+k*b))%15,(s*b)%15,
                      (s*(c+k*d))%3,(s*d)%3)
                     for s in [-1,1] for k in range(15))
classes=set(map(orbit,labels)); classes=sorted(classes,key=lambda s:sorted(s))
index={x:i for i,cs in enumerate(classes) for x in cs}
def gal(t,x):
    a,b,c,d=x
    return ((t*a)%15,b,(t*c)%3,d)
def cycles(t):
    perm=[index[gal(t,next(iter(cs)))] for cs in classes]
    left=set(range(len(classes)));result=[]
    while left:
        i=min(left);cyc=[]
        while i not in cyc:
            cyc.append(i);left.remove(i);i=perm[i]
        result.append(cyc)
    return result
print('cusp count',len(classes),'widths',Counter(len(cs)//2 for cs in classes))
for t in [7,4]:
    cy=cycles(t)
    print('Galois action',t,'orbit counts',dict(Counter(map(len,cy))))
    print('orbits',cy)
print('CUSP_ORBITS_3_15_DONE')
