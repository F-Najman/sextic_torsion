"""Exact primitive-column cusp/Frobenius audit of X1(65)/<4>.

For d=gcd(y,65), a Gamma1 cusp is (x mod d,y mod65), x a unit mod d,
modulo simultaneous sign. The cyclotomic Frobenius2 acts on x by2;
diamond a acts as (a*x,a^-1*y). This is the usual integral cusp
description. This checks the finite permutations independently of Magma.
"""
from collections import Counter
from math import gcd, comb
import json

N = 65

def canon(x, y):
    y %= N
    d = gcd(y, N)
    return min((x % d, y), ((-x) % d, (-y) % N))

cusps = {canon(x, y) for y in range(N)
         for x in range(gcd(y,N)) if gcd(x,gcd(y,N)) == 1}
assert len(cusps) == 96

def diamond(c, a):
    x,y = c
    return canon(a*x, pow(a,-1,N)*y)

def frob(c):
    x,y = c
    return canon(2*x,y)

def orbit(x, f):
    result = set()
    while x not in result:
        result.add(x)
        x = f(x)
    return frozenset(result)

source_orbits = {orbit(c,frob) for c in cusps}
assert Counter(map(len,source_orbits)) == {1:24,4:6,6:4,12:2}
h_orbits = {orbit(c,lambda x:diamond(x,4)) for c in cusps}
assert len(h_orbits) == 32 and {len(o) for o in h_orbits} == {3}

def quotient_frob(o):
    return frozenset(frob(c) for c in o)

target_orbits = {orbit(o,quotient_frob) for o in h_orbits}
assert Counter(map(len,target_orbits)) == {1:8,2:4,4:4}
source_to_target = {c:o for o in h_orbits for c in o}
rows=[]
for target in target_orbits:
    source = [s for s in source_orbits
              if source_to_target[next(iter(s))] in target]
    ds = {gcd(c[1],N) for s in source for c in s}
    assert len(ds) == 1
    d = ds.pop()
    degrees = {len(s) for s in source}
    assert len(degrees) == 1
    degree = degrees.pop()
    assert degree % len(target) == 0
    rows.append((d,degree,len(target),degree//len(target),len(source)))

counts=Counter(rows)
assert counts == {(1,1,1,1,3):8,(5,4,4,1,3):2,
                  (13,12,4,3,1):2,(65,6,2,3,1):4}
admissible=comb(13,6)+2*comb(9,2)+4
assert admissible == 1792

# The two explicit AS lifts lie above independent modular diamonds;
# their common square is the relative-sign involution.
def target_diamond(o,a):
    return frozenset(diamond(c,a) for c in o)

assert all(target_diamond(target_diamond(o,27),27) ==
           target_diamond(target_diamond(o,41),41) for o in h_orbits)
assert any(target_diamond(target_diamond(o,27),27) != o for o in h_orbits)
assert all(target_diamond(target_diamond(target_diamond(target_diamond(o,a),a),a),a)==o
           for o in h_orbits for a in (27,41))
perms=set()
for a in range(4):
    for b in range(4):
        images=tuple(sorted((tuple(sorted(o)),
                            tuple(sorted(target_diamond(o,pow(27,a,N)*pow(41,b,N)%N))))
                           for o in h_orbits))
        perms.add(images)
assert len(perms)==8

print(json.dumps(dict(source_signature=sorted(Counter(map(len,source_orbits)).items()),
                      target_signature=sorted(Counter(map(len,target_orbits)).items()),
                      cusp_types=[dict(d=k[0],source_degree=k[1],target_degree=k[2],
                                       coefficient_step=k[3],source_places=k[4],
                                       target_places=v) for k,v in sorted(counts.items())],
                      admissible_degree6_divisors=admissible,
                      residual_diamond_group_order=len(perms),exact=True),indent=2))
