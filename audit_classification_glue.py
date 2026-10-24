#!/usr/bin/env python3
"""Independent finite-group glue audit for the current48-group lower bound.

Uses maximal index-prime subgroup types, independently of candidates.py's
full divisor-pair enumeration. Standard library only; read-only input.
The individual geometric exclusion certificates are audited separately.
"""
import argparse
import json
from pathlib import Path


def prime_divisors(n):
    p=2;out=[]
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out


def phi(n):
    result=n
    for p in prime_divisors(n):result=result//p*(p-1)
    return result


def divisors(n):
    return [d for d in range(1,n+1) if n%d==0]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('candidates_json',type=Path)
    a=p.parse_args()
    data=json.loads(a.candidates_json.read_text())
    primes={2,3,5,7,11,13,17,19,37}
    # If phi(m)|6 then every prime p|m has p-1|6, hence p=2,3,7.
    # Prime-power exponents are a<=2,b<=2,c<=1 respectively.
    ms=sorted({2**a*3**b*7**c for a in range(3) for b in range(3)
               for c in range(2) if 6%phi(2**a*3**b*7**c)==0})
    assert ms==[1,2,3,4,6,7,9,14,18]
    inf=({(1,n) for n in range(1,31) if n not in (23,25,29)}
         |{(2,2*n) for n in range(1,11)}
         |{(3,3*n) for n in range(1,5)}|{(4,4),(4,8),(6,6)})
    known=inf|{(1,25),(1,37),(4,12),(2,24)}
    assert len(inf)==44 and len(known)==48
    assert inf=={tuple(g) for g in data['phi_infinite']}
    assert known=={tuple(g) for g in data['known_lower_bound']}
    for m,n in known:
        assert all((x,y) in known for x in divisors(m) for y in divisors(n) if y%x==0)
    bound=max(n for m,n in known if m==1)*max(primes)
    assert bound==1369
    minimal=[]
    for m in ms:
        for n in range(m,bound+1,m):
            if (m,n) in known or not set(prime_divisors(n))<=primes:continue
            maximal=[]
            for ell in prime_divisors(n):
                if m%ell==0:maximal.append((m//ell,n))
                if (n//ell)%m==0:maximal.append((m,n//ell))
            if all(g in known for g in maximal):minimal.append((m,n))
    assert len(minimal)==66
    assert sorted(minimal)==[tuple(g) for g in data['minimal_outside_known']]
    after37=[g for g in minimal if g[1]%37]
    assert len(after37)==57
    assert sorted(after37)==[tuple(g) for g in data['after_37_classification']]
    excluded=set()
    for k in ('reduction_exclusions','noncyclic_exclusions','sieve_exclusions',
              'isogeny_exclusions','operator_exclusions','divisor_class_exclusions'):
        excluded.update(tuple(g) for g in data[k])
    global_ns={g['n'] for g in data['global_exclusions']}
    assert global_ns=={119,121,125,143,169,187,209,221,289,323,361}
    for row in data['global_exclusions']:
        n,aa=row['n'],row['a'];a0=aa%n;power=a0;order=1
        while power not in (1,n-1):power=power*a0%n;order+=1
        assert order==row['diamond_order'] and 6%order!=0
        assert 325*row['index']>6*32768*row['k3']
        assert n%3 and row['cm_minimum']>6
    excluded.update((1,n) for n in global_ns)
    remaining=sorted(g for g in after37 if g not in excluded)
    assert remaining==[tuple(g) for g in data['remaining_after_verified_exclusions']]
    print(json.dumps({'independent_maximal_subgroup_enumeration':True,
                      'infinite_groups':len(inf),'known_groups':len(known),
                      'minimal_obstructions':len(minimal),'after_37':len(after37),
                      'inherited_global_rows':len(global_ns),'remaining':remaining},sort_keys=True))


if __name__=='__main__':main()
