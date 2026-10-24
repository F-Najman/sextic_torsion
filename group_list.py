#!/usr/bin/env python3
"""Lemma 3.1: the minimal groups G = Z/m x Z/n (m | n) outside B with phi(m) | 6 and every
prime divisor of n in S(6).  Prints the list and its size, and checks that it has 66 elements.
Run: python3 group_list.py      (closing marker GROUPLIST_DONE)"""
from math import gcd

S6 = {2, 3, 5, 7, 11, 13, 17, 19, 37}

def phi(m):
    return sum(1 for k in range(1, m + 1) if gcd(k, m) == 1)

def primes(n):
    out, p = set(), 2
    while p * p <= n:
        while n % p == 0:
            out.add(p); n //= p
        p += 1
    if n > 1:
        out.add(n)
    return out

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

PHI_INF = {(1, n) for n in range(1, 31) if n not in (23, 25, 29)}
PHI_INF |= {(2, 2 * n) for n in range(1, 11)} | {(3, 3 * n) for n in range(1, 5)}
PHI_INF |= {(4, 4), (4, 8), (6, 6)}
assert len(PHI_INF) == 44
B = PHI_INF | {(1, 25), (1, 37), (4, 12), (2, 24)}

def proper_subgroups(m, n):
    return [(a, b) for a in divisors(m) for b in divisors(n) if b % a == 0 and (a, b) != (m, n)]

def minimal_outside(Bset):
    out = []
    for m in [k for k in range(1, 19) if 6 % phi(k) == 0]:
        # if m > 1 then Z/n is a proper subgroup, so n <= 37; for m = 1, n <= 37^2
        nmax = 37 * 37 if m == 1 else 37
        for n in range(m, nmax + 1, m):
            if (m, n) in Bset or not primes(n) <= S6:
                continue
            if all(s in Bset for s in proper_subgroups(m, n)):
                out.append((m, n))
    return out

def name(g):
    m, n = g
    return f"Z/{n}" if m == 1 else f"Z/{m} x Z/{n}"

if __name__ == "__main__":
    G = minimal_outside(B)
    print("closure of B under subgroups:", all(s in B for g in B for s in proper_subgroups(*g)))
    print("m values allowed:", [k for k in range(1, 19) if 6 % phi(k) == 0])
    print("minimal groups outside B:", len(G))
    print("cyclic:", [n for m, n in G if m == 1])
    print("non-cyclic:", [name(g) for g in G if g[0] > 1])
    assert len(G) == 66
    print("GROUPLIST_DONE")
