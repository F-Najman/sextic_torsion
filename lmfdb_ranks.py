#!/usr/bin/env python3
"""The analytic ranks and character values used in Sections 4, 8 and 9.

For every level n below, we list the weight 2 newforms of level dividing n with positive
analytic rank, check that their ranks are recorded in the LMFDB as proved, and evaluate their
characters (given by Conrey labels) at the elements a used in the paper.

Run: python3 lmfdb_ranks.py           (offline, uses data/lmfdb_newforms_cache.json)
     python3 lmfdb_ranks.py --fetch   (refetch from the LMFDB mirror and rewrite the cache)
Closing marker: RANKS_DONE OK"""
import json, sys
from fractions import Fraction as Fr
from math import gcd

CACHE = "data/lmfdb_newforms_cache.json"
LEVELS = [119, 121, 125, 143, 169, 187, 209, 221, 289, 323, 361,
          77, 85, 91, 95, 133, 247, 51, 57, 63, 65, 55, 75]

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def fetch():
    import psycopg2
    conn = psycopg2.connect(host="devmirror.lmfdb.xyz", dbname="lmfdb", user="lmfdb",
                            password="lmfdb", connect_timeout=60)
    cur = conn.cursor()
    levels = sorted({d for n in LEVELS for d in divisors(n) if d >= 11})
    cur.execute("""select label, level, char_orbit_label, conrey_index, char_order, dim,
                   analytic_rank, analytic_rank_proved from mf_newforms
                   where weight=2 and level = ANY(%s) order by level, label""", (levels,))
    rows = [dict(label=r[0], level=r[1], char_orbit=r[2], conrey=int(r[3]), char_order=int(r[4]),
                 dim=int(r[5]), analytic_rank=(None if r[6] is None else int(r[6])),
                 rank_proved=bool(r[7])) for r in cur.fetchall()]
    with open(CACHE, "w") as f:
        json.dump(rows, f, indent=1)
    print("fetched", len(rows), "newforms")

# Conrey characters for odd n: chi_n(c, .) = product over p^e || n of the characters given by
# the smallest primitive root modulo p^e.
def primitive_root(q, p):
    phi = q - q // p
    fac = [r for r in range(2, phi + 1) if phi % r == 0 and all(r % s for s in range(2, r))]
    for g in range(2, q):
        if gcd(g, q) == 1 and all(pow(g, phi // r, q) != 1 for r in fac):
            return g

def dlog(g, x, q):
    y, k = 1, 0
    while y != x % q:
        y = y * g % q; k += 1
    return k

def factor(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def conrey_value(n, c, a):
    """chi_n(c, a) = exp(2 pi i t), returns t in [0,1)."""
    t = Fr(0)
    for p, e in factor(n).items():
        q = p ** e; phi = q - q // p
        g = primitive_root(q, p)
        t += Fr(dlog(g, c, q) * dlog(g, a, q), phi)
    return t % 1

# (level n, elements a that must lie in the kernel of every positive-rank character)
KERNEL = {119: [3], 121: [2], 125: [3], 143: [67], 169: [3], 187: [122], 209: [34], 221: [3],
          289: [3], 323: [101], 361: [3], 77: [3], 85: [2], 91: [3], 247: [12, 56],
          57: [2, 8], 65: [8, 4]}

def main():
    if "--fetch" in sys.argv:
        fetch()
    rows = json.load(open(CACHE))
    ok = True
    assert all(r["analytic_rank"] is not None for r in rows)
    for n in LEVELS:
        pos = [r for r in rows if n % r["level"] == 0 and r["analytic_rank"] > 0]
        print(f"n={n}: positive rank:", [(r["label"], r["conrey"], r["char_order"]) for r in pos])
        ok &= all(r["rank_proved"] for r in pos)
        for a in KERNEL.get(n, []):
            vals = [conrey_value(r["level"], r["conrey"], a % r["level"]) for r in pos]
            print(f"   values at {a}:", [str(v) for v in vals])
            ok &= all(v == 0 for v in vals)
    # special values used in the text
    special = [(95, 8, 31, Fr(1, 2)), (133, 39, 6, Fr(0)), (133, 31, 3, Fr(0)),
               (133, 39, 3, Fr(2, 3)), (133, 31, 6, Fr(1, 6)), (63, 47, 8, Fr(1, 2))]
    for n, c, a, expected in special:
        v = conrey_value(n, c, a)
        print(f"chi_{n}({c},{a}) = exp(2 pi i * {v})")
        ok &= v == expected
    pos63 = [r["label"] for r in rows if 63 % r["level"] == 0 and r["analytic_rank"] > 0]
    pos95 = [r["label"] for r in rows if 95 % r["level"] == 0 and r["analytic_rank"] > 0]
    pos57 = [r["label"] for r in rows if 57 % r["level"] == 0 and r["analytic_rank"] > 0]
    pos65 = [r["label"] for r in rows if 65 % r["level"] == 0 and r["analytic_rank"] > 0]
    pos51 = [r["label"] for r in rows if 51 % r["level"] == 0 and r["analytic_rank"] > 0]
    ok &= pos63 == ["63.2.s.a"] and pos95 == ["95.2.l.a"] and pos57 == ["57.2.a.a"]
    ok &= pos65 == ["65.2.a.a"] and pos51 == []
    print("RANKS_DONE", "OK" if ok else "FAILED")

if __name__ == "__main__":
    main()
