#!/usr/bin/env python3
"""Check rank-zero Cartier image and every exported jet projection over F2.

This replays exact finite linear algebra independently of Magma. It does
not reconstruct the function field, the Cartier operator, or the local
expansions. replay_fi_jets.py separately enumerates every cusp divisor.
"""
import ast
from pathlib import Path
import re
import sys

source = Path(sys.argv[1]).read_text()
assert 'FI_JET_DATA_READY' in source
def integer(name):
    ms = re.findall(rf'\b{name}:=(\d+);', source)
    assert len(ms) == 1, name
    return int(ms[0])
def matrix(pattern):
    ms = re.findall(pattern + r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);', source, re.S)
    assert len(ms) == 1, pattern
    nr, nc = map(int, ms[0][:2]); es = ast.literal_eval(ms[0][2])
    assert len(es) == nr*nc and all(e in (0,1) for e in es)
    return [sum(es[nc*i+j] << j for j in range(nc)) for i in range(nr)], nc
def echelon(rows):
    piv = {}
    for x in rows:
        while x:
            b = x.bit_length()-1
            if b in piv: x ^= piv[b]
            else:
                piv[b] = x
                break
    return piv
def contains(piv, x):
    while x:
        b = x.bit_length()-1
        if b not in piv: return False
        x ^= piv[b]
    return True
def apply(v, rows):
    out = 0
    while v:
        b = v & -v; out ^= rows[b.bit_length()-1]; v ^= b
    return out

N, D, g, u = (integer(k) for k in ('N','D','gx','udim'))
ds = ast.literal_eval(re.search(r'\bds:=(\[.*?\]);', source, re.S)[1])
CM, nc = matrix(r'\bCartierMatrix'); assert len(CM) == nc == g
U, nc = matrix(r'\bImageBasis'); assert len(U) == u and nc == g
im = [row ^ (1 << i) for i,row in enumerate(CM)]
if 'DifferenceImageBasis:=' in source:
    DU, nc = matrix(r'\bDifferenceImageBasis'); assert nc == g
    if 'MT:=' in source and 'MU:=' in source:
        MT, nt = matrix(r'\bMT'); MU, nu = matrix(r'\bMU')
        assert len(MT) == len(MU) == nt == nu == g
        expected_du = [row ^ (1 << i) for M in (MT,MU) for i,row in enumerate(M)]
        dp = echelon(expected_du)
        assert len(dp) == len(echelon(DU))
        assert all(contains(dp,row) for row in DU)
    im = [apply(v,im) for v in DU]
else:
    assert re.search(r'\bAs:=\[\];', source)
piv = echelon(im)
assert len(piv) == u == len(echelon(U))
assert all(contains(piv,row) for row in U)
checked = 0
for i, degree in enumerate(ds,1):
    B, nc = matrix(rf'\bbasisjets\[{i}\]')
    assert len(B) == g and nc == (D//degree)*degree
    projection = [apply(row,B) for row in U]
    for m in range(1,D//degree+1):
        J, width = matrix(rf'\blocaljets\[<{i},{m}>\]')
        assert len(J) == u and width == m*degree
        mask = (1 << width)-1
        assert J == [row & mask for row in projection], (i,m)
        checked += 1
print({'level':N,'genus':g,'cartier_image_dimension':u,
       'cusps':len(ds),'projected_jet_matrices':checked,
       'independent_projection_replay':True})
