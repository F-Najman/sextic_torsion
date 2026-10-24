#!/usr/bin/env python3
"""Independent GF(2) replay of the exact integral JOne49 certificate."""
import ast
from pathlib import Path
import re

source = (Path(__file__).parent / "data/two_torsion_49_T2_matrices.m").read_text()
def matrix(name):
    matches = re.findall(rf"\b{name}:=Matrix\(Integers\(\),(\d+),(\d+),(\[.*?\])\);", source, re.S)
    assert len(matches) == 1
    nr, nc = map(int, matches[0][:2])
    es = ast.literal_eval(matches[0][2]); assert len(es) == nr*nc
    return [sum((es[nc*i+j] % 2) << j for j in range(nc)) for i in range(nr)]
def kernel(rows):
    pivots, deps = {}, []
    for i, row in enumerate(rows):
        word = 1 << i
        while row:
            k = row.bit_length()
            if k not in pivots:
                pivots[k] = row, word
                break
            p, q = pivots[k]; row ^= p; word ^= q
        if not row: deps.append(word)
    return deps
def act(v, rows):
    out = 0
    while v:
        b = v & -v; out ^= rows[b.bit_length()-1]; v ^= b
    return out
A5, T2 = matrix("A5Z"), matrix("T2Z")
assert len(A5) == len(T2) == 138
K = kernel(A5)
assert len(K) == 2 and all(act(v,A5) == 0 for v in K)
images = [act(v,T2) for v in K]
assert len(kernel(images)) == 0
root_file=Path(__file__).parent / "data/two_torsion_49_T2_root_matrices.m"
if root_file.exists():
    root_source=root_file.read_text()
    for name in ("A5Z","T2Z"):
        pattern=rf"\b{name}:=Matrix\(Integers\(\),138,138,(\[.*?\])\);"
        assert re.findall(pattern,source,re.S)==re.findall(pattern,root_source,re.S)
    print("TWO49_ROOT_EXPORTED_MATRICES_EQUAL true")
print({"N":49,"integral_homology_rank":138,"kernel_A5":2,
       "kernel_A5_T2":0,"independent_binary_replay":True})
