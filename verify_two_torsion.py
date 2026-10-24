#!/usr/bin/env python3
"""Independent GF(2) replay of exported integral Betti operators.

This checks the finite matrix certificate, not Magma's construction of
JOne or its integral homology lattice. The script that constructs those
operators explicitly uses JOne, sign0, and IntegralMatrix on one lattice.
"""
import ast
from pathlib import Path
import re
import sys


def matrix_rows(source, name):
    found = re.findall(rf"\b{name}:=Matrix\(Integers\(\),(\d+),(\d+),(\[.*?\])\);",
                       source, re.S)
    assert len(found) == 1, (name, "missing or repeated matrix")
    rows, columns = map(int, found[0][:2])
    entries = ast.literal_eval(found[0][2])
    assert len(entries) == rows * columns
    return [sum((entries[i * columns + j] % 2) << j for j in range(columns))
            for i in range(rows)]


def left_kernel(rows):
    # Preserve the linear combination of the original rows during elimination.
    # Every dependency gives a row vector v with v*A=0.
    pivots, dependencies = {}, []
    for index, row in enumerate(rows):
        combination = 1 << index
        while row:
            pivot = row.bit_length()
            if pivot not in pivots:
                pivots[pivot] = row, combination
                break
            other, other_combination = pivots[pivot]
            row ^= other
            combination ^= other_combination
        if not row:
            dependencies.append(combination)
    assert len(pivots) + len(dependencies) == len(rows)
    return dependencies


def apply(vector, rows):
    result = 0
    while vector:
        bit = vector & -vector
        result ^= rows[bit.bit_length() - 1]
        vector ^= bit
    return result


source = Path(sys.argv[1]).read_text()
level = int(re.search(r"\bN:=(\d+);", source)[1])
dimension = int(re.search(r"\bd:=(\d+);", source)[1])
diamond = matrix_rows(source, "D2Z")
assert len(diamond) == dimension
operators = re.findall(r"\b(Aq\d+Z):=", source)
assert operators
kernel = None
for name in operators:
    operator = matrix_rows(source, name)
    if kernel is None:
        kernel = left_kernel(operator)
    else:
        images = [apply(vector, operator) for vector in kernel]
        combinations = left_kernel(images)
        kernel = [apply(vector, kernel) for vector in combinations]
    assert all(apply(vector, operator) == 0 for vector in kernel)
    assert all(apply(vector, diamond) == vector for vector in kernel), "Diamond is nontrivial"
if level == 57:
    assert dimension == 170 and len(kernel) == 2
    for index in range(dimension):
        vector = 1 << index
        for _ in range(9):
            vector = apply(vector, diamond)
        assert vector == 1 << index
print(dict(level=level, integral_homology_rank=dimension,
           kernel_dimension=len(kernel), diamond_two_trivial=True,
           independent_binary_replay=True))
