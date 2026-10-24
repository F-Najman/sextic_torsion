#!/usr/bin/env python3
"""Replay exported FI jets using independent Python GF(2) arithmetic.

The exporter computes the geometric input in Magma.  This reader checks
matrix dimensions/prefix consistency, enumerates every weighted divisor,
and verifies one representative of each diamond orbit.  It cannot replace
the geometric construction or its independent DifferentialSpace checks.
"""
import argparse
import ast
import json
from pathlib import Path
import re


def binary_rank(columns):
    basis = {}
    for column in columns:
        while column:
            pivot = column.bit_length()
            if pivot not in basis:
                basis[pivot] = column
                break
            column ^= basis[pivot]
    return len(basis)


def weighted_vectors(weights, degree):
    def extend(index, remaining, prefix):
        if index == len(weights):
            if remaining == 0:
                yield tuple(prefix)
            return
        for coefficient in range(remaining // weights[index] + 1):
            prefix.append(coefficient)
            yield from extend(index + 1, remaining - coefficient * weights[index], prefix)
            prefix.pop()
    yield from extend(0, degree, [])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jet_file", type=Path)
    parser.add_argument("--limit", type=int, default=0,
                        help="Stop after this many orbit representatives; zero is unlimited.")
    parser.add_argument("--stop-after-exceptions", type=int, default=0)
    args = parser.parse_args()
    source = args.jet_file.read_text()
    assert source.count('print "FI_JET_DATA_READY";') == 1, "Incomplete or appended export"
    values = {}
    for name in ("N", "D", "gx", "o", "udim", "ds"):
        matches = re.findall(rf"\b{name}:=(.*?);", source)
        assert len(matches) == 1, (name, "missing or repeated header")
        values[name] = ast.literal_eval(matches[0].strip())
    degree, order, dimension = values["D"], values["o"], values["udim"]
    factor_match = re.search(r"\bJetFactor:=(\d+);", source)
    jet_factor = int(factor_match[1]) if factor_match else 2
    weights = values["ds"]
    step_match = re.search(r"\bCoefficientSteps:=(.*?);", source, re.S)
    steps = ast.literal_eval(step_match[1]) if step_match else [1]*len(weights)
    assert len(steps) == len(weights) and all(s >= 1 for s in steps)
    if re.search(r"\bAs:=", source):
        values["As"] = ast.literal_eval(re.search(r"\bAs:=(.*?);", source)[1])
        permutations = ast.literal_eval(re.search(r"\bperms:=(.*?);", source, re.S)[1])
    else:
        values["As"] = [ast.literal_eval(re.search(r"\ba:=(.*?);", source)[1])]
        permutations = [ast.literal_eval(re.search(r"\bperm:=(.*?);", source, re.S)[1])]
    permutations = [tuple(i - 1 for i in permutation) for permutation in permutations]
    for permutation in permutations:
        assert sorted(permutation) == list(range(len(weights)))
        assert all(weights[i] == weights[permutation[i]] for i in range(len(weights)))
        assert all(steps[i] == steps[permutation[i]] for i in range(len(weights)))
    identity = tuple(range(len(weights)))
    group, queue = {identity}, [identity]
    while queue:
        power = queue.pop()
        for permutation in permutations:
            product = tuple(permutation[i] for i in power)
            if product not in group:
                group.add(product)
                queue.append(product)
                assert len(group) <= order
    assert len(group) == order

    jets = {}
    pattern = r"localjets\[<(\d+),(\d+)>\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);"
    for match in re.finditer(pattern, source, flags=re.S):
        index, multiplicity, rows, columns = map(int, match.group(1, 2, 3, 4))
        index -= 1
        key = index, multiplicity
        assert key not in jets and rows == dimension
        assert columns == jet_factor * multiplicity * weights[index]
        entries = ast.literal_eval(match[5])
        assert len(entries) == rows * columns and set(entries) <= {0, 1}
        jets[key] = tuple(sum(entries[row * columns + column] << row
                              for row in range(rows)) for column in range(columns))
    for index, weight in enumerate(weights):
        previous = ()
        for multiplicity in range(1, degree // weight + 1):
            current = jets[index, multiplicity]
            assert current[:len(previous)] == previous
            previous = current
    assert len(jets) == sum(degree // weight for weight in weights)

    def move(vector, permutation):
        result = [0] * len(weights)
        for index, value in enumerate(vector):
            result[permutation[index]] = value
        return tuple(result)

    tested = covered = enumerated = 0
    exceptions = []
    complete = True
    for raw_vector in weighted_vectors([w*s for w,s in zip(weights,steps)], degree):
        vector = tuple(v*s for v,s in zip(raw_vector,steps))
        enumerated += 1
        orbit = {move(vector, permutation) for permutation in group}
        if vector != min(orbit):
            continue
        if args.limit and tested >= args.limit:
            complete = False
            break
        columns = [column for index, multiplicity in enumerate(vector) if multiplicity
                   for column in jets[index, multiplicity]]
        assert len(columns) == jet_factor * degree
        rank = binary_rank(columns)
        tested += 1
        covered += len(orbit)
        if rank < jet_factor * degree:
            item = dict(index=tested, rank=rank, vector=vector, orbit_size=len(orbit))
            exceptions.append(item)
            print("FI_REPLAY_EXCEPTION " + json.dumps(item), flush=True)
            if args.stop_after_exceptions and len(exceptions) >= args.stop_after_exceptions:
                complete = False
                break
    if complete:
        assert covered == enumerated
        if values["N"] == 57 and degree == 6:
            assert (tested, covered) == (9104 if order == 18 else 18080, 162672)
        if values["N"] == 65 and degree == 6:
            if values["gx"] == 41:
                assert (tested, covered) == (1792,1792)
            elif values["gx"] == 21:
                assert (tested, covered) == ((31,106) if step_match else (161,560))
            else:
                assert (tested, covered) == (19921 if order == 24 else 39779, 476824)
    print("FI_REPLAY_DONE " + json.dumps(dict(level=values["N"], diamonds=values["As"],
          degree=degree, jet_factor=jet_factor, representatives=tested, covered=covered, enumerated=enumerated,
          exceptions=len(exceptions), complete=complete)), flush=True)


if __name__ == "__main__":
    main()
