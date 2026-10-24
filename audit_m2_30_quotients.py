"""Check coverage of the completed first and other-two quotient filters."""
from pathlib import Path
import re
import sys

base = Path(__file__).resolve().parent
prefix = (base / "logs/n30/first_quotient_prefix.log").read_text(errors="replace")
mixed = (base / "logs/n30/first_quotient_mixed_continuation.log").read_text(errors="replace")
# One parser error occurred between the completed single-place loop and
# the unstarted mixed loops. It is retained explicitly, rather than
# treating the initial process as a completed full run.
assert prefix.count("User error: bad syntax") == 1
assert prefix.rstrip().endswith("User error: bad syntax\n>")
assert not re.search(r"(?:Runtime error|Assertion failed|Error:)", prefix)
assert prefix.index("5006 18 2704.650") < prefix.index("User error: bad syntax")
assert not re.search(r"(?:Runtime error|User error|Assertion failed|Error:)", mixed)
assert mixed.count("M2_30_MIXED_RESUME_SETUP retained_single_survivors=20 single_types_complete=5") == 1
first = prefix + "\n" + mixed

expected = {
    (2,): 6, (3,): 48, (4,): 186, (5,): 684, (6,): 5006,
    (2, 2): 144, (3, 2): 1152, (4, 2): 4464, (3, 3): 9216,
    (2, 2, 2): 1800,
}
found = {}
pattern = r"M2_30_QUOTIENT_TYPE \[ ([\d, ]+) \]\s+(\d+)\s+(\d+)\s+[\d.]+"
for kind, count, passed in re.findall(pattern, first):
    kind = tuple(int(x.strip()) for x in kind.split(","))
    assert kind not in found
    found[kind] = (int(count), int(passed))
assert {k: v[0] for k, v in found.items()} == expected
passed_first = sum(v[1] for v in found.values())
marker = re.findall(r"M2_30_QUOTIENT_RESUME_DONE survivors=(\d+) class_maps=(\d+) single_types=5 mixed_types=5 complete=true", mixed)
assert len(marker) == 1 and int(marker[0][0]) == passed_first
assert len(re.findall(r"M2_30_QUOTIENT_SURVIVOR \[", first)) == passed_first
assert {k for k in found if len(k) == 1} == {(2,), (3,), (4,), (5,), (6,)}
if "--first-only" in sys.argv:
    print(f"M2_30_FIRST_QUOTIENT_COVERAGE complete=true patterns={len(expected)} tests={sum(expected.values())} survivors={passed_first} retained_single_survivors=20 resumed_mixed_survivors={passed_first - 20}")
    sys.exit(0)

other = (base / "logs/n30/other_quotients.log").read_text(errors="replace")
assert not re.search(r"(?:Runtime error|User error|Assertion failed|Error:)", other)

tests = re.findall(r"M2_30_OTHER_QUOTIENT_TEST (\d+) \[ [\d, ]+ \]\s+([23])\s+(true|false)\s+[\d.]+", other)
by_index = {}
for index, variant, result in tests:
    index, variant = int(index), int(variant)
    row = by_index.setdefault(index, {})
    assert variant not in row
    row[variant] = result == "true"
assert set(by_index) == set(range(1, passed_first + 1))
for index, row in by_index.items():
    assert 2 in row
    assert (3 in row) == row[2]
passed2 = sum(row[2] for row in by_index.values())
passed3_indices = {i for i, row in by_index.items() if row.get(3, False)}
survivor_indices = [int(i) for i in re.findall(r"M2_30_THREE_QUOTIENT_SURVIVOR (\d+) ", other)]
assert len(survivor_indices) == len(set(survivor_indices))
assert set(survivor_indices) == passed3_indices
marker = re.findall(r"M2_30_OTHER_QUOTIENTS_DONE input=(\d+) pass2=(\d+) pass3=(\d+) remaining=(\d+) class_maps=(\d+) complete=true", other)
assert len(marker) == 1
input_count, p2, p3, remaining, maps = map(int, marker[0])
assert (input_count, p2, p3, remaining) == (passed_first, passed2, len(passed3_indices), len(passed3_indices))
print(f"M2_30_QUOTIENT_COVERAGE complete=true patterns={len(expected)} tests={sum(expected.values())} first_survivors={passed_first} pass2={passed2} remaining={remaining}")
