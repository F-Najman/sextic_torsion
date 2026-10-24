"""Independent integer-only A3 obstruction, without a subgroup solver.

Project coordinates 7 and 8 of the recorded Picard group to C51 and
C1326. Every relation and every rational cusp generator maps to zero;
all ten candidate divisor classes have nonzero image.
"""
from pathlib import Path
import ast

base = Path(__file__).resolve().parent
source = (base / "data/A3_m2_30_arithmetic.txt").read_text()
source = source.replace("[*", "[").replace("*]", "]")
source = source.replace("false", "False").replace("true", "True")
invariants, relations, cusp_generators, cusp_degrees, cusp_classes, records = ast.literal_eval(source)
assert invariants == [2, 2, 2, 2, 2, 2, 102, 10608, 0]
assert len(records) == 10

def project(vector):
    assert len(vector) == 9 and vector[-1] == 0
    return vector[6] % 51, vector[7] % 1326

assert all(project(v) == (0, 0) for v in relations + cusp_generators)
assert [row[0] for row in records] == list(range(1, 11))
for index, kind, group_class, passed, parts in records:
    assert sum(kind) == 6 and not passed
    assert project(group_class) != (0, 0)
    assert sum(degree * multiplicity for degree, multiplicity, _, _, _ in parts) == 6
    assert all(degree == target_degree * residue_multiplicity for degree, _, target_degree, residue_multiplicity, _ in parts)
    print(f"A3_PROJECTION_WITNESS index={index} type={kind} image={project(group_class)}")
print("A3_ELEMENTARY_PROJECTION_REPLAY_DONE input=10 nonzero=10 complete=true")
