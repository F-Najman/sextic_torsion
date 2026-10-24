#!/usr/bin/env python3
"""Audit C49 coverage: completed A5/A2 tests and count-derived empty types.

No completion is inferred from the interrupted A5 degree-six header.
The combined record distinguishes executed tests from empty-type deductions.
"""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
ERROR=re.compile(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',re.I|re.M)

def partitions(n,largest=None):
    if n==0:
        yield ()
        return
    for first in range(min(n,largest or n),0,-1):
        for rest in partitions(n-first,first):yield (first,)+rest

def audit():
    prefix=(HERE/'logs/sieve/X1_49_p3_prefix.log').read_text()
    a2=(HERE/'logs/sieve/X1_49_p3_A2.log').read_text()
    assert not ERROR.search(prefix) and not ERROR.search(a2)
    assert 'NONCUSP counts deg 1..6: [ 0, 0, 0, 21, 63, 154 ]' in prefix
    assert 'REPS up to diamond: [ 0, 0, 0, 1, 3, 8 ]' in prefix
    done=set()
    for match in re.finditer(r'TYPE \[([^\]]+)\]: no divisors|RESULT type=\[([^\]]+)\] total=\d+ survivors=0',prefix):
        done.add(tuple(map(int,(match.group(1) or match.group(2)).replace(',',' ').split())))
    low={p for d in range(1,6) for p in partitions(d)}
    assert done==low and len(low)==18
    assert 'RESULT type=[ 4 ] total=1 survivors=0' in prefix
    assert 'RESULT type=[ 5 ] total=3 survivors=0' in prefix
    tests=re.findall(r'C49_A2_TEST i=(\d+) principal=(\w+)',a2)
    assert tests==[(str(i),'false') for i in range(1,9)]
    assert 'C49_A2_PLACES total=154 representatives=8' in a2
    assert 'C49_A2_DONE total=8 survivors=0 data=[]' in a2
    counts=[0,0,0,21,63,154]
    empty={p for p in partitions(6) if any(counts[d-1]==0 for d in p)}
    assert len(empty)==10 and set(partitions(6))-empty=={(6,)}
    expected={p for d in range(1,7) for p in partitions(d)}
    verified=low|empty|{(6,)}
    assert verified==expected and len(expected)==29
    return {'complete':True,'expected_partition_types':29,
            'verified_partition_types':29,'missing_types':[],
            'executed_A5_types':sorted(low),'executed_A2_types':[[6]],
            'derived_empty_types':sorted(empty),'noncuspidal_place_counts':counts,
            'A2_representatives_tested':8,'annihilator':'A2',
            'prefix_sha256':hashlib.sha256(prefix.encode()).hexdigest(),
            'a2_sha256':hashlib.sha256(a2.encode()).hexdigest()}

def main():
    record=audit()
    prefix=(HERE/'logs/sieve/X1_49_p3_prefix.log').read_text()
    a2=(HERE/'logs/sieve/X1_49_p3_A2.log').read_text()
    combined=('COMBINED C49 RECORD: 18 completed A5 types, one completed A2 type, ten count-derived empty types.\n'
              'The interrupted A5 [6] header is not a completed test.\n\n'
              '=== ORIGINAL PREFIX ===\n'+prefix+'\n=== COMPLETED A2 TEST ===\n'+a2+
              '\n=== DERIVED EMPTY DEGREE-SIX TYPES ===\n'+json.dumps(record['derived_empty_types'])+'\n'
              'Each listed type contains a part1,2,or3, where no noncuspidal places exist.\n'
              'COMBINED_C49_A2_DONE full_partition_coverage=true survivors=0\n')
    (HERE/'logs/sieve/X1_49_p3.log').write_text(combined)
    (HERE/'data/sieve_49_coverage.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
