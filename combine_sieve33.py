#!/usr/bin/env python3
"""Validate exhaustive partition coverage across the interrupted and resumed run.

The prefix proves all completed types before interruption. A filtered resume
checks the remaining types with cached Hecke images. A completed type requires
either 'no divisors' or an explicit zero-survivor RESULT. A TYPE header alone
does not prove anything. No original output or completion marker is invented.
"""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n,largest or n),0,-1):
        for rest in partitions(n-first,first): yield (first,)+rest
def completed(text):
    out=set()
    for mat in re.finditer(r'TYPE \[([^\]]+)\]: no divisors|RESULT type=\[([^\]]+)\] total=\d+ survivors=0',text):
        out.add(tuple(map(int,(mat.group(1) or mat.group(2)).replace(',',' ').split())))
    return out
prefix=(HERE/'logs/sieve/X1_33_p5_prefix.log').read_text()
resume=(HERE/'logs/sieve/X1_33_p5_resume.log').read_text()
for text in (prefix,resume):
    assert not re.search(r'Runtime error|User error|Assertion failed|syntax error|^ERROR',text,re.M|re.I)
expected={p for d in range(1,7) for p in partitions(d)}
done=completed(prefix)|completed(resume)
complete=expected<=done and 'SIEVE_RESUME_DONE FromDegree=6 MaxPart=3' in resume
record={'complete':complete,'expected_partition_types':len(expected),
        'verified_partition_types':len(expected&done),
        'missing_types':sorted(expected-done),
        'prefix_sha256':hashlib.sha256(prefix.encode()).hexdigest(),
        'resume_sha256':hashlib.sha256(resume.encode()).hexdigest()}
combined=('COMBINED COMPUTATIONAL RECORD: original prefix plus cached continuation.\n'
          'The continuation is filtered; completeness is established by partition coverage.\n\n'
          '=== ORIGINAL PREFIX ===\n'+prefix+'\n=== CACHED CONTINUATION ===\n'+resume)
if complete:combined+='\nCOMBINED_SIEVE_DONE coveragevalidated=true all_survivors=0\n'
(HERE/'logs/sieve/X1_33_p5.log').write_text(combined)
(HERE/'data/sieve_33_coverage.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
