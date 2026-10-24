#!/usr/bin/env python3
"""Audit and combine C42 prefix, cubic pairs, and four fixed-input shards.

No partially executed degree-six prefix is used. Every sextic representative
must instead have exactly one completed, indexed nonprincipality witness.
"""
from pathlib import Path
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent
ERROR=re.compile(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',re.I|re.M)
MASTER_SHA='fa129d97164e224e5366e335594c412bd5835e7b7ebc40cc9828478caf30c6cf'

def sha(name):return hashlib.sha256((HERE/name).read_bytes()).hexdigest()

def read(name):
    s=(HERE/name).read_text()
    assert not ERROR.search(s),(name,'Magma error')
    return s

def partitions(n,k=None):
    if n==0:
        yield ();return
    for first in range(min(n,k or n),0,-1):
        for rest in partitions(n-first,first):yield (first,)+rest

def audit():
    prefix=read('logs/sieve/X1_42_p5_prefix_session.log')
    assert sha('data/c42_fixed_representatives.m')==MASTER_SHA
    assert 'NONCUSP counts deg 1..6: [ 0, 0, 4, 183, 648, 2286 ]' in prefix
    assert 'REPS up to diamond: [ 0, 0, 2, 33, 108, 392 ]' in prefix
    for degree,total in ((3,2),(4,33),(5,108)):
        assert len(re.findall(rf'RESULT type=\[ {degree} \] total={total} survivors=0 time=',prefix))==1
    assert [int(i) for i in re.findall(r'C42_EXPORT_INPUT (\d+)',prefix)]==list(range(1,393))
    assert len(re.findall(r'C42_EXPORT_STEP 4 ',prefix))==392
    assert prefix.count('C42_FIXED_REPRESENTATIVES_EXPORTED count=392 roundtrip=true')==1
    assert re.findall(r'C42_CUBIC_PAIR_TEST (\d+) (true|false) ',prefix)==[(str(i),'true') for i in range(1,9)]
    assert prefix.count('C42_CUBIC_PAIR_DONE total=8 survivors=0 complete=true')==1
    assert 'SURVIVOR type=' not in prefix and 'C42_CUBIC_PAIR_SURVIVOR' not in prefix
    witnesses=[];shards=[]
    for shard in range(4):
        name=f'logs/sieve/X1_42_p5_shard{shard}.log';s=read(name)
        assert s.count(f'C42_FIXED_SHARD_SETUP shard={shard}/4 count=98 total=392 genus=25 p=5 qlist=[11,13,17]')==1
        expected=[i for i in range(1,393) if i%4==shard]
        actual=[]
        for match in re.finditer(r'C42_FIXED_CASE shard=(\d+) index=(\d+) nonprincipal_q=(\d+)',s):
            sh,index,q=map(int,match.groups())
            assert sh==shard and q in (11,13,17)
            actual.append(index);witnesses.append({'index':index,'q':q,'shard':shard})
        assert actual==expected,(shard,'missing/duplicate/out-of-order index',len(actual))
        assert 'C42_FIXED_SURVIVOR' not in s
        assert len(re.findall(rf'C42_FIXED_SHARD_DONE shard={shard}/4 tested=98 survivors=0 indices=\[\s*\] complete=true',s))==1
        shards.append({'shard':shard,'indices':actual,'log_sha256':sha(name)})
    assert sorted(w['index'] for w in witnesses)==list(range(1,393))
    return {'complete':True,'expected_partition_types':29,'verified_partition_types':29,
            'nonempty_partition_counts':{'3':2,'4':33,'5':108,'6':392,'3,3':8},
            'prefix_sha256':sha('logs/sieve/X1_42_p5_prefix_session.log'),'master_sha256':MASTER_SHA,
            'shard_source_sha256':sha('sieve_42_fixed_shard.m'),
            'shards':shards,'sextic_witnesses':sorted(witnesses,key=lambda w:w['index'])}

def main():
    coverage=audit()
    (HERE/'data/sieve_42_coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
    counts={(3,):2,(4,):33,(5,):108,(6,):392,(3,3):8}
    lines=['Audited C42 combination: lower-degree prefix and cubic pairs; four disjoint fixed-input sextic shards.',
           'SETUP N=42 M=1 p=5 genus=25 D=6']
    for d in range(1,7):
        for part in partitions(d):
            label='[ '+', '.join(map(str,part))+' ]'
            if part in counts:lines.append(f'RESULT type={label} total={counts[part]} survivors=0')
            else:lines.append(f'TYPE {label}: no divisors')
    lines+=['SUMMARY N=42 M=1 p=5 D=6 shard=0/1 total_survivors=0',
            'COMBINED_C42_DONE complete=true indexed_sextic_witnesses=392', 'SIEVE_DONE']
    (HERE/'logs/sieve/X1_42_p5.log').write_text('\n'.join(lines)+'\n')
    print('C42_COMBINED_AUDIT_DONE complete=true partitions=29 sextic_witnesses=392 survivors=0')

if __name__=='__main__':main()
