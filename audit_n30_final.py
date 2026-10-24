#!/usr/bin/env python3
"""Audit all stages and immutable arithmetic inputs of the last exclusion."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ERROR=re.compile(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',re.I|re.M)

PATHS={'m2_30_quotient_prefix.log':'logs/n30/first_quotient_prefix.log',
       'm2_30_mixed_continuation.log':'logs/n30/first_quotient_mixed_continuation.log',
       'm2_30_other_quotients.log':'logs/n30/other_quotients.log',
       'm2_30_quotient_survivors.m':'data/m2_30_quotient_survivors.m',
       'm2_30_three_quotient_survivors.m':'data/m2_30_three_quotient_survivors.m',
       'filter_a3_m2_30.log':'logs/n30/a3_filter.log',
       'A3_m2_30_arithmetic.txt':'data/A3_m2_30_arithmetic.txt',
       'test_a3_m2_30_map_v4.log':'logs/n30/a3_map_test.log',
       'm2_30_normalization_guard.log':'logs/n30/normalization_guard.log',
       'A3_n30_arithmetic_replay_root.log':'logs/n30/a3_arithmetic_replay.log',
       'replay_a3_projection.log':'logs/n30/a3_projection_replay.log',
       'A3_n30_arithmetic_replay_root.json':'data/A3_n30_arithmetic_replay_root.json',
       'n30_final_coverage.json':'data/n30_final_coverage.json'}
def sha(name):return hashlib.sha256((HERE/PATHS.get(name,name)).read_bytes()).hexdigest()

def read(name,marker):
    s=(HERE/PATHS.get(name,name)).read_text()
    assert not ERROR.search(s),(name,'error')
    assert marker in s,(name,'incomplete')
    return s

def audit():
    run=subprocess.run([sys.executable,str(HERE/'audit_m2_30_quotients.py')],capture_output=True,text=True,check=True)
    assert run.stdout.strip()=='M2_30_QUOTIENT_COVERAGE complete=true patterns=10 tests=22706 first_survivors=68 pass2=10 remaining=10'
    read('filter_a3_m2_30.log','A3_M2_30_FILTER_DONE input=10 remaining=0 class_maps=8 complete=true')
    read('test_a3_m2_30_map_v4.log','A3_M2_30_MAP_TEST_DONE exact=true')
    read('m2_30_normalization_guard.log','M2_30_NORMALIZATION_GUARD_DONE global_a1_zero=true global_a3_zero=true')
    read('A3_n30_arithmetic_replay_root.log','A3_N30_INDEPENDENT_ARITHMETIC_DONE')
    read('replay_a3_projection.log','A3_ELEMENTARY_PROJECTION_REPLAY_DONE input=10 nonzero=10 complete=true')
    rr=json.loads((HERE/PATHS['A3_n30_arithmetic_replay_root.json']).read_text())
    assert rr['source_sha256']==sha('A3_m2_30_arithmetic.txt')
    assert rr['input_divisors_sha256']==sha('m2_30_three_quotient_survivors.m')
    assert rr['independent_integer_lattice_replay'] and rr['divisors_tested']==10 and rr['remaining']==0
    assert (rr['finite_group_order'],rr['rational_cusp_group_order'],rr['quotient_order'])==(69249024,256,270504)
    assert rr['projection_moduli']==[51,1326]
    assert len(rr['nonzero_projection_images'])==10 and all(v!=[0,0] for v in rr['nonzero_projection_images'])
    raw=['m2_30_quotient_prefix.log','m2_30_mixed_continuation.log',
         'm2_30_other_quotients.log','m2_30_quotient_survivors.m',
         'm2_30_three_quotient_survivors.m','filter_a3_m2_30.log','A3_m2_30_arithmetic.txt']
    return {'complete':True,'support_patterns':10,'first_tests':22706,'first_survivors':68,
            'three_map_survivors':10,'remaining':0,'sha256':{name:sha(name) for name in raw}}

def main():
    data=audit();(HERE/PATHS['n30_final_coverage.json']).write_text(json.dumps(data,indent=2)+'\n')
    print('N30_FINAL_AUDIT_DONE complete=true support_patterns=10 tests=22706 stages=68,10,0 remaining=0')

if __name__=='__main__':main()
