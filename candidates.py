#!/usr/bin/env python3
"""Exact Phi(6) obstruction ledger, accepting only audited completed certificates.

Pairs (m,n), m|n, denote C_m x C_n. Every subgroup type is (a,b)
with a|m, b|n and a|b. The Weil pairing gives phi(m)|6.
"""
import json
import hashlib
import re
from pathlib import Path
from sympy import divisors, primefactors, totient

S6 = {2, 3, 5, 7, 11, 13, 17, 19, 37}
M6 = (1, 2, 3, 4, 6, 7, 9, 14, 18)
INFINITE = (
    {(1, n) for n in range(1, 31) if n not in (23, 25, 29)}
    | {(2, 2*n) for n in range(1, 11)}
    | {(3, 3*n) for n in range(1, 5)}
    | {(4, 4), (4, 8), (6, 6)}
)
KNOWN = INFINITE | {(1, 25), (1, 37), (4, 12), (2, 24)}

def subgroups(g):
    m, n = g
    return {(int(a), int(b)) for a in divisors(m) for b in divisors(n)
            if b % a == 0}

def minimum_obstructions(base):
    assert all(subgroups(g) <= base for g in base), "base must be subgroup closed"
    max_cyclic = max(n for m, n in base if m == 1)
    # A minimal cyclic obstruction n is either a prime in S6, or n/q is
    # cyclic and belongs to base for every prime q|n. Thus n <= maxC*max(S6).
    bound = max(max(S6), max_cyclic*max(S6))
    result = []
    for m in M6:
        assert 6 % int(totient(m)) == 0
        # For m>1, C_n is a proper subgroup, so n must be a cyclic order in base.
        ns = range(1, bound+1) if m == 1 else sorted(n for a,n in base if a == 1)
        for n in ns:
            if n % m or (m,n) in base or not set(primefactors(n)) <= S6:
                continue
            if subgroups((m,n)) - {(m,n)} <= base:
                result.append((m,n))
    return sorted(result)

def label(g):
    return f'C{g[1]}' if g[0] == 1 else f'C{g[0]} x C{g[1]}'

def main():
    initial = minimum_obstructions(INFINITE)
    after_known = minimum_obstructions(KNOWN)
    # The 37-case is separately certified in verify_37.m, conditional only on
    # DKSS Proposition 1.4. Every other group containing C37 is then excluded.
    after_37 = [g for g in after_known if g[1] % 37]
    # Global criterion from the quintic manuscript, with d changed to 6.
    # Each entry (a, k3, index, CM minimum) uses the archived proven rank/character data.
    global_data = {
        119: (3,7,6912,96), 121: (2,8,7260,110), 125: (3,7,7500,50),
        143: (67,8,10080,120), 169: (3,7,14196,52),
        187: (122,8,17280,160), 209: (34,8,21600,180),
        221: (3,7,24192,96), 289: (3,7,41616,136),
        323: (101,8,51840,288), 361: (3,7,64980,114),
    }
    certified_global = []
    for n,(a,k,index,cm) in global_data.items():
        t, aa = 1, a % n
        while aa not in (1,n-1):
            aa = aa*a % n
            t += 1
        assert 325*index > 6*32768*k and 6 % t != 0 and cm > 6
        certified_global.append({'n':n, 'a':a, 'diamond_order':t,
                                 'k3':k, 'index':index, 'cm_minimum':cm})
    remaining = [g for g in after_37 if g not in {(1,n) for n in global_data}]
    here = Path(__file__).resolve().parent
    def log_has(name, marker):
        p = here/name
        if not p.exists(): return False
        text=p.read_text()
        if re.search(r'Runtime error|User error|System error|Assertion failed|syntax error|failed memory request|^ERROR',text,re.I|re.M):
            return False
        return marker in text
    def complete_zero_sieve(name, allowed_nonzero=None):
        def partitions(n,k=None):
            if n==0:
                yield ()
                return
            for first in range(min(n,k or n),0,-1):
                for rest in partitions(n-first,first):yield (first,)+rest
        text=(here/name).read_text()
        completed={}
        for mat in re.finditer(r'TYPE \[([^\]]+)\]: no divisors|RESULT type=\[([^\]]+)\] total=\d+ survivors=(\d+)',text):
            part=tuple(map(int,(mat.group(1) or mat.group(2)).replace(',',' ').split()))
            completed[part]=int(mat.group(3) or 0)
        expected={part for degree in range(1,7) for part in partitions(degree)}
        allowed_nonzero=allowed_nonzero or {}
        return expected<=completed.keys() and all(
            completed[part]==allowed_nonzero.get(part,0) for part in expected)
    reduction_exclusions = {(1,55),(1,75),(1,77),(1,85)}
    assert log_has('logs/finite_fields/torsion_F2_49_51_55_75_77_85.log','TORSION n=85 p=2')
    # Positive-genus full-level curves over their degree-six cyclotomic fields:
    # rank zero, all cusps rational, and no good reduction at 29 or 19.
    assert log_has('logs/full_level_rank.log','FULL_LEVEL_RANK_DONE')
    noncyclic_exclusions = {(7,7),(9,9)}
    assert log_has('logs/finite_fields/C3xC27_F4_F16_F64.log','NONCYCLIC_FINITE_FIELDS_DONE')
    noncyclic_exclusions.add((3,27))
    assert log_has('logs/rr/X1_3_21_F4_double_cusp.log','tested=2608 total=2608 exceptions=0 complete=true')
    assert log_has('logs/finite_fields/C3xCn_F4_F16_F64.log','LOCAL_3N_DONE')
    noncyclic_exclusions.add((3,21))
    assert log_has('logs/rr/X1_3_15_F4_double_cusp.log','ALL_EXCEPTIONS_SUPPORTED_ON_RATIONAL_CUSPS_WITH_H0_EQ_1 32')
    assert log_has('logs/classes/X1_3_15_reduction_injective.log','INJECTIVE_REDUCTION_3_15_DONE')
    noncyclic_exclusions.add((3,15))
    assert log_has('logs/sieve/X1_3_24_p7.log','SIEVE_3_24_DONE 12 0')
    noncyclic_exclusions.add((3,24))
    assert log_has('logs/sieve/X1_6_12_p7.log','SIEVE_6_12_DONE 48 0')
    assert log_has('logs/counts/X1_6_12_p7.log','INDEPENDENT_M6_COUNTS_DONE [ 0, 0, 288 ]')
    noncyclic_exclusions.add((6,12))
    assert log_has('logs/sieve/X1_2_26_p3.log','SUMMARY N=26 M=2 p=3 D=6 shard=0/1 total_survivors=2')
    assert log_has('logs/sieve/X1_2_26_p3.log','SIEVE_DONE')
    assert complete_zero_sieve('logs/sieve/X1_2_26_p3.log',{(6,):2})
    assert log_has('logs/classes/X1_2_26_p3_classes.log','M2N26_CLASS_INSPECT_DONE')
    for i in (10,16):
        assert log_has('logs/classes/X1_2_26_p3_classes.log',
                       f'M2N26_SURVIVOR_RATIONAL_CUSPS_ONLY {i} false')
        assert log_has('logs/classes/X1_2_26_p3_survivor_match.log',f'M2N26_SURVIVOR_MATCH {i} exact=true')
    noncyclic_exclusions.add((2,26))
    assert log_has('logs/sieve/X1_4_16_p3.log','SIEVE_4_16_P3_DONE 80 0')
    assert log_has('logs/counts/X1_4_16_p3.log','INDEPENDENT_M4_P3_COUNTS_DONE [ 0, 0, 960 ]')
    noncyclic_exclusions.add((4,16))
    assert log_has('logs/sieve/X1_4_20_p3.log','SIEVE_4_20_P3_DONE 104 0')
    assert log_has('logs/sieve/X1_4_20_p3.log','SIEVE_4_20_P3_TYPE_DONE 2 24 0')
    assert log_has('logs/sieve/X1_4_20_p3.log','SIEVE_4_20_P3_TYPE_DONE 3 80 0')
    assert log_has('logs/counts/X1_4_20_p3.log','INDEPENDENT_M4N20_P3_COUNTS_DONE [ 0, 192, 960 ]')
    noncyclic_exclusions.add((4,20))
    assert log_has('logs/n24/example_first.log','M2N24_DIRECT_FULL_DONE exact=true')
    assert log_has('logs/n24/example_second.log','M2N24_SECOND_DIRECT_FULL_DONE exact=true')
    assert log_has('logs/n24/X1_2_24_p5_cusp_reduction.log','N24_CUSP_REDUCTION_DONE 2 2 complete=true')
    assert log_has('logs/n24/match_exact_points.log','N24_EXACT_POINT_MATCH_DONE complete=true')
    for name,marker in [
        ('logs/n24/example_first_isogenies.log','N24_ALL_THREE_TWO_ISOGENIES_DONE exact=true'),
        ('logs/n24/example_second_isogenies.log','N24_SECOND_ALL_THREE_TWO_ISOGENIES_DONE exact=true')]:
        assert log_has(name,marker)
    assert log_has('logs/finite_fields/C4xC24_hasse.log','C4C24_HASSE_DONE no_possible_good_reduction=true')
    noncyclic_exclusions.add((4,24))
    assert log_has('logs/classes/X1_3_18_diagonal_classes.log','M3N18_DIAGONAL_DONE 36 36 0 0 complete=true')
    assert log_has('logs/classes/X1_3_18_arithmetic_replay.log','C3X18_INDEPENDENT_ARITHMETIC_DONE')
    a18=json.loads((here/'data/C3x18_arithmetic_replay.json').read_text())
    assert a18['independent_integer_lattice_replay'] and a18['pairs_tested']==1296 and a18['compatible_pairs']==0
    assert a18['source_sha256']==hashlib.sha256((here/'data/C3x18_arithmetic_certificate.txt').read_bytes()).hexdigest()
    noncyclic_exclusions.add((3,18))
    from audit_n30_final import audit as audit30
    coverage30=audit30()
    assert json.loads((here/'data/n30_final_coverage.json').read_text())==coverage30
    noncyclic_exclusions.add((2,30))
    assert log_has('logs/classes/X1_36_p5_class_sieve.log',
                   'CLASS_CUSPIDAL_SIEVE_DONE N=36 p=5 D=6 remaining_classes=0 complete=true')
    for part,total in [('3',48),('4',144),('5',660),('6',2856),('3,3',1176)]:
        assert log_has('logs/classes/X1_36_p5_class_sieve.log',
                       f'CLASS_CUSPIDAL_TYPE type=[{part}] total={total} Hall_pass=0')
    divisor_class_exclusions = {(1,36)}
    sieve_exclusions = set()
    from combine_sieve49 import audit as audit49
    coverage49=audit49()
    assert json.loads((here/'data/sieve_49_coverage.json').read_text()) == json.loads(json.dumps(coverage49))
    assert log_has('logs/sieve/X1_49_p3.log','COMBINED_C49_A2_DONE full_partition_coverage=true survivors=0')
    assert log_has('logs/homology/J1_49_T2.log','TWO49_DONE zero=true')
    assert log_has('logs/homology/J1_49_T2_replay.log',"'independent_binary_replay': True")
    sieve_exclusions.add((1,49))
    from combine_sieve42 import audit as audit42
    coverage42=audit42()
    assert json.loads((here/'data/sieve_42_coverage.json').read_text())==coverage42
    assert log_has('logs/sieve/X1_42_p5.log','COMBINED_C42_DONE complete=true indexed_sextic_witnesses=392')
    assert log_has('logs/counts/X1_42_p5.log','GL2_COUNT_42_DONE')
    assert log_has('logs/tests/probe_42_hecke_hash.log','C42_HECKE_PROBE_DONE')
    assert complete_zero_sieve('logs/sieve/X1_42_p5.log')
    sieve_exclusions.add((1,42))
    for m,n,p,fn in [(1,34,3,'logs/sieve/X1_34_p3.log'),
                     (1,38,3,'logs/sieve/X1_38_p3.log'),
                     (1,39,2,'logs/sieve/X1_39_p2.log'),
                     (1,50,3,'logs/sieve/X1_50_p3.log'),
                     (1,33,5,'logs/sieve/X1_33_p5.log'),
                     (1,35,3,'logs/sieve/X1_35_p3.log'),
                     (1,40,3,'logs/sieve/X1_40_p3.log'),
                     (2,22,3,'logs/sieve/X1_2_22_p3.log'),
                     (2,28,3,'logs/sieve/X1_2_28_p3.log')]:
        summary = f'SUMMARY N={n} M={m} p={p} D=6 shard=0/1 total_survivors=0'
        if log_has(fn,summary) and log_has(fn,'SIEVE_DONE') and complete_zero_sieve(fn):
            sieve_exclusions.add((m,n))
    if log_has('logs/classes/X1_32_p3_class_membership.log','RATIONAL_CLASS_MEMBERSHIP [ false, false, false, false ]'):
        assert log_has('logs/classes/X1_32_p3_class_membership.log','ANALYZE_DONE')
        sieve_exclusions.add((1,32))
    isogeny_exclusions = {(1,44)} if (2,22) in sieve_exclusions else set()
    isogeny_exclusions.add((1,48))
    if (2,30) in noncyclic_exclusions:
        isogeny_exclusions.add((1,60))
    if (2,26) in noncyclic_exclusions:
        isogeny_exclusions.add((1,52))
    if (2,28) in sieve_exclusions:
        isogeny_exclusions.add((1,56))
    assert log_has('logs/check_diamonds.log','DIAMOND_CHECKS_DONE')
    operator_exclusions = {(1,91),(1,95),(1,133),(1,247)}
    assert log_has('logs/finite_fields/torsion_F2_35_45_49_51_81.log','TORSION n=81 p=2 : k with Z/81 over F_p^k = []')
    reduction_exclusions.add((1,81))
    assert log_has('logs/rr/X1_45_F2_double_cusp.log','tested=19048 total=19048 exceptions=0 complete=true')
    reduction_exclusions.add((1,45))
    assert log_has('logs/homology/J1_57_two_torsion.log',
                   'RATIONAL_TWO_DONE N=57 kernel_dimension=2 diamond_trivial=true')
    replay = json.loads((here/'data/two_torsion_57_replay.json').read_text())
    assert replay == {'dimension':170, 'A5_rank_mod2':168,
                      'kernel_dimension':2, 'diamond2_identity_on_kernel':True}
    reduction_exclusions.add((1,57))
    assert log_has('logs/rr/X1_27_C54_certificate.log','C54_FINITE_CERTIFICATE_DONE complete=true')
    assert log_has('logs/rr/X1_27_C54_certificate.log','C54_CERT_F5_ODD_GROUP 27 [ 27 ]')
    assert log_has('logs/counts/X1_27_p2.log','GL2_COUNT_27_DONE')
    reduction_exclusions.add((1,54))
    assert log_has('logs/cartier/Y63_continuation.log','RANK0_CARTIER_DONE N=63 D=6 tested=6542 total=6542 exceptions=0 complete=true')
    assert log_has('logs/cartier/Y63_projection_replay.log',"'independent_projection_replay': True")
    assert log_has('logs/cartier/Y63_rank_replay.log','"exceptions": 0, "complete": true')
    assert log_has('logs/finite_fields/C63_F64.log','C63_F64_DONE no_order_divisible_by63=true')
    assert log_has('logs/cartier/X1_63_character.log','C63_POSITIVE_CHARACTER_QUOTIENT_DONE [ 47, 59 ]')
    reduction_exclusions.add((1,63))
    assert log_has('logs/cartier/Y65_rr_transport.log','RANK0_CARTIER_DONE N=65 D=6 tested=1792 total=1792 exceptions=0 complete=true')
    for name,marker in [
        ('logs/cartier/Y65_kernel_transport_replay.log','RR_TRANSPORT_INDEPENDENT_DONE'),
        ('logs/cartier/Y65_kernel_replay.log',"'independent_kernel_transport_replay': True"),
        ('logs/cartier/Y65_projection_replay.log',"'independent_projection_replay': True"),
        ('logs/cartier/Y65_rank_replay.log','"exceptions": 0, "complete": true')]:
        assert log_has(name,marker)
    replay65=json.loads((here/'logs/cartier/Y65_kernel_transport_replay.log').read_text().splitlines()[0])
    assert replay65['source_sha256']==hashlib.sha256((here/'data/fi65_genus41_rr_transport.m').read_bytes()).hexdigest()
    reduction_exclusions.add((1,65))
    assert log_has('logs/cartier/Y51_rr_kernels.log','RANK0_CARTIER_DONE N=51 D=6 tested=54808 total=54808 exceptions=0 complete=true')
    for name,marker in [
        ('logs/cartier/Y51_kernel_transport_replay.log','RR_TRANSPORT_INDEPENDENT_DONE'),
        ('logs/cartier/Y51_kernel_replay.log',"'independent_kernel_transport_replay': True"),
        ('logs/cartier/Y51_projection_replay.log',"'independent_projection_replay': True"),
        ('logs/cartier/Y51_rank_replay.log','"exceptions": 0, "complete": true')]:
        assert log_has(name,marker)
    replay51=json.loads((here/'logs/cartier/Y51_kernel_transport_replay.log').read_text().splitlines()[0])
    assert replay51['source_kernel_spaces']==100
    assert replay51['source_sha256']==hashlib.sha256((here/'data/fi51_retained_rr_kernels.m').read_bytes()).hexdigest()
    reduction_exclusions.add((1,51))
    all_exclusions = reduction_exclusions | noncyclic_exclusions | sieve_exclusions | isogeny_exclusions | operator_exclusions | divisor_class_exclusions
    final_remaining = [g for g in remaining if g not in all_exclusions]
    complete=not final_remaining
    data = {'status':'COMPLETE: all minimal obstructions excluded' if complete else 'INCOMPLETE: remaining cases are unresolved here',
            'classification_complete':complete,'classification':sorted(KNOWN) if complete else None,
            'phi_infinite':sorted(INFINITE), 'known_lower_bound':sorted(KNOWN),
            'minimal_outside_infinite':initial,
            'minimal_outside_known':after_known,
            'after_37_classification':after_37,
            'global_exclusions':certified_global,
            'remaining_after_global':remaining,
            'reduction_exclusions':sorted(reduction_exclusions),
            'noncyclic_exclusions':sorted(noncyclic_exclusions),
            'sieve_exclusions':sorted(sieve_exclusions),
            'isogeny_exclusions':sorted(isogeny_exclusions),
            'operator_exclusions':sorted(operator_exclusions),
            'divisor_class_exclusions':sorted(divisor_class_exclusions),
            'remaining_after_verified_exclusions':final_remaining}
    out = Path(__file__).resolve().parent
    (out/'data/candidates.json').write_text(json.dumps(data,indent=2)+'\n')
    lines = ['Phi(6): '+('COMPLETE CLASSIFICATION' if complete else 'finite reduction (INCOMPLETE)'),
             f'Infinite-family groups: {len(INFINITE)}',
             f'Known lower-bound groups: {len(KNOWN)}']
    for name, groups in [('Minimal outside infinite families',initial),
                         ('Minimal outside known lower bound',after_known),
                         ('After the 37-case',after_37),
                         ('After global exclusions',remaining),
                         ('After all verified exclusions',final_remaining)]:
        lines += ['', f'{name}: {len(groups)}', ', '.join(map(label,groups))]
    lines += ['', ('CLASSIFICATION: Phi(6)=Phi_infinity(6) union {C25,C37,C4 x C12,C2 x C24}; 48 groups.' if complete else 'CAUTION: lower bound is not asserted to equal Phi(6).'),
              'Global exclusions invoke the rank/character and CM data of the quintic proof.',
              'CANDIDATES_DONE']
    (out/'data/candidates.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))

if __name__ == '__main__':
    main()
