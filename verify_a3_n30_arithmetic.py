#!/usr/bin/env python3
"""Independent integer-lattice replay of the last ten C2 x C30 cases."""
import ast
import hashlib
import json
from pathlib import Path
import re
from math import prod
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

HERE=Path(__file__).resolve().parent
source=HERE/'data/A3_m2_30_arithmetic.txt'
s=source.read_text().replace('[*','[').replace('*]',']')
s=re.sub(r'\bfalse\b','False',s);s=re.sub(r'\btrue\b','True',s)
invariants,relations,hgenerators,cuspdegrees,cuspclasses,records=ast.literal_eval(s)
assert invariants==[2,2,2,2,2,2,102,10608,0]
n=len(invariants)-1
for row in relations+hgenerators+cuspclasses:
    assert len(row)==n+1 and row[-1]==0
R=Matrix([row[:-1] for row in relations]).T
assert R.shape==(n,n) and abs(int(R.det()))==prod(invariants[:-1])==69249024
snf=smith_normal_form(R,domain=ZZ)
assert [abs(int(snf[i,i])) for i in range(n)]==invariants[:-1]
B=hermite_normal_form(Matrix([row[:-1] for row in relations+hgenerators]).T)
assert B.shape==(n,n)
quotient_order=abs(int(B.det()));horder=abs(int(R.det()))//quotient_order
assert (horder,quotient_order)==(256,270504)
bsnf=smith_normal_form(B,domain=ZZ)
assert [abs(int(bsnf[i,i])) for i in range(n) if abs(int(bsnf[i,i]))!=1]==[102,2652]
Binv=B.inv();Rinv=R.inv()

def member(inverse,vector):
    assert len(vector)==n+1 and vector[-1]==0
    return all(c.q==1 for c in inverse*Matrix(vector[:-1]))

assert all(member(Binv,v) for v in hgenerators)
assert sorted(cuspdegrees)==[1]*8+[4]*2
assert len(cuspclasses)==len(cuspdegrees)
assert len(records)==10
passed=[]
for index,row in enumerate(records,1):
    i,lam,pushed,claimed,parts=row
    assert i==index and sum(lam)==6 and claimed is False
    assert len(pushed)==n+1 and pushed[-1]==0
    assert sorted(lam)==([6] if index<=4 else [2,2,2])
    summed=[0]*(n+1);expanded=[]
    for sd,mult,td,resmult,weighted in parts:
        assert sd==td*resmult and sd>=1 and td>=1 and mult>=1
        assert len(weighted)==n+1 and weighted[-1]==0
        expanded.extend([sd]*mult)
        for k,v in enumerate(weighted):summed[k]+=mult*v
    assert sorted(expanded)==sorted(lam)
    assert member(Rinv,[a-b for a,b in zip(summed,pushed)])
    passes=member(Binv,pushed)
    assert passes==claimed,(index,'membership mismatch')
    if passes:passed.append(index)
assert not passed
projection=lambda row:[row[6]%51,row[7]%1326]
assert all(projection(row)==[0,0] for row in relations+hgenerators)
images=[projection(row[2]) for row in records]
assert all(v!=[0,0] for v in images)
result={'independent_integer_lattice_replay':True,'finite_group_order':69249024,
        'rational_cusp_group_order':horder,'quotient_order':quotient_order,
        'divisors_tested':10,'remaining':0,
        'projection_moduli':[51,1326],'nonzero_projection_images':images,
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'input_divisors_sha256':hashlib.sha256((HERE/'data/m2_30_three_quotient_survivors.m').read_bytes()).hexdigest()}
(HERE/'data/A3_n30_arithmetic_replay_root.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
print('A3_N30_INDEPENDENT_ARITHMETIC_DONE')
