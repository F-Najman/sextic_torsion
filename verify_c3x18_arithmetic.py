#!/usr/bin/env python3
"""Independent integer-lattice replay of the two-prime C3x18 obstruction.

Reconstruct the diagonal cusp subgroup from exact relation matrices;
Hermite normal forms, not Magma group-coordinate assumptions, test all pairs.
This does not recompute the geometry, class maps, or Hecke correspondences.
"""
import ast,hashlib,json,math
from pathlib import Path
from sympy import Matrix,ilcm
from sympy.matrices.normalforms import hermite_normal_form
HERE=Path(__file__).resolve().parent
source=(HERE/'data/C3x18_arithmetic_certificate.txt').read_text()
d=ast.literal_eval(source.replace('[*','[').replace('*]',']'))
assert len(d)==13
primes,zetas,signatures,degrees,base,inv7,rel7,inv13,rel13,cusp7,cusp13,state7,state13=d
assert primes==[7,13] and all((z*z+z+1)%p==0 and z%p!=1 for p,z in zip(primes,zetas))
assert len(signatures)==24 and len(set(map(tuple,signatures)))==24
assert sorted(degrees)==[1]*18+[3]*6 and degrees[base-1]==1
assert len(cusp7)==len(cusp13)==24 and len(state7)==len(state13)==36

def finite(inv,rels,cusps,states):
    assert inv[-1]==0 and all(n>0 for n in inv[:-1])
    n=len(inv)
    assert all(len(r)==n and r[-1]==0 for r in rels+cusps+states)
    return [r[:-1] for r in rels],[r[:-1] for r in cusps],[r[:-1] for r in states]
r7,c7,s7=finite(inv7,rel7,cusp7,state7)
r13,c13,s13=finite(inv13,rel13,cusp13,state13)
n7,n13=len(inv7)-1,len(inv13)-1

def basis(rows,n):
    H=hermite_normal_form(Matrix(rows).T)
    assert H.shape==(n,n) and H.det()!=0
    return H

g7=abs(int(basis(r7,n7).det()));g13=abs(int(basis(r13,n13).det()))
assert g7==math.prod(inv7[:-1]) and g13==math.prod(inv13[:-1])
h7=g7//abs(int(basis(r7+c7,n7).det()))
h13=g13//abs(int(basis(r13+c13,n13).det()))
rows=[r+[0]*n13 for r in r7]+[[0]*n7+r for r in r13]+[a+b for a,b in zip(c7,c13)]
H=basis(rows,n7+n13)
index=abs(int(H.det()));assert (g7*g13)%index==0
hdiag=(g7*g13)//index
assert h7==h13==hdiag==20253807
Hinv=H.inv();den=int(ilcm(*[x.q for x in Hinv]));M=Hinv*den
assert all(x.q==1 for x in M)
def transform(v):return [int(x)%den for x in M*Matrix(v)]
left=[transform([9*x for x in a]+[0]*n13) for a in s7]
right=[transform([0]*n7+[9*x for x in b]) for b in s13]
compatible=[(i+1,j+1) for i,a in enumerate(left) for j,b in enumerate(right)
            if all((x+y)%den==0 for x,y in zip(a,b))]
assert compatible==[]
record={'independent_integer_lattice_replay':True,'states':[36,36],
        'cusp_group_order':hdiag,'multiplier':9,'pairs_tested':1296,
        'compatible_pairs':0,'source_sha256':hashlib.sha256(source.encode()).hexdigest()}
(HERE/'data/C3x18_arithmetic_replay.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
print('C3X18_INDEPENDENT_ARITHMETIC_DONE')
