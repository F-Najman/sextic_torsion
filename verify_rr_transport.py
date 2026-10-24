#!/usr/bin/env python3
"""Independently replay F2 differential kernels and diamond transports.

The function field, differential-space kernels and geometric cusp action are
exact inputs computed in Magma. This verifies their exported linear algebra,
all transport edges, and the consistency of the artificial local restrictions.
The separate FI replay checks every admissible divisor, without orbit reduction.
"""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

path=Path(sys.argv[1]);source=path.read_text()
assert source.count('print "FI_JET_DATA_READY";') == 1

def value(name):
    matches=re.findall(rf'\b{name}:=(.*?);',source,re.S)
    assert len(matches)==1,(name,len(matches))
    return ast.literal_eval(matches[0])

def parse_matrix(groups):
    nr,nc=int(groups[0]),int(groups[1]);es=ast.literal_eval(groups[2])
    assert len(es)==nr*nc and set(es)<={0,1}
    return [sum(es[nc*i+j]<<j for j in range(nc)) for i in range(nr)],nc

def matrix(name):
    ms=re.findall(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',source,re.S)
    assert len(ms)==1,(name,len(ms))
    return parse_matrix(ms[0])

def echelon(rows):
    out={}
    for x in rows:
        while x:
            b=x.bit_length()-1
            if b not in out:
                out[b]=x;break
            x^=out[b]
    return out

def contains(basis,x):
    while x:
        b=x.bit_length()-1
        if b not in basis:return False
        x^=basis[b]
    return True

def same_span(a,b):
    aa,bb=echelon(a),echelon(b)
    return len(aa)==len(bb) and all(contains(aa,x) for x in b)

def apply(v,M):
    answer=0
    while v:
        low=v&-v;answer^=M[low.bit_length()-1];v^=low
    return answer

def multiply(A,B):return [apply(a,B) for a in A]

def inverse(A):
    n=len(A);left=A[:];right=[1<<i for i in range(n)]
    for i in range(n):
        choices=[j for j in range(i,n) if (left[j]>>i)&1]
        assert choices,'Singular transport matrix'
        j=choices[0];left[i],left[j]=left[j],left[i];right[i],right[j]=right[j],right[i]
        for j in range(n):
            if j!=i and ((left[j]>>i)&1):left[j]^=left[i];right[j]^=right[i]
    assert left==[1<<i for i in range(n)]
    return right

def columns(A,width):
    return [sum(((row>>j)&1)<<i for i,row in enumerate(A)) for j in range(width)]

g=value('gx');D=value('D');ds=value('ds');N=value('N')
if 'TransportOrigins:=' not in source:
    count=0
    for i,d in enumerate(ds,1):
        A,nc=matrix(rf'\bbasisjets\[{i}\]')
        assert len(A)==g and nc==(D//d)*d
        for m in range(1,D//d+1):
            W,width=matrix(rf'\bkernelspaces\[<{i},{m}>\]');assert width==g
            mask=(1<<(m*d))-1;prefix=[row&mask for row in A]
            assert all(apply(w,prefix)==0 for w in W)
            assert len(echelon(W))==g-len(echelon(prefix))
            count+=1
    result={'level':N,'genus':g,'cusps':len(ds),'source_kernel_spaces':count,
            'transport_edges':0,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'independent_kernel_transport_replay':True}
    print(json.dumps(result,sort_keys=True));print('RR_TRANSPORT_INDEPENDENT_DONE')
    sys.exit(0)
origins=value('TransportOrigins');ps=value('TransportPermutations')
assert len(origins)==len(ds)
MT,nc=matrix(r'\bMT');MU,nc2=matrix(r'\bMU')
assert len(MT)==len(MU)==nc==nc2==g
I=[1<<i for i in range(g)]
assert multiply(MT,MU)==multiply(MU,MT)
assert multiply(MT,MT)==multiply(MU,MU)
assert multiply(multiply(MT,MT),multiply(MT,MT))==I
matgroup={tuple(I)};queue=[I]
while queue:
    A=queue.pop()
    for B in (MT,MU):
        AB=multiply(A,B)
        if tuple(AB) not in matgroup:matgroup.add(tuple(AB));queue.append(AB)
assert len(matgroup)==8
CM,nc=matrix(r'\bCartierMatrix');assert len(CM)==nc==g
for M in (MT,MU):assert multiply(M,CM)==multiply(CM,M)
pergroup={tuple(range(len(ds)))};queue=list(pergroup)
for p in ps:
    assert sorted(p)==list(range(1,len(ds)+1))
while queue:
    p=queue.pop()
    for q in ps:
        pq=tuple(q[i]-1 for i in p)
        if pq not in pergroup:pergroup.add(pq);queue.append(pq)
assert len(pergroup)==8
B={};Ts={}
for i in range(1,len(ds)+1):
    B[i],nc=matrix(rf'\bbasisjets\[{i}\]');assert len(B[i])==g and nc==(D//ds[i-1])*ds[i-1]
    Ts[i],nc=matrix(rf'\bTransportMatrices\[{i}\]');assert len(Ts[i])==nc==g
    inverse(Ts[i]);rep=origins[i-1]
    assert 1<=rep<=len(ds) and origins[rep-1]==rep and ds[rep-1]==ds[i-1]
for i in B:
    assert B[i]==multiply(Ts[i],B[origins[i-1]])
source_count=0
for rep in sorted(set(origins)):
    assert Ts[rep]==I
    for m in range(1,D//ds[rep-1]+1):
        W,nc=matrix(rf'\bkernelspaces\[<{rep},{m}>\]');assert nc==g
        mask=(1<<(m*ds[rep-1]))-1
        A=[row&mask for row in B[rep]]
        assert all(apply(w,A)==0 for w in W)
        assert len(echelon(W))==g-len(echelon(A))
        source_count+=1
transport_edges=0
for k,M in enumerate((MT,MU)):
    inv=inverse(M)
    for i in B:
        j=ps[k][i-1];assert ds[i-1]==ds[j-1]
        transported=multiply(inv,B[i])
        for m in range(1,D//ds[i-1]+1):
            width=m*ds[i-1]
            assert same_span(columns(transported,width),columns(B[j],width)),(i,j,m)
        transport_edges+=1
result={'level':N,'genus':g,'cusps':len(ds),'source_kernel_spaces':source_count,
        'transport_edges':transport_edges,'diamond_group_order':len(matgroup),
        'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'independent_kernel_transport_replay':True}
print(json.dumps(result,sort_keys=True))
print('RR_TRANSPORT_INDEPENDENT_DONE')
