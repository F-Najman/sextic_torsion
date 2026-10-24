#!/usr/bin/env python3
"""Replay exact differential-space annihilators and optional diamond transport.

The function-field DifferentialSpace outputs remain geometric Magma input.
This independently checks every kernel equality, prefix, transport identity,
and generated matrix/permutation group over F2. It does not recompute RR.
"""
import ast
from pathlib import Path
import re
import sys

s=Path(sys.argv[1]).read_text()
assert s.count('print "FI_JET_DATA_READY";') == 1

def value(name):
    m=re.findall(rf'\b{name}:=(.*?);',s,re.S)
    assert len(m)==1,name
    return ast.literal_eval(m[0])

def matrix(pattern):
    ms=re.findall(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S)
    assert len(ms)==1,pattern
    nr,nc=map(int,ms[0][:2]);es=ast.literal_eval(ms[0][2])
    assert len(es)==nr*nc and set(es)<={0,1}
    return [sum(es[nc*i+j]<<j for j in range(nc)) for i in range(nr)],nc

def pivots(rows):
    p={}
    for x in rows:
        while x:
            k=x.bit_length()-1
            if k in p:x^=p[k]
            else:p[k]=x;break
    return p

def apply(x,rows):
    y=0
    while x:
        b=x&-x;y^=rows[b.bit_length()-1];x^=b
    return y

def mul(A,B):return [apply(x,B) for x in A]

N,D,g,ds=(value(k) for k in ('N','D','gx','ds'))
Bs={}
for i,d in enumerate(ds,1):
    A,nc=matrix(rf'\bbasisjets\[{i}\]')
    assert len(A)==g and nc==(D//d)*d
    Bs[i]=A
keys=[tuple(map(int,m)) for m in re.findall(r'\bkernelspaces\[<(\d+),(\d+)>\]:=',s)]
assert len(keys)==len(set(keys))
for i,m in keys:
    W,nc=matrix(rf'\bkernelspaces\[<{i},{m}>\]')
    assert nc==g and len(pivots(W))==len(W)
    A=[row&((1<<(m*ds[i-1]))-1) for row in Bs[i]]
    assert all(apply(w,A)==0 for w in W),(i,m,'orthogonality')
    assert len(pivots(A))+len(W)==g,(i,m,'kernel dimension')
transported=0
if 'TransportOrigins:=' in s:
    origins=value('TransportOrigins');ps=value('TransportPermutations')
    assert len(origins)==len(ds) and len(ps)==2
    steps=value('CoefficientSteps')
    for p in ps:
        assert sorted(p)==list(range(1,len(ds)+1))
        assert all(ds[i]==ds[p[i]-1] and steps[i]==steps[p[i]-1] for i in range(len(ds)))
    MT,nt=matrix(r'\bMT');MU,nu=matrix(r'\bMU');assert nt==nu==g
    I=[1<<i for i in range(g)]
    assert mul(MT,MU)==mul(MU,MT)
    assert mul(mul(MT,MT),mul(MT,MT))==I
    assert mul(mul(MU,MU),mul(MU,MU))==I
    invs=[mul(mul(M, M),M) for M in (MT,MU)]
    group={tuple(I):tuple(range(1,len(ds)+1))};queue=[tuple(I)]
    while queue:
        cur=queue.pop();cp=group[cur]
        for L,p in zip(invs,ps):
            new=tuple(mul(L,cur));np=tuple(p[j-1] for j in cp)
            if new in group:assert group[new]==np
            else:group[new]=np;queue.append(new)
    assert len(group)==8
    for i,origin in enumerate(origins,1):
        L,nc=matrix(rf'\bTransportMatrices\[{i}\]')
        assert len(L)==nc==g and tuple(L) in group
        assert group[tuple(L)][origin-1]==i
        assert mul(L,Bs[origin])==Bs[i],(i,'transport identity')
        assert all((origin,m) in keys for m in range(1,D//ds[i-1]+1))
        transported+=1
else:
    assert set(keys)=={(i,m) for i,d in enumerate(ds,1) for m in range(1,D//d+1)}
print({'level':N,'genus':g,'exact_kernel_equalities':len(keys),
       'transported_blocks':transported,'independent_kernel_transport_replay':True})
