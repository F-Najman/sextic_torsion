/* Exact double-cuspidal RR test via intersections of differential spaces.
   For every cusp c and admissible multiplicity m, precompute Omega(2mc)
   inside the common space of holomorphic differentials.  Then
       l(2C) = 2D + 1 - g + dim(intersection_c Omega(2m_c c)).
   This avoids repeating function-field RR calculations for every C.
   The first DirectCheck cases are checked independently with ordinary RR;
   every exceptional case is also checked with ordinary RR. */
SetColumns(0);
SetMemoryLimit(2*10^9);
if not assigned N then N:=45;
elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned D then D:=6;
elif Type(D) eq MonStgElt then D:=StringToInteger(D); end if;
if not assigned Limit then Limit:=0;
elif Type(Limit) eq MonStgElt then Limit:=StringToInteger(Limit); end if;
if not assigned DirectCheck then DirectCheck:=5;
elif Type(DirectCheck) eq MonStgElt then DirectCheck:=StringToInteger(DirectCheck); end if;
if not assigned UseJets then UseJets:=1;
elif Type(UseJets) eq MonStgElt then UseJets:=StringToInteger(UseJets); end if;
AttachSpec("mdmagma/v2/mdmagma.spec");
t0:=Cputime();
X:=MDX1(N,GF(2));
C:=Curve(X);
g:=Genus(X);
cusps:=[c:c in Cusps(X)|Degree(c) le D];
degs:=[Degree(c):c in cusps];
printf "CUSP_DOUBLE_LOCAL_REF_SETUP N=%o genus=%o cusp_degrees=%o time=%o\n",
 N,g,degs,Cputime(t0);

function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else []; end if;
 result:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do
   Append(~result,[k] cat v);
  end for;
 end for;
 return result;
end function;
vs:=WeightedVectors(degs,D);
printf "CUSP_DOUBLE_FAST_COUNT N=%o effective_divisors=%o time=%o\n",N,#vs,Cputime(t0);
Omega,om:=SpaceOfHolomorphicDifferentials(C);
assert Dimension(Omega) eq g;
W:=VectorSpace(GF(2),g);
localspaces:=AssociativeArray();
omega0:=om(Basis(Omega)[1]);
fixedratios:=[om(b)/omega0:b in Basis(Omega)];
ratios:=fixedratios;
for i in [1..#cusps] do
 if UseJets eq 1 then
  printf "CUSP_DOUBLE_CUSP_BEGIN N=%o cusp=%o time=%o\n",N,i,Cputime(t0);
  original_order:=Valuation(omega0,cusps[i]);
  ref:=1;
  ratios:=fixedratios;
  if original_order ne 0 then
   found:=false;
   for j in [2..g] do
    localomega:=om(Basis(Omega)[j]);
    if Valuation(localomega,cusps[i]) eq 0 then
     ref:=j;
     found:=true;
     break;
    end if;
   end for;
   assert found;
   ratios:=[om(b)/localomega:b in Basis(Omega)];
  end if;
  k:=0;
  printf "CUSP_DOUBLE_LOCAL_REFERENCE N=%o cusp=%o basis_index=%o original_order=%o\n",N,i,ref,original_order;
  maxmult:=D div degs[i];
  maxord:=2*maxmult;
  S,ex:=Completion(FunctionField(C),cusps[i]:Precision:=maxord+10);
  t:=S.1;
  expansions:=[t^k*ex(f):f in ratios];
  assert &and[AbsolutePrecision(s) ge maxord and
   (IsZero(s) or Valuation(s) ge 0):s in expansions];
  rows:=[];
  for s in expansions do
   row:=[];
   for j in [0..maxord-1] do
    cs:=Eltseq(Coefficient(s,j));
    assert #cs eq degs[i];
    row cat:=[GF(2)!c:c in cs];
   end for;
   Append(~rows,row);
  end for;
  A:=Matrix(GF(2),rows);
  for m in [1..maxmult] do
   B:=Submatrix(A,1,1,g,2*m*degs[i]);
   U:=Nullspace(B);
   assert Degree(U) eq g;
   localspaces[<i,m>]:=sub<W|[W!Eltseq(b):b in Basis(U)]>;
  end for;
 else
  for m in [1..D div degs[i]] do
   U,um:=DifferentialSpace(2*m*Divisor(cusps[i]));
   rows:=[W!Eltseq(Omega!um(b)):b in Basis(U)];
   localspaces[<i,m>]:=sub<W|rows>;
   assert Dimension(localspaces[<i,m>]) eq Dimension(U);
  end for;
 end if;
 printf "CUSP_DOUBLE_FAST_LOCAL N=%o cusp=%o/%o time=%o\n",N,i,#cusps,Cputime(t0);
end for;
bad:=[];
tested:=0;
for v in vs do
 if Limit gt 0 and tested ge Limit then break; end if;
 I:=W;
 for i in [1..#cusps] do
  if v[i] gt 0 then I meet:=localspaces[<i,v[i]>]; end if;
 end for;
 ell:=2*D+1-g+Dimension(I);
 assert ell ge 1;
 if tested lt DirectCheck or ell ne 1 then
  E:=&+[v[i]*Divisor(cusps[i]):i in [1..#cusps]];
  assert Degree(E) eq D;
  U:=RiemannRochSpace(2*E);
  assert ell eq Dimension(U);
 end if;
 if ell ne 1 then
  Append(~bad,v);
  printf "CUSP_DOUBLE_FAST_EXCEPTION N=%o index=%o dim=%o vector=%o\n",N,tested+1,ell,v;
 end if;
 tested+:=1;
 if tested mod 1000 eq 0 then
  printf "CUSP_DOUBLE_FAST_PROGRESS N=%o tested=%o exceptions=%o time=%o\n",N,tested,#bad,Cputime(t0);
 end if;
end for;
printf "CUSP_DOUBLE_FAST_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",
 N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);
quit;
