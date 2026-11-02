/* Exact A3 quotient model and oriented finite-place map.
   Load after the v2 mdmagma spec and cyclic common-zero inverse.
   Library without a log of its own: it is loaded by filter_a3_m2_30.m and test_a3_m2_30_map.m,
   and the genus 9 of W is printed in logs/n30/a3_filter.log (marker A3_M2_30_MODEL_READY). */
A3p:=7;A3X15:=MDX1(15,GF(A3p));A3C15:=Curve(A3X15);
A3F15:=FunctionField(A3C15);
A3A15,A3ma:=AlgorithmicFunctionField(A3F15);
assert Genus(A3A15) eq 1;
A3E15:=EllipticCurve([A3ma(a):a in A3X15`_E]);
A3P15:=A3E15![A3ma(a):a in A3X15`_P];assert Order(A3P15) eq 15;
A3disc15:=Discriminant(A3E15);A3j15:=jInvariant(A3E15);
A3R<A3v>:=PolynomialRing(A3A15);A3f:=A3v^2-A3disc15;
assert IsIrreducible(A3f) and GCD(A3f,Derivative(A3f)) eq 1;
A3K<A3root>:=FunctionField(A3f);assert Degree(A3K) eq 2 and Genus(A3K) eq 9;
function A3PlacesAboveBase(pp)
 aa,bb:=TwoGenerators(pp);
 raw:=Setseq(Seqset(Decomposition(A3K,pp)));
 // In a tower, Decomposition may return places over other intermediate
 // places above the same rational base place. Select the actual base
 // place using both ideal generators before using residue degrees.
 return [qq:qq in raw|Valuation(A3K!aa,qq) gt 0 and Valuation(A3K!bb,qq) gt 0];
end function;
A3basecs:=Poles(A3j15);A3cs:=[];
for pp in A3basecs do
 assert IsOdd(Valuation(A3disc15,pp));
 qq:=A3PlacesAboveBase(pp);assert #qq eq 1 and Degree(qq[1]) eq Degree(pp);
 assert Valuation(A3K!A3j15,qq[1]) lt 0;Append(~A3cs,qq[1]);
end for;
assert #Seqset(A3cs) eq #A3cs and {*Degree(q):q in A3cs*} eq {*1^^8,4^^2*};
assert &+[-Valuation(A3K!A3j15,q)*Degree(q):q in A3cs] eq 192;
// Every base cusp is totally ramified, with unchanged residue field.
// Consequently each rational base cusp orbit has its unique rational
// orbit upstairs. This is exactly one half its conorm divisor.
A3baseorbits:=CuspOrbitsQ(A3X15);A3orbits:=[];
for OO in A3baseorbits do
 DD:=FunctionFieldDivisor(OO);assert FunctionField(DD) eq A3A15;
 ss,mm:=Support(DD);UU:=DivisorGroup(A3K)!0;
 for i in [1..#ss] do
  qq:=A3PlacesAboveBase(ss[i]);assert #qq eq 1 and Degree(qq[1]) eq Degree(ss[i]);
  UU+:=mm[i]*Divisor(qq[1]);
 end for;
 assert Degree(UU) eq Degree(OO);Append(~A3orbits,UU);
end for;
assert Sort([Degree(OO):OO in A3orbits]) eq [1,1,1,1,2,2,4,4];
A3baseindex:=Index([Degree(OO):OO in A3orbits],1);
A3c0:=Support(A3orbits[A3baseindex])[1];assert Degree(A3c0) eq 1;
A3coordinates:=[A3K!A3ma(A3F15.1),A3K!A3ma(A3F15.2),A3root];
A3fields:=[GF(A3p^d):d in [1..6]];

function A3FrobeniusTupleKey(vals,common,dd)
 L:=Parent(vals[1]);assert &and[Parent(v) eq L:v in vals];
 assert Characteristic(L) eq A3p and dd mod Degree(L) eq 0 and #common eq A3p^dd;
 if L cmpeq common then vv:=vals;
 elif Degree(L) eq 1 then vv:=[common!(Integers()!a):a in vals];
 else
  ff:=DefiningPolynomial(L,PrimeField(L));assert Degree(ff) eq Degree(L);
  fk:=PolynomialRing(common)![common!(Integers()!a):a in Coefficients(ff)];
  rr:=Roots(fk);assert #rr eq Degree(L);hh:=hom<L->common|rr[1][1]>;
  vv:=[hh(a):a in vals];
 end if;
 keys:=[];
 for i in [1..dd] do
  key:=&cat[[Integers()!a:a in Eltseq(v,GF(A3p))]:v in vv];
  assert #key eq #vals*dd;Append(~keys,key);vv:=[v^A3p:v in vv];
 end for;
 return Sort(keys)[1];
end function;

function A3OrientedValues(E,P2,Q30)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));
 assert a1 eq 0 and a3 eq 0 and a6 eq 0;
 assert Order(P2) eq 2 and Order(Q30) eq 30 and P2 ne 15*Q30;
 P:=2*Q30;assert Order(P) eq 15;
 e1:=P2[1];e2:=(15*Q30)[1];e3:=-a2-e1-e2;
 assert e1 ne e2 and e2 ne e3 and e3 ne e1;
 rr:=4*(e1-e2)*(e2-e3)*(e3-e1);assert rr^2 eq Discriminant(E);
 xx:=P[1];yy:=P[2];assert yy ne 0;
 slope:=(3*xx^2+2*a2*xx+a4)/(2*yy);
 A:=3*xx+a2-slope^2;assert A ne 0;uu:=2*yy/A;
 bc:=TateNormalForm_bc(E,P);ET:=EllipticCurve([1-bc[2],-bc[1],-bc[1],0,0]);
 assert Discriminant(ET)*uu^12 eq Discriminant(E);
 oriented:=rr/uu^6;assert oriented^2 eq Discriminant(ET);
 xy:=TateNormalForm_xy(E,P);
 return xy cat [oriented];
end function;

function A3PlaceFromLevel(E,P2,Q30)
 dd:=Degree(BaseRing(E));assert dd in [1..6];
 vals:=A3OrientedValues(E,P2,Q30);
 key:=A3FrobeniusTupleKey(vals,A3fields[dd],dd);
 P:=2*Q30;
 y:=ModuliPointCyclicCommonZeros(A3X15,E,rec<recformat<P:PtEll>|P:=P>);
 dy:=Degree(y);assert dd mod dy eq 0 and not IsCusp(A3X15,y);
 baseplace:=FunctionFieldPlace(y);assert FunctionField(baseplace) eq A3A15;
 candidates:=A3PlacesAboveBase(baseplace);matched:=[];
 for q in candidates do
  if dd mod Degree(q) ne 0 then continue;end if;
  vv:=[Evaluate(f,q):f in A3coordinates];
  if A3FrobeniusTupleKey(vv,A3fields[dd],dd) eq key then Append(~matched,q);end if;
 end for;
 assert #matched eq 1;
 q:=matched[1];assert dd mod Degree(q) eq 0 and q notin A3cs;
 return q,dd div Degree(q),key;
end function;

printf "A3_M2_30_MODEL_READY genus=9 jdegree=192 rational_cusp_orbit_degrees=%o finite_cusp_degrees=%o\n",
 [Degree(OO):OO in A3orbits],[Degree(q):q in A3cs];
