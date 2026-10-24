// Propositions 10.1-10.4 (the sporadic curves): fields, torsion, uniqueness statements,
// the three 2-isogenous curves of E_{2,24}, and the facts in the last sentence of Prop. 10.3.
// Run from the repository root:  magma -n -b verify_sporadic_curves.m
// Closing marker: SPORADIC_DONE
SetColumns(0); SetMemoryLimit(8*10^9);
R<x> := PolynomialRing(Rationals());

// ---- Z/37 (Proposition 10.1) ----
E0 := EllipticCurve([-10395, 444150]);
assert jInvariant(E0) eq -9317 and Conductor(E0) eq 1225;
assert IsIsomorphic(E0, EllipticCurve([1,1,1,-8,6]));           // 1225.b2
fac := Factorization(DivisionPolynomial(E0, 37));
assert Sort([Degree(f[1]) : f in fac]) eq [6,6,6,666] and &and[f[2] eq 1 : f in fac];
f37 := x^6-3150*x^5+796635*x^4-75770100*x^3+3111596775*x^2-44606598750*x-85333003875;
assert exists{f : f in fac | f[1] eq f37};
K37<w> := NumberField(f37);
assert IsNormal(K37) and IsCyclic(GaloisGroup(K37)) and Signature(K37) eq 6;
assert Discriminant(Integers(K37)) eq 300125;                    // 6.6.300125.1
for f in fac do if Degree(f[1]) eq 6 then assert #Roots(f[1], K37) eq 6; end if; end for;
delta := w^3-10395*w+444150;
E37 := EllipticCurve([-10395*delta^2, 444150*delta^3]);
P37 := E37![delta*w, delta^2];
assert Order(P37) eq 37;
curves := [];
for f in fac do if Degree(f[1]) eq 6 then
  for r in Roots(f[1], K37) do
    d := r[1]^3-10395*r[1]+444150;
    Append(~curves, EllipticCurve([-10395*d^2, 444150*d^3]));
  end for;
end if; end for;
assert #curves eq 18 and &and[IsIsomorphic(C, E37) : C in curves];
for p in [29, 41] do
  Fp := GF(p);
  w0 := [r[1] : r in Roots(ChangeRing(f37, Fp))][1];
  d0 := w0^3-10395*w0+444150;
  assert #EllipticCurve([-10395*d0^2, 444150*d0^3]) eq 37;
end for;
print "Z37: one field (cyclic, 6.6.300125.1), one curve up to isomorphism, torsion order 37";

// ---- Z/25 and Z/4 x Z/12 (Proposition 10.2) ----
L<a> := NumberField(x^3-x^2+1); S<z> := PolynomialRing(L);
K1 := ext<L | z^2+(a^2-2*a-1)*z+a>; K25 := AbsoluteField(K1);
aa := K25!(K1!a); bb := K25!(K1.1);
r := (aa^2*bb-aa*bb+bb-1)/(aa*(aa*bb-1)); s := (aa*bb-bb+1)/(aa*bb);
E25 := EllipticCurve([1-s*(r-1), -r*s*(r-1), -r*s*(r-1), 0, 0]);
assert Order(E25![0,0]) eq 25 and Invariants(TorsionSubgroup(E25)) eq [25];
assert Discriminant(Integers(K25)) eq 28037 and Signature(K25) eq 2;   // 6.2.28037.1
print "Z25: torsion [25] over 6.2.28037.1";
L<al> := NumberField(x^3-3*x-4);
K412 := AbsoluteField(ext<L | PolynomialRing(L).1^2+1>);
E412 := BaseChange(EllipticCurve([1,-1,1,4,-1]), K412);               // 162.d2
assert Invariants(TorsionSubgroup(E412)) eq [4,12];
assert IsIsomorphic(K412, NumberField(x^6-2*x^3+9*x^2+6*x+2));       // 6.0.419904.2
print "Z4xZ12: torsion [4,12] over 6.0.419904.2";

// ---- Z/2 x Z/24 (Propositions 10.3 and 10.4) ----
K<w> := NumberField(x^6-2*x^4+2*x^2-3);                               // 6.2.231852.1
assert Discriminant(Integers(K)) eq 231852 and Signature(K) eq 2;
b := (3*w^5+7*w^4-24*w^3-15*w^2+21*w+18)/32;
c := (-7*w^5+5*w^4+24*w^3-5*w^2-17*w-2)/16;
E := EllipticCurve([0, c, 0, (1-b-c)*b, 0]);
P := E![0,0]; Q := E![b,b];
assert Order(P) eq 2 and Order(Q) eq 24 and P ne 12*Q;
T, mT := TorsionSubgroup(E);
assert Invariants(T) eq [2,24];
A := Automorphisms(K);
assert #A eq 2;
for sig in A do
  Es := EllipticCurve([sig(cc) : cc in aInvariants(E)]);
  assert IsIsomorphic(Es, E);
end for;
assert not HasComplexMultiplication(E);
mj := MinimalPolynomial(jInvariant(E));
assert Degree(mj) eq 3;
k := NumberField(mj);
assert Discriminant(Integers(OptimizedRepresentation(k))) eq -139;
assert Norm(Conductor(E)) eq 6;
print "Z2xZ24: torsion [2,24]; #Aut(K)=2 and E^sigma ~ E; j cubic (disc -139); conductor norm 6";
// the two other models found by the first computation give the same curve
fs := [x^6-3*x^5+1/2*x^4+4*x^3-x^2-3/2*x+1/2, x^6+4*x^5+41/6*x^4+31/6*x^3+4/3*x^2-1/3*x-1/6];
hs := [3-3*x-13*x^2+2*x^3+10*x^4-4*x^5, (-92-15*x+542*x^2+1074*x^3+744*x^4+216*x^5)/17];
for i in [1,2] do
  Ki<v> := NumberField(fs[i]); u := Evaluate(hs[i], v);
  q := (u+v+2)/(u-v);
  t := (-u^3-u^2*v-2*u^2+u*v^2+v^3+2*v^2)/(u^3+u^2*v+2*u^2+u*v^2+2*u*v+v^3+2*v^2-2);
  Ei := EllipticCurve([0, t^2-2*q*t-2, 0, -(t^2-1)*(q*t+1)^2, 0]);
  iso, phi := IsIsomorphic(Ki, K); assert iso;
  assert IsIsomorphic(EllipticCurve([phi(cc) : cc in aInvariants(Ei)]), E);
end for;
print "Z2xZ24: both points of the Derickx-Sutherland model give this curve";
// level structures: 16 up to sign
pts := [mT(g) : g in T];
lev := #[<p,q> : p in pts, q in pts | Order(p) eq 2 and Order(q) eq 24 and p ne 12*q] div 2;
assert lev eq 16;
// the three 2-isogenous curves (Proposition 10.4)
R2 := [p : p in pts | Order(p) eq 2];
tors := [];
for p in R2 do
  E2 := IsogenyFromKernel(E, PolynomialRing(K).1 - p[1]);
  Append(~tors, Invariants(TorsionSubgroup(E2)));
end for;
assert Sort(tors) eq Sort([[24],[2,12],[24]]);
print "Z48: the 2-isogenous curves have torsion", tors;
print "SPORADIC_DONE";
quit;
