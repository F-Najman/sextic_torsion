// Independent check of the Cartier formal-immersion certificates for
// Y63 = X1(63)/<8> (Case:="63") and Y51 = X1(51)/<16> (Case:="51").
// Models are built from scratch: universal Tate normal form over X1(9)
// (resp. X1(17), computed here from x(8P)=x(9P)), then a root of the
// 7- (resp. 3-) division polynomial is adjoined.
SetColumns(0);
SetMemoryLimit(12*10^9);
if not assigned Case then Case := "63"; end if;
t0 := Cputime();
F2 := GF(2);
if Case eq "63" then
  P1<t> := FunctionField(F2);
  b := t^2*(t-1)*(t^2-t+1); c := t^2*(t-1);
  E := EllipticCurve([1-c, -b, -b, 0, 0]);
  P := E![0,0];
  assert 9*P eq E!0 and 3*P ne E!0;
  psi := DivisionPolynomial(E, 7);
  base := P1; ell := 7;
else
  // X1(17) in char 2 from the Tate normal form b=rs(r-1), c=s(r-1)
  A2<r,s> := PolynomialRing(F2, 2);
  Q2 := FieldOfFractions(A2);
  bb := Q2!(r*s*(r-1)); cc := Q2!(s*(r-1));
  EE := EllipticCurve([1-cc, -bb, -bb, 0, 0]);
  PP := EE![0,0];
  P8 := 8*PP; P9 := 9*PP;
  rel := Numerator(P8[1]/P8[3] - P9[1]/P9[3]);
  fac := Factorization(rel);
  printf "factors of x(8P)-x(9P): %o\n", [<TotalDegree(f[1]), f[2]> : f in fac];
  // the X1(17) factor: genus 5
  found := false;
  for fpair in fac do
    f := fpair[1];
    if Degree(f, s) eq 0 or Degree(f, r) eq 0 then continue; end if;
    P1<x1> := FunctionField(F2);
    RS<S1> := PolynomialRing(P1);
    g := Evaluate(f, [x1, S1]);
    if not IsIrreducible(g) then continue; end if;
    Ftry := FunctionField(g/LeadingCoefficient(g));
    gen := Genus(Ftry);
    printf "factor of bidegree (%o,%o): genus %o\n", Degree(f,r), Degree(f,s), gen;
    if gen eq 5 then FX := Ftry; found := true; break; end if;
  end for;
  assert found;
  rr := FX!x1; ss := FX.1;
  b := rr*ss*(rr-1); c := ss*(rr-1);
  E := EllipticCurve([1-c, -b, -b, 0, 0]);
  P := E![0,0];
  assert 17*P eq E!0 and P ne E!0;
  psi := DivisionPolynomial(E, 3);
  base := FX; ell := 3;
end if;
psi := psi/LeadingCoefficient(psi);
printf "psi_%o: degree %o, irreducible %o, separable %o\n", ell, Degree(psi), IsIrreducible(psi), IsSeparable(psi);
F<u> := FunctionField(psi);
g := Genus(F);
printf "genus %o  (%o s)\n", g, Cputime(t0);
jF := F!jInvariant(E);
printf "degree of j: %o\n", Degree(jF);
cusps := Poles(jF);
degs := [Degree(Q) : Q in cusps];
printf "cusp degrees: %o\n", degs;

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
B := [h(V.i) : i in [1..g]];
M2 := KMatrixSpace(F2, g, g);
MC := M2 ! &cat[Eltseq(Cartier(B[k]) @@ h) : k in [1..g]];
I := IdentityMatrix(F2, g);
U := RowSpace(MC - I);
printf "dim U = image of (Cartier-1): %o  (%o s)\n", Dimension(U), Cputime(t0);

Wsp := AssociativeArray();
for i in [1..#cusps] do
  for m in [1..(6 div degs[i])] do
    Vm, hm := DifferentialSpace(m*cusps[i]);
    Wsp[<i,m>] := sub< VectorSpace(F2,g) | [Vector(F2, Eltseq(hm(Vm.k) @@ h)) : k in [1..Dimension(Vm)]] >;
  end for;
end for;
printf "kernel spaces: %o (%o s)\n", #Keys(Wsp), Cputime(t0);

procedure rec(i, rem, cur, ~list)
  if rem eq 0 then Append(~list, cur); return; end if;
  if i gt #cusps then return; end if;
  for m in [0..(rem div degs[i])] do
    rec(i+1, rem - m*degs[i], cur cat (m gt 0 select [<i,m>] else []), ~list);
  end for;
end procedure;
divs := [];
rec(1, 6, [], ~divs);
printf "degree-6 cuspidal divisors: %o\n", #divs;
bad := 0;
for D in divs do
  K := VectorSpace(F2,g);
  for pr in D do K meet:= Wsp[pr]; end for;
  rk := Dimension(U) - Dimension(U meet K);
  if rk ne 6 then bad +:= 1; if bad le 20 then printf "FAILURE %o rank %o\n", D, rk; end if; end if;
end for;
printf "CARTIER_CHECK_DONE case %o: divisors %o failures %o (%o s)\n", Case, #divs, bad, Cputime(t0);
quit;
