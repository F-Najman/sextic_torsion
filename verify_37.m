/* Complete the 37-primary case, using DKSS Proposition 1.4 for j=-9317.
   No assertion about the remaining sextic groups is made here. */
SetColumns(0);
SetSeed(20261021);
SetMemoryLimit(3*10^9);
Q := Rationals(); R<x> := PolynomialRing(Q);
A := -10395; B := 444150;
E0 := EllipticCurve([Q|0,0,0,A,B]);
assert jInvariant(E0) eq -9317;
f := x^6-3150*x^5+796635*x^4-75770100*x^3
      +3111596775*x^2-44606598750*x-85333003875;
assert IsIrreducible(f);
psi := DivisionPolynomial(E0,37);
assert psi mod f eq 0;
fac := Factorisation(psi);
printf "DIVISION_FACTOR_DEGREES %o\n", [<Degree(t[1]),t[2]> : t in fac];
small := [t[1] : t in fac | Degree(t[1]) le 6];
assert #small eq 3 and &and[Degree(g) eq 6 : g in small];
for ii in [1..#small] do
f := small[ii];
printf "CASE %o\nFIELD_POLYNOMIAL %o\n",ii,f;
K<w> := NumberField(f);
delta := w^3+A*w+B;
E := EllipticCurve([K|0,0,0,A*delta^2,B*delta^3]);
P := E![delta*w,delta^2,1];
assert P ne E!0 and 37*P eq E!0;
printf "FIELD_DEGREE %o\nPOINT_ORDER 37\n", Degree(K);
counts := [];
for p in PrimesUpTo(200) do
 if p le 3 then continue; end if;
 F := GF(p); RF<z> := PolynomialRing(F);
 fp := RF!f;
 if Discriminant(fp) eq 0 then continue; end if;
 for rt in Roots(fp) do
  r := rt[1]; dd := r^3+F!A*r+F!B;
  if dd eq 0 then continue; end if;
  Ep := EllipticCurve([F|0,0,0,F!A*dd^2,F!B*dd^3]);
  if Discriminant(Ep) eq 0 then continue; end if;
  nn := #Ep;
  assert nn mod 37 eq 0;
  Append(~counts,<p,Integers()!r,nn>);
  printf "GOOD_REDUCTION p=%o root=%o order=%o\n",p,r,nn;
  break;
 end for;
end for;
assert #counts ge 2;
primes := Seqset([c[1] : c in counts]);
for c in counts do
 primes join:= Seqset(PrimeDivisors(c[3]));
end for;
bound := 1;
for ell in Sort(Setseq(primes)) do
 vals := [Valuation(c[3],ell) : c in counts | c[1] ne ell];
 assert #vals gt 0;
 bound *:= ell^Minimum(vals);
end for;
printf "TORSION_ORDER_UPPER_BOUND %o\n",bound;
assert bound eq 37;
printf "FULL_TORSION C37\nCASE_DONE %o\n",ii;
end for;
printf "VERIFY_37_DONE\n";
