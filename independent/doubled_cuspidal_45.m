// Independent check of the C45 doubled-cuspidal certificate:
// for every effective cuspidal divisor C of degree 6 on X1(45)_F2,
// h^0(2C) = 1, i.e. restriction of H^0(Omega) to 2C has rank 12.
// Model: universal Tate normal form over X1(9), adjoin a point of order 5.
SetColumns(0);
SetMemoryLimit(12*10^9);
t0 := Cputime();
F2 := GF(2);
P1<t> := FunctionField(F2);
b := t^2*(t-1)*(t^2-t+1); c := t^2*(t-1);
E9 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 9*E9![0,0] eq E9!0 and 3*E9![0,0] ne E9!0;
psi := DivisionPolynomial(E9, 5);
psi := psi/LeadingCoefficient(psi);
printf "psi_5 degree %o irreducible %o\n", Degree(psi), IsIrreducible(psi);
F1<x5> := FunctionField(psi);
a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E9)]);
S<Y> := PolynomialRing(F1);
q := Y^2 + a1*x5*Y + a3*Y - (x5^3 + a2*x5^2 + a4*x5 + a6);
printf "y-equation irreducible %o\n", IsIrreducible(q);
F<y5> := FunctionField(q);
EF := BaseChange(E9, F);
P5 := EF![F!x5, y5];
P9 := EF![0,0];
assert 5*P5 eq EF!0;
g := Genus(F);
printf "genus %o (%o s)\n", g, Cputime(t0);
jF := F!jInvariant(E9);
// collect poles of j by decomposing the places of P1 where j has a pole
jP1 := jInvariant(E9);
basepl := [p : p in Support(PoleDivisor(Divisor(jP1)))];
cusps := [];
for p0 in basepl do
  for q0 in Decomposition(F1, p0) do
    cusps cat:= Decomposition(F, q0);
  end for;
end for;
cuspsU := []; for P in cusps do if not P in cuspsU then Append(~cuspsU, P); end if; end for; cusps := cuspsU;
degs := [Degree(P) : P in cusps];
printf "cusp degrees: %o\n", Sort(degs);
printf "degree of j = %o\n", &+[-Valuation(jF,P)*Degree(P) : P in cusps];
small := [i : i in [1..#cusps] | degs[i] le 6];
printf "cusps of degree <= 6: %o\n", Sort([degs[i] : i in small]);

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
Wsp := AssociativeArray();
for i in small do
  for m in [1..(12 div degs[i])] do
    Vm, hm := DifferentialSpace(m*cusps[i]);
    Wsp[<i,m>] := sub< VectorSpace(F2,g) | [Vector(F2, Eltseq(hm(Vm.k) @@ h)) : k in [1..Dimension(Vm)]] >;
  end for;
end for;
printf "kernel spaces: %o (%o s)\n", #Keys(Wsp), Cputime(t0);
procedure rec(j, rem, cur, ~list)
  if rem eq 0 then Append(~list, cur); return; end if;
  if j gt #small then return; end if;
  i := small[j];
  for m in [0..(rem div degs[i])] do
    rec(j+1, rem - m*degs[i], cur cat (m gt 0 select [<i,m>] else []), ~list);
  end for;
end procedure;
divs := [];
rec(1, 6, [], ~divs);
printf "degree-6 cuspidal divisors: %o\n", #divs;
bad := 0;
for D in divs do
  K := VectorSpace(F2,g);
  for pr in D do K meet:= Wsp[<pr[1], 2*pr[2]>]; end for;
  if g - Dimension(K) ne 12 then bad +:= 1; if bad le 20 then printf "FAILURE %o codim %o\n", D, g-Dimension(K); end if; end if;
end for;
printf "C45_DOUBLED_CHECK_DONE divisors %o failures %o (%o s)\n", #divs, bad, Cputime(t0);
quit;
