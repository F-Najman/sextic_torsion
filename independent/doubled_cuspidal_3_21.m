// Independent check of the C3 x C21 doubled-cuspidal certificate:
// on the fixed-pairing component X of X1(3,21) over F4, h^0(2C)=1 for every
// effective cuspidal divisor C of degree 3 (equivalently, restriction of
// H^0(Omega) to 2C has rank 6).
// Model: the 3-division field of the universal Tate normal form over X1(7)_F4.
SetColumns(0);
SetMemoryLimit(12*10^9);
t0 := Cputime();
F4<zeta> := GF(4);
P1<t> := FunctionField(F4);
b := t^3-t^2; c := t^2-t;
E7 := EllipticCurve([1-c, -b, -b, 0, 0]);
assert 7*E7![0,0] eq E7!0;
psi := DivisionPolynomial(E7, 3); psi := psi/LeadingCoefficient(psi);
printf "psi3 irreducible %o\n", IsIrreducible(psi);
F1<x3> := FunctionField(psi);
a1,a2,a3,a4,a6 := Explode([F1!a : a in aInvariants(E7)]);
S1<Y> := PolynomialRing(F1);
qy := Y^2 + a1*x3*Y + a3*Y - (x3^3 + a2*x3^2 + a4*x3 + a6);
printf "y-eq irreducible %o\n", IsIrreducible(qy);
F2f<y3> := FunctionField(qy);
S2<X> := PolynomialRing(F2f);
psiF := Polynomial([F2f!F1!c : c in Eltseq(psi)]);
cub := psiF div (X - F2f!x3);
fc := Factorization(cub);
printf "remaining cubic factors: %o\n", [Degree(f[1]) : f in fc];
assert #fc eq 1 and Degree(fc[1][1]) eq 3;
F<xp> := FunctionField(fc[1][1]);
g := Genus(F);
printf "genus %o (%o s)\n", g, Cputime(t0);
printf "constant field %o\n", ConstantField(F);
// cusps: places above the poles of j on P1
jP1 := jInvariant(E7);
basepl := [p : p in Support(PoleDivisor(Divisor(jP1)))];
cusps := [];
for p0 in basepl do
  for q1 in Decomposition(F1, p0) do
    for q2 in Decomposition(F2f, q1) do
      cusps cat:= Decomposition(F, q2);
    end for;
  end for;
end for;
cuspsU := []; for P in cusps do if not P in cuspsU then Append(~cuspsU, P); end if; end for; cusps := cuspsU;
degs := [Degree(P) : P in cusps];
printf "cusp degrees: %o (total %o)\n", Sort(degs), &+degs;
V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
small := [i : i in [1..#cusps] | degs[i] le 3];
Wsp := AssociativeArray();
for i in small do
  for m in [1..(6 div degs[i])] do
    Vm, hm := DifferentialSpace(m*cusps[i]);
    Wsp[<i,m>] := sub< VectorSpace(F4,g) | [Vector(F4, Eltseq(hm(Vm.k) @@ h)) : k in [1..Dimension(Vm)]] >;
  end for;
end for;
printf "kernel spaces %o (%o s)\n", #Keys(Wsp), Cputime(t0);
procedure rec(j, rem, cur, ~list)
  if rem eq 0 then Append(~list, cur); return; end if;
  if j gt #small then return; end if;
  i := small[j];
  for m in [0..(rem div degs[i])] do
    rec(j+1, rem - m*degs[i], cur cat (m gt 0 select [<i,m>] else []), ~list);
  end for;
end procedure;
divs := [];
rec(1, 3, [], ~divs);
printf "degree-3 cuspidal divisors: %o\n", #divs;
bad := 0;
for D in divs do
  K := VectorSpace(F4,g);
  for pr in D do K meet:= Wsp[<pr[1], 2*pr[2]>]; end for;
  if g - Dimension(K) ne 6 then bad +:= 1; if bad le 20 then printf "FAILURE %o codim %o\n", D, g-Dimension(K); end if; end if;
end for;
printf "C3C21_DOUBLED_CHECK_DONE divisors %o failures %o (%o s)\n", #divs, bad, Cputime(t0);
quit;
