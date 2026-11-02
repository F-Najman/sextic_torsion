// Independent check of the C65 Cartier/formal-immersion certificate,
// built only from the equations printed in Lemma 9.4 of the paper.
SetColumns(0);
SetMemoryLimit(8*10^9);
t0 := Cputime();
F2 := GF(2);
P1<t> := FunctionField(F2);
R<Z> := PolynomialRing(P1);
j5 := (t^12+t^8+t^4+1)/(t^7+t^6+t^5);
num13 := Z^16+Z^4+1;
den13 := Z^28+Z^25+Z^24+Z^21+Z^20+Z^17+Z^16+Z^13;
pol := den13*j5 - num13;
pol := pol/LeadingCoefficient(pol);
printf "fibre product polynomial irreducible: %o\n", IsIrreducible(pol);
F0<z> := FunctionField(pol);
printf "genus Y0 = %o  (%o s)\n", Genus(F0), Cputime(t0);
A5 := (F0!t)^2/((F0!t)+1)^3;
A13 := (z^7+z^6+z^2+z)/(z^4+z+1);
S<W> := PolynomialRing(F0);
F<w> := FunctionField(W^2+W+A5+A13);
g := Genus(F);
printf "genus Y = %o  (%o s)\n", g, Cputime(t0);
jF := F!j5;
basepl := Zeros(P1!t) cat Poles(P1!t) cat Zeros(P1!(t^2+t+1));
cusps := [];
for p0 in basepl do
  for q0 in Decomposition(F0, p0) do
    cusps cat:= Decomposition(F, q0);
  end for;
end for;
cuspsU := []; for P in cusps do if not P in cuspsU then Append(~cuspsU, P); end if; end for; cusps := cuspsU;
printf "cusp degrees: %o\n", [Degree(P) : P in cusps];
printf "j valuations at cusps: %o\n", [Valuation(jF, P) : P in cusps];
printf "degree of j = %o  (%o s)\n", &+[-Valuation(jF,P)*Degree(P) : P in cusps], Cputime(t0);

// diamonds
tt := F!t; zz := F!z;
phT1 := hom< P1 -> F | 1/tt >;
phT0 := hom< F0 -> F | phT1, zz >;
tauT := hom< F -> F | phT0, w + tt/(tt+1) >;
phU1 := hom< P1 -> F | tt >;
phU0 := hom< F0 -> F | phU1, zz+1 >;
tauU := hom< F -> F | phU0, w + zz >;
// check that they respect the defining equations
polT := Polynomial([phT1(c) : c in Eltseq(DefiningPolynomial(F0))]);
polU := Polynomial([phU1(c) : c in Eltseq(DefiningPolynomial(F0))]);
printf "tauT respects eqs: %o %o\n", Evaluate(polT, tauT(zz)) eq 0, tauT(w)^2+tauT(w)+tauT(F!A5)+tauT(F!A13) eq 0;
printf "tauU respects eqs: %o %o\n", Evaluate(polU, tauU(zz)) eq 0, tauU(w)^2+tauU(w)+tauU(F!A5)+tauU(F!A13) eq 0;
printf "orders/commute: %o %o %o\n", tauT(tauT(tauT(tauT(w)))) eq w, tauU(tauU(tauU(tauU(w)))) eq w, tauT(tauU(w)) eq tauU(tauT(w));

V, h := DifferentialSpace(DivisorGroup(F)!0);
assert Dimension(V) eq g;
dt := Differential(tt);
B := [h(V.i) : i in [1..g]];
fs := [b/dt : b in B];
function pullback(sig, k)
  return sig(fs[k]) * Differential(sig(tt));
end function;
M2 := KMatrixSpace(F2, g, g);
MT := M2 ! &cat[Eltseq(pullback(tauT,k) @@ h) : k in [1..g]];
MU := M2 ! &cat[Eltseq(pullback(tauU,k) @@ h) : k in [1..g]];
MC := M2 ! &cat[Eltseq(Cartier(B[k]) @@ h) : k in [1..g]];
I := IdentityMatrix(F2, g);
printf "Cartier rank %o, rank(C-1) %o\n", Rank(MC), Rank(MC-I);
Ssp := RowSpace(MT-I) + RowSpace(MU-I);
printf "dim diamond image = %o\n", Dimension(Ssp);
U := sub< VectorSpace(F2,g) | [b*(MC-I) : b in Basis(Ssp)] >;
printf "dim U = %o  (%o s)\n", Dimension(U), Cputime(t0);

// kernels of restriction: coefficient spaces of H^0(Omega(-mP))
Wsp := AssociativeArray();
for i in [1..#cusps] do
  P := cusps[i];
  for m in [1..(6 div Degree(P))] do
    Vm, hm := DifferentialSpace(m*P);
    Wsp[<i,m>] := sub< VectorSpace(F2,g) | [Vector(F2, Eltseq(hm(Vm.k) @@ h)) : k in [1..Dimension(Vm)]] >;
  end for;
end for;
printf "kernel spaces computed (%o s)\n", Cputime(t0);
abovez01 := [ Valuation(zz, P) gt 0 or Valuation(zz+1, P) gt 0 : P in cusps ];
printf "quartic cusps above z=0,1: %o\n", [i : i in [1..#cusps] | Degree(cusps[i]) eq 4 and abovez01[i]];

// enumerate all effective cuspidal divisors of degree 6
degs := [Degree(P) : P in cusps];
procedure rec(i, rem, cur, ~list)
  if rem eq 0 then Append(~list, cur); return; end if;
  if i gt #cusps then return; end if;
  for m in [0..(rem div degs[i])] do
    rec(i+1, rem - m*degs[i], cur cat (m gt 0 select [<i,m>] else []), ~list);
  end for;
end procedure;
divs := [];
rec(1, 6, [], ~divs);
printf "total degree-6 cuspidal divisors: %o\n", #divs;
admissible := function(D)
  for pr in D do
    i := pr[1]; m := pr[2];
    if degs[i] eq 2 and m mod 3 ne 0 then return false; end if;
    if degs[i] eq 4 and not abovez01[i] and m mod 3 ne 0 then return false; end if;
  end for;
  return true;
end function;
nadm := 0; badadm := 0; badall := 0;
for D in divs do
  K := VectorSpace(F2,g);
  for pr in D do K meet:= Wsp[pr]; end for;
  r := Dimension(U) - Dimension(U meet K);
  adm := admissible(D);
  if adm then nadm +:= 1; end if;
  if r ne 6 then
    badall +:= 1;
    if adm then badadm +:= 1; printf "ADMISSIBLE FAILURE %o rank %o\n", D, r; end if;
  end if;
end for;
printf "admissible: %o, admissible failures: %o, failures among all: %o\n", nadm, badadm, badall;
printf "C65_CHECK_DONE %o s\n", Cputime(t0);
quit;
