/* sieve_M2_strong.m -- degree-D Hecke sieve on X_1(N) (M=1) or X_1(2,N) (M=2) over F_p. In this repository
   it is run with D:=6 for the sextic paper: jobs.txt gives the parameters of each run, and
   README.md the statements it supports.
   Variant for M=2 that uses the sufficient cusp condition q in <p> modulo 2N, valid for every
   cuspidal complement (used for the first run on X_1(2,22)).
   It implements the Hecke sieve of Derickx-Najman (J. reine angew. Math., Proposition 5.5;
   Proposition 5.4 in arXiv:2412.16016v2). For a non-cuspidal part D' of degree d' < D the
   cuspidal complement C has degree k = D-d', and [A_q C] = 0 is guaranteed if q lies in <p>.H_c
   for every Galois orbit c of cusps whose F_p-reduction has degree <= k, where
   H_c = Gal(Q(zeta_N)/F_c) and F_c is the field of definition of c.
   Hypotheses checked elsewhere: rk J(Q) = 0, and the totally cuspidal reduction is excluded by
   Proposition 5.1 of the paper.
   The fast Hecke operators of fast_hecke.m are compared with the reference correspondence in
   check_char2_hecke.m and verify_35.m, and the X_1(2,N) routines of fast_hecke_11.m in
   test_fast_hecke_*.m and test_fast_diamonds_22.m.
   Usage (example): magma -n -b N:=22 p:=3 Qs:=5,7 D:=6 M:=2 sieve_M2_strong.m
   Log markers: SETUP CUSPS CUSPFIELDS CC NONCUSP REPS QLIST TYPE PROGRESS SURVIVOR SURVIVOR_DATA
                RESULT SUMMARY SIEVE_DONE ERROR */
SetColumns(0);
if not assigned N then N := 26; elif Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned p then p := 7; elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned D then D := 5; elif Type(D) eq MonStgElt then D := StringToInteger(D); end if;
if not assigned M then M := 1; elif Type(M) eq MonStgElt then M := StringToInteger(M); end if;
if not assigned Qs then Qs := [3,5,7,11,13]; elif Type(Qs) eq MonStgElt then Qs := [StringToInteger(s) : s in Split(Qs, ",")]; end if;
if not assigned MemGB then MemGB := 8; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned Shard then Shard := 0; elif Type(Shard) eq MonStgElt then Shard := StringToInteger(Shard); end if;
if not assigned NShard then NShard := 1; elif Type(NShard) eq MonStgElt then NShard := StringToInteger(NShard); end if;
SetMemoryLimit(MemGB*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");   // v6: fast Hecke operators and diamond orbits (X_1(N) only)
function HeckeOp(X, q, D) if M eq 1 then return HeckeOperatorFast(X, q, D); else return HeckeOperator(X, q, D); end if; end function;
function PUD(X, S) if #S eq 0 then return S; elif M eq 1 then return PlacesUpToDiamondFast(X, S); else return PlacesUpToDiamond(X, S); end if; end function;
t0 := Cputime();
if M eq 1 then X := MDX1(N, GF(p)); else X := MDX11(2, N, GF(p)); end if;
printf "SETUP N=%o M=%o p=%o genus=%o time=%o\n", N, M, p, Genus(X), Cputime(t0);
cusps := Cusps(X);
printf "CUSPS number=%o degrees=%o\n", #cusps, Sort([Degree(c) : c in cusps]);

// Strong sufficient cusp condition: all cusps lie in Q(zeta_(2N)).
// q in <p> modulo 2N kills every cuspidal complement class after reduction.
assert M eq 2;
lev := 2*N; Zlev := Integers(lev);
Hp := {Zlev!1}; hh := Zlev!p;
while hh notin Hp do Include(~Hp,hh); hh *:= p; end while;
Qs := [q : q in Qs | IsPrime(q) and q ne p and (2*N) mod q ne 0];
QsCC := [[q : q in Qs | Zlev!q in Hp] : k in [1..D-1]];
printf "STRONG_CUSP_CONDITION modulus=%o subgroup=%o usable=%o\n",lev,Hp,QsCC[1];
printf "QLIST all=%o\n",Qs;

// Which degrees d can carry non-cuspidal places?  A non-cuspidal place of degree d gives an
// elliptic curve over F_{p^d} with a point of order N (M=1) or a subgroup Z/2 x Z/N (M=2).
// Enumerate all E/F_{p^d} (all j, all twists) -- unconditional -- and skip empty degrees.
function DegreeHasPoints(d)
  F := GF(p^d);
  for j in F do
    E := EllipticCurveFromjInvariant(j);
    for Et in Twists(E) do
      inv := Invariants(AbelianGroup(Et)); if #inv eq 0 then inv := [1]; end if;
      if M eq 1 then
        if inv[#inv] mod N eq 0 then return true; end if;
      else
        if #inv eq 2 and inv[2] mod N eq 0 and inv[1] mod 2 eq 0 then return true; end if;
      end if;
    end for;
  end for;
  return false;
end function;
hasdeg := [ (p^d le 400000) select DegreeHasPoints(d) else true : d in [1..D] ];
printf "DEGREES with points of the required type over F_p^d (d=1..%o): %o  time=%o\n", D, hasdeg, Cputime(t0);
places := [ hasdeg[d] select NoncuspidalPlaces(X, d) else [] : d in [1..D] ];
printf "NONCUSP counts deg 1..%o: %o  time=%o\n", D, [#places[d] : d in [1..D]], Cputime(t0);
reps := [ PUD(X, places[d]) : d in [1..D] ];
printf "REPS up to diamond: %o  time=%o\n", [#reps[d] : d in [1..D]], Cputime(t0);

function DivisorsOfType(lambda)
  d1 := lambda[1]; rest := lambda[2..#lambda];
  degs := Sort(Setseq(Seqset(rest)));
  result := [];
  for e in degs do if #places[e] eq 0 then return result; end if; end for;
  if #reps[d1] eq 0 then return result; end if;
  choices := [* *];
  for e in degs do
    k := #[x : x in rest | x eq e];
    Append(~choices, [ &+[Divisor(x) : x in S] : S in Multisets(Seqset(places[e]), k) ]);
  end for;
  for R in reps[d1] do
    partial := [Divisor(R)];
    for ch in choices do partial := [ a + b : a in partial, b in ch ]; end for;
    result cat:= partial;
  end for;
  return result;
end function;

function Describe(Dv)
  a, b := Support(Dv);
  return [< b[i], Degree(a[i]), jInvariant(X, a[i]) > : i in [1..#a]];
end function;
// reconstructible description of a divisor: for each place [* deg, mult, defpoly, coords *]
// where defpoly is the defining polynomial (coefficient list over F_p) of the field K of a
// representative point and coords the Eltseq's of its projective coordinates in K;
// reconstruct with K := ext<GF(p)|Polynomial(GF(p),defpoly)>; Places(Curve(X)(K)![K!c : c in coords]).
function DescribeData(Dv)
  a, b := Support(Dv);
  out := [* *];
  for i in [1..#a] do
    P := a[i]; dg := Degree(P);
    try
      pt := RepresentativePoint(P);
      K := Ring(Parent(pt));
      defpoly := [Integers()!c : c in Eltseq(DefiningPolynomial(K, GF(p)))];
      coords := [[Integers()!e : e in Eltseq(c, GF(p))] : c in Eltseq(pt)];
      Append(~out, [* dg, b[i], defpoly, coords *]);
    catch e
      Append(~out, [* dg, b[i], "ERROR", Sprint(e`Object) *]);
    end try;
  end for;
  return out;
end function;

allsurv := [];
for dd in [1..D] do
  qlist := (dd eq D) select Qs else QsCC[D-dd];
  for lambda in Partitions(dd) do
    divs := DivisorsOfType(lambda);
    if #divs eq 0 then printf "TYPE %o: no divisors\n", lambda; continue; end if;
    if NShard gt 1 then divs := [divs[i] : i in [1..#divs] | (i mod NShard) eq Shard]; end if;
    printf "TYPE %o: %o divisors (up to diamond%o), qlist=%o\n", lambda, #divs, NShard gt 1 select Sprintf(", shard %o of %o", Shard, NShard) else "", qlist;
    if #qlist eq 0 then printf "ERROR type=%o: no admissible q -- all %o divisors count as survivors\n", lambda, #divs; end if;
    surv := 0; cnt := 0;
    for Dv in divs do
      cnt +:= 1; killed := false; err := false;
      for q in qlist do
        try
          AqD := HeckeOp(X, q, Dv) - q*DiamondOperator(X, q, Dv) - Dv;
          if not IsPrincipal(AqD) then killed := true; break; end if;
        catch e
          err := true; printf "ERROR type=%o q=%o divisor=%o msg=%o\n", lambda, q, Describe(Dv), e`Object;
        end try;
      end for;
      if not killed then
        surv +:= 1; Append(~allsurv, <lambda, Dv>);
        try printf "SURVIVOR type=%o err=%o divisor=%o\n", lambda, err, Describe(Dv); catch e printf "SURVIVOR type=%o err=%o divisor=(describe failed)\n", lambda, err; end try;
        printf "SURVIVOR_DATA type=%o data=%o\n", lambda, DescribeData(Dv);
      end if;
      if cnt mod 100 eq 0 then printf "PROGRESS type=%o %o/%o survivors=%o time=%o\n", lambda, cnt, #divs, surv, Cputime(t0); end if;
    end for;
    printf "RESULT type=%o total=%o survivors=%o time=%o\n", lambda, #divs, surv, Cputime(t0);
  end for;
end for;
printf "SUMMARY N=%o M=%o p=%o D=%o shard=%o/%o total_survivors=%o time=%o\n", N, M, p, D, Shard, NShard, #allsurv, Cputime(t0);
print "SIEVE_DONE";
quit;
