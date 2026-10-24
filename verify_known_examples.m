/* Verify the other two known sporadic sextic groups.
   C25 point: van Hoeij's LowDegreePlaces, entry N=25, degv=6.
   C4 x C12 point: rational curve 162d1 over Q(alpha,i), alpha^3-3alpha-4=0. */
SetColumns(0);
SetMemoryLimit(3*10^9);
R<x> := PolynomialRing(Rationals());
f := x^3-x^2+1;
assert IsIrreducible(f);
L<a> := NumberField(f);
S<t> := PolynomialRing(L);
g := t^2+(a^2-2*a-1)*t+a;
assert IsIrreducible(g);
K<b> := ext<L|g>;
assert AbsoluteDegree(K) eq 6;
xx := K!a; yy := b;
rr := (xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));
ss := (xx*yy-yy+1)/(xx*yy);
bb := rr*ss*(rr-1); cc := ss*(rr-1);
E := EllipticCurve([1-cc,-bb,-bb,K!0,K!0]);
P := E![0,0,1];
assert 25*P eq E!0 and 5*P ne E!0;
counts := [];
for p in PrimesUpTo(200) do
 if p le 3 or p eq 5 then continue; end if;
 F := GF(p); RF<z> := PolynomialRing(F);
 found := false;
 for rt in Roots(RF!f) do
  aa := rt[1];
  for uu in Roots(z^2+(aa^2-2*aa-1)*z+aa) do
   yy := uu[1]; xx := aa;
   if xx*yy*(xx*yy-1) eq 0 then continue; end if;
   rr := (xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));
   ss := (xx*yy-yy+1)/(xx*yy);
   bb := rr*ss*(rr-1); cc := ss*(rr-1);
   Ep := EllipticCurve([1-cc,-bb,-bb,F!0,F!0]);
   if Discriminant(Ep) eq 0 then continue; end if;
   nn := #Ep; assert nn mod 25 eq 0;
   Append(~counts,<p,nn>); found := true;
   printf "C25_GOOD_REDUCTION p=%o order=%o\n",p,nn;
   break;
  end for;
  if found then break; end if;
 end for;
end for;
primes := Seqset([c[1] : c in counts]);
for c in counts do primes join:= Seqset(PrimeDivisors(c[2])); end for;
bound := 1;
for ell in primes do
 vals := [Valuation(c[2],ell) : c in counts | c[1] ne ell];
 assert #vals gt 0;
 bound *:= ell^Minimum(vals);
end for;
printf "C25_EXAMPLE field_degree=%o torsion_order_upper_bound=%o\n",AbsoluteDegree(K),bound;
assert bound eq 25;
L<a> := NumberField(x^3-3*x-4);
S<t> := PolynomialRing(L);
assert IsIrreducible(t^2+1);
K<i> := ext<L|t^2+1>;
E0 := EllipticCurve("162d1");
E := BaseChange(E0,K);
A,h := TorsionSubgroup(E);
printf "C4_C12_EXAMPLE rational_model=%o field_degree=%o torsion=%o\n",aInvariants(E0),AbsoluteDegree(K),Invariants(A);
assert AbsoluteDegree(K) eq 6 and Invariants(A) eq [4,12];
print "KNOWN_EXAMPLES_DONE";
