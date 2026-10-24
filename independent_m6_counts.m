SetColumns(0);
function ExactOrderPoints(inv,n)
 return &+[MoebiusMu(n div d)*(&*[GCD(d,a):a in inv]):d in Divisors(n)];
end function;
counts:=[];
for d in [1..3] do
 F:=GF(7^d);total:=0;eligible:=0;
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   inv:=Invariants(AbelianGroup(E));
   if #inv ne 2 or inv[1] mod 6 ne 0 then continue;end if;
   nn:=ExactOrderPoints(inv,12);if nn eq 0 then continue;end if;
   // For each Q of exact order12 there are exactly6 choices of P
   // with e6(P,2Q) equal to a specified primitive sixth root.
   aut:=#Automorphisms(E);summand:=6*nn/aut;
   assert Denominator(summand) eq 1;
   total+:=Integers()!summand;eligible+:=1;
  end for;
 end for;
 Append(~counts,total);
 print "INDEPENDENT_M6_COUNT",d,total,eligible;
end for;
assert counts eq [0,0,288];
print "INDEPENDENT_M6_COUNTS_DONE",counts;
AttachSpec("mdmagma/v2/mdmagma.spec");
X:=MDX11(3,12,GF(7));assert Genus(X) eq 3;
print "M3N12_GENUS_CHECK",Genus(X);
quit;
