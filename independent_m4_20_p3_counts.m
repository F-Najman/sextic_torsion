SetColumns(0);
function ExactOrderPoints(inv,n)
 return &+[MoebiusMu(n div d)*(&*[GCD(d,a):a in inv]):d in Divisors(n)];
end function;
counts:=[];
for d in [1..3] do
 F:=GF(9^d);total:=0;eligible:=0;
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   inv:=Invariants(AbelianGroup(E));
   if #inv ne 2 or inv[1] mod 4 ne 0 then continue;end if;
   nn:=ExactOrderPoints(inv,20);if nn eq 0 then continue;end if;
   summand:=4*nn/#Automorphisms(E);
   assert Denominator(summand) eq 1;
   total+:=Integers()!summand;eligible+:=1;
  end for;
 end for;
 Append(~counts,total);
 print "INDEPENDENT_M4N20_P3_COUNT",d,total,eligible;
end for;
assert counts eq [0,192,960];
print "INDEPENDENT_M4N20_P3_COUNTS_DONE",counts;
quit;
