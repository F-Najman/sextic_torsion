SetColumns(0);
function ExactOrderPoints(inv,n)
 return &+[MoebiusMu(n div d)*(&*[GCD(d,a):a in inv]):d in Divisors(n)];
end function;
for N in [18,24] do
 counts:=[];
 for d in [1..3] do
  F:=GF(7^d);total:=0;
  for j in F do
   for E in Twists(EllipticCurveFromjInvariant(j)) do
    inv:=Invariants(AbelianGroup(E));
    if #inv ne 2 or inv[1] mod 3 ne 0 then continue;end if;
    nn:=ExactOrderPoints(inv,N);
    aut:=#Automorphisms(E);
    summand:=3*nn/aut;
    assert Denominator(summand) eq 1;
    total+:=Integers()!summand;
   end for;
  end for;
  Append(~counts,total);
 end for;
 expected:=(N eq 18) select [0,54,432] else [0,0,144];
 assert counts eq expected;
 print "INDEPENDENT_M3_COUNTS",N,counts;
end for;
print "INDEPENDENT_M3_COUNTS_DONE";
quit;
