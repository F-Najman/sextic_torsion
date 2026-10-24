/* All j-invariants and all twists over F64. Hasse handles smaller fields. */
SetColumns(0);F:=GF(64);orders:={};curves:=0;
for j in F do
 E:=EllipticCurveFromjInvariant(j);
 for A in [E] cat Twists(E) do
  n:=#A;assert n mod 63 ne 0;Include(~orders,n);curves+:=1;
 end for;
end for;
print "C63_F64_ALL_TWISTS",curves,Sort(Setseq(orders));
print "C63_F64_DONE no_order_divisible_by63=true";
quit;
