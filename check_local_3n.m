SetColumns(0);
for N in [15,21,27] do
 for q in [4,16,64] do
  F:=GF(q); structures:={};
  for j in F do
   for E in Twists(EllipticCurveFromjInvariant(j)) do
    inv:=Invariants(AbelianGroup(E));
    Include(~structures,inv);
    assert not (#inv eq 2 and inv[1] mod 3 eq 0 and inv[2] mod N eq 0);
   end for;
  end for;
  print "NO_C3_CN",N,q,#structures;
 end for;
end for;
print "LOCAL_3N_DONE";
quit;
