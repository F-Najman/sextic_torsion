/* Finite-field obstruction for C3 x C27 over a cubic extension of Q(zeta_3).
   Over primes above 2 the residue sizes can only be 4,16,64. */
SetColumns(0);
SetMemoryLimit(10^9);
for q in [4,16,64] do
 F := GF(q); invs := {};
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   inv := Invariants(AbelianGroup(E));
   if #inv eq 0 then inv := [1]; end if;
   Include(~invs,inv);
   assert not (#inv eq 2 and inv[1] mod 3 eq 0 and inv[2] mod 27 eq 0);
  end for;
 end for;
 printf "NO_C3_C27 q=%o checked_group_structures=%o\n",q,#invs;
end for;
print "NONCYCLIC_FINITE_FIELDS_DONE";
