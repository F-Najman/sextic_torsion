/* Check the characteristic-two calculation against the reference correspondence. */
SetColumns(0);
SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
X := MDX1(39,GF(2));
places := NoncuspidalPlaces(X,6);
reps := PlacesUpToDiamondFast(X,places);
assert #places eq 12 and #reps eq 1;
for q in [5,7] do
 Dfast := HeckeOperatorFast(X,q,reps[1]);
 Dref := HeckeOperator(X,q,reps[1]);
 assert Dfast eq Dref;
 printf "CHAR2_HECKE q=%o exact_divisors_agree=true\n",q;
end for;
print "CHAR2_HECKE_DONE";
