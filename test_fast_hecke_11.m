SetColumns(0);
SetMemoryLimit(1500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
X := MDX11(2,10,GF(13));
xs := NoncuspidalPlaces(X,1);
print "POINTS",#xs;
assert #xs gt 0;
x:=xs[1];
for q in [3,7,11] do
 time fast := HeckeOperatorFast(X,q,x);
 time reference := HeckeOperator(X,q,x);
 assert fast eq reference;
 print "HECKE_FAST_REFERENCE_EQUAL",2,10,13,1,q,Degree(fast);
end for;
time fastreps := PlacesUpToDiamondFast(X,xs);
time reference := PlacesUpToDiamond(X,xs);
assert #fastreps eq #reference;
for x in xs do
 assert #[y : y in fastreps | y in DiamondOrbit(X,x)] eq 1;
end for;
print "DIAMOND_FAST_REFERENCE_EQUAL",#xs,#fastreps;
print "FAST_HECKE_11_TEST_DONE";
quit;
