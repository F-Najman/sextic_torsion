SetColumns(0);
SetMemoryLimit(1000000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
X:=MDX11(2,22,GF(3));
for d in [5,6] do
 time xs:=NoncuspidalPlaces(X,d);
 time reps:=PlacesUpToDiamondFast(X,xs);
 print "FAST_DIAMOND_COUNTS",d,#xs,#reps;
 assert <d,#xs,#reps> in [<5,60,12>,<6,100,20>];
end for;
print "FAST_DIAMOND_22_TEST_DONE";
quit;
