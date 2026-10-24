SetColumns(0);
SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("hecke_two_commonzero.m");
t0:=Cputime();
for np in [<35,3>,<33,5>] do
 N,p:=Explode(np);X:=MDX1(N,GF(p));
 xs:=NoncuspidalPlaces(X,3);assert #xs gt 0;
 for x in xs[1..Min(2,#xs)] do
  actual:=HeckeTwoCommonZeros(X,x);expected:=HeckeOperator(X,2,x);
  assert actual eq expected;
  E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
  assert ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=2*P>) eq DiamondOperator(X,2,x);
  printf "HECKE_TWO_COMMONZERO_EQUAL N=%o p=%o degree=%o time=%o\n",N,p,Degree(x),Cputime(t0);
 end for;
end for;
print "HECKE_TWO_COMMONZERO_TEST_DONE";
quit;
