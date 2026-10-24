SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_hecke_cyclic_commonzero.m");
load "a3_m2_30_model.m";
X:=MDX11(2,30,GF(7));assert X`_E[1] eq 0 and X`_E[3] eq 0;
for dd in [2,3] do
 xs:=NoncuspidalPlaces(X,dd);
 for idx in [1..Minimum(4,#xs)] do
  x:=xs[idx];E:=EllipticCurve(X,x);LL:=LevelStructure(X,x);
  q,mult,key:=A3PlaceFromLevel(E,LL`P,LL`Q);
  assert Degree(q)*mult eq dd;
  // The A3 even cycle of the three nonzero2-points fixes the oriented
  // point, while a transposition sends it to the deck conjugate.
  R2:=15*LL`Q;S15:=16*LL`Q;
  qeven,meven,keven:=A3PlaceFromLevel(E,R2,LL`P+R2+S15);
  assert qeven eq q and meven eq mult and keven eq key;
  qodd,modd,kodd:=A3PlaceFromLevel(E,R2,LL`P+S15);
  assert modd*Degree(qodd) eq dd;
  vv:=A3OrientedValues(E,LL`P,LL`Q);vv[3]:=-vv[3];
  assert A3FrobeniusTupleKey(vv,A3fields[dd],dd) eq kodd;
  print "A3_M2_30_MAP_TEST",dd,idx,Degree(q),mult,"even_cycle=true odd_deck=true",Cputime(t0);
 end for;
end for;
print "A3_M2_30_MAP_TEST_DONE exact=true",Cputime(t0);
quit;
