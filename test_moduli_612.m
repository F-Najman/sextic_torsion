SetColumns(0);
SetMemoryLimit(2500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
Attach("model_6_12.m");
t0:=Cputime();X:=MDX612(GF(7));
print "M6_MODEL_SETUP",Genus(X),Cputime(t0);
for d in [1,2,3] do
 xs:=NoncuspidalPlaces(X,d);
 print "M6_MODEL_PLACE_COUNTS",d,#xs,Cputime(t0);
 for x in xs[1..Min(12,#xs)] do
  Ls:=LevelStructure612(X,x);E:=Curve(Ls`P);
  assert ModuliPoint612(X,E,Ls) eq x;
  K:=BaseRing(E);
  u:=PrimitiveElement(K);r:=K!1;s:=K!1;tt:=K!1;
  a1,a2,a3,a4,a6:=Explode(aInvariants(E));
  bs:=[(a1+2*s)/u,(a2-s*a1+3*r-s^2)/u^2,
       (a3+r*a1+2*tt)/u^3,
       (a4-s*a3+2*r*a2-(tt+r*s)*a1+3*r^2-2*s*tt)/u^4,
       (a6+r*a4+r^2*a2+r^3-tt*a3-r*tt*a1-tt^2)/u^6];
  E2:=EllipticCurve(bs);
  images:=[E2![(A[1]-r)/u^2,(A[2]-tt-s*(A[1]-r))/u^3]:A in [Ls`P,Ls`Q]];
  Ls2:=rec<recformat<P:PtEll,Q:PtEll>|P:=images[1],Q:=images[2]>;
  assert ModuliPoint612(X,E2,Ls2) eq x;
  print "M6_MODEL_INVERSE_EQUAL",d,Cputime(t0);
 end for;
end for;
rs:=Representatives612(X,xs);
print "M6_MODEL_DIAMOND_ORBITS",#xs,#rs,Cputime(t0);
print "M6_MODEL_TEST_DONE",Cputime(t0);
quit;
