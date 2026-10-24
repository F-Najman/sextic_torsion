SetColumns(0);
SetMemoryLimit(1500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
for nd in [<6,7,2>,<9,7,2>,<12,13,2>,<15,7,2>,<9,4,3>] do
 N,p,d:=Explode(nd);
 X:=MDX11(3,N,GF(p));
 time xs:=NoncuspidalPlaces(X,d);
 print "PLACES",N,p,d,#xs;
 assert #xs gt 0;
 for x in xs[1..Min(#xs,12)] do
  E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
  y:=ModuliPoint3(X,E,Ls);
  assert y eq x;
  K:=BaseRing(E);
  for coeffs in [[PrimitiveElement(K),K!1,K!1,K!1],[K!1,K!0,K!0,K!0]] do
   u,r,s,t0:=Explode(coeffs);
   a1,a2,a3,a4,a6:=Explode(aInvariants(E));
   bs:=[(a1+2*s)/u,(a2-s*a1+3*r-s^2)/u^2,
        (a3+r*a1+2*t0)/u^3,
        (a4-s*a3+2*r*a2-(t0+r*s)*a1+3*r^2-2*s*t0)/u^4,
        (a6+r*a4+r^2*a2+r^3-t0*a3-r*t0*a1-t0^2)/u^6];
   E2:=EllipticCurve(bs);
   images:=[E2![(A[1]-r)/u^2,(A[2]-t0-s*(A[1]-r))/u^3] : A in [Ls`P,Ls`Q]];
   Ls2:=rec<recformat<P:PtEll,Q:PtEll>|P:=images[1],Q:=images[2]>;
   assert ModuliPoint3(X,E2,Ls2) eq x;
  end for;
 end for;
 print "MODULI_POINT3_ROUNDTRIP",N,p,d,Min(#xs,12);
end for;
print "MODULI_POINT3_TEST_DONE";
quit;
