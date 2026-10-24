/* Identify the reductions of both exact Q-points with the two
   representatives in the exhaustive finite-class certificate. */
SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_diamond_hash_m2.m");
X:=MDX11(2,24,GF(5));C:=Curve(X);
xs:=NoncuspidalPlaces(X,6);rs:=PlacesUpToDiamondM2HashDirect(X,xs);
assert #xs eq 2712 and #rs eq 684;
R<V>:=PolynomialRing(Rationals());
fs:=[V^6-3*V^5+V^4/2+4*V^3-V^2-3*V/2+1/2,
 V^6+4*V^5+41*V^4/6+31*V^3/6+4*V^2/3-V/3-1/6];
hs:=[3-3*V-13*V^2+2*V^3+10*V^4-4*V^5,
 (-92-15*V+542*V^2+1074*V^3+744*V^4+216*V^5)/17];
indices:=[87,577];Fp:=GF(5);S<T>:=PolynomialRing(Fp);
for j in [1,2] do
 f:=S!fs[j];h:=S!hs[j];assert IsIrreducible(f);
 K<a>:=ext<Fp|f>;pt:=C(K)![Evaluate(h,a),a,K!1];
 ps:=Places(pt);assert #ps eq 1 and Degree(ps[1]) eq 6;
 assert ps[1] eq rs[indices[j]];
 print "N24_EXACT_POINT_SURVIVOR_MATCH",j,indices[j],"exact=true";
end for;
print "N24_EXACT_POINT_MATCH_DONE complete=true",Cputime(t0);quit;
