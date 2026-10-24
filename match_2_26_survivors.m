/* Match the reconstructed survivors in the exhaustive sieve to the exact
   representatives used by the independent rational-cusp class calculation. */
SetColumns(0);
SetMemoryLimit(2*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
X:=MDX11(2,26,GF(3));
xs:=NoncuspidalPlaces(X,6);
rs:=PlacesUpToDiamondFast(X,xs);
assert #xs eq 122 and #rs eq 26;
defs:=[[2,0,0,1,2,2,1],[2,0,2,0,1,1,1]];
us:=[[2,1,0,0,2,1],[2,0,2,0,2,1]];
indices:=[10,16];
R<t>:=PolynomialRing(GF(3));
for k in [1,2] do
 K<a>:=ext<GF(3)|R!defs[k]>;
 u:=Evaluate(R!us[k],a);
 pt:=Curve(X)(K)![u,a,1];
 places:=Places(pt);
 assert #places eq 1;
 assert places[1] eq rs[indices[k]];
 printf "M2N26_SURVIVOR_MATCH %o exact=true\n",indices[k];
end for;
print "M2N26_SURVIVOR_MATCH_DONE";
quit;
