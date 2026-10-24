/* Independent Frobenius-class point counts and direct Hecke comparison. */
SetColumns(0);
SetMemoryLimit(2*10^9);
AttachSpec("gl2/magma.spec");
n:=35;
G:=GL(2,Integers(n));
H:=sub<G | [G![1,0,1,1]] cat
 [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1]>;
assert GL2Level(H) eq n and #H eq n*EulerPhi(n);
assert GL2Genus(H) eq 25;
counts:=[&+GL2jCounts(H,3^d):d in [1..6]];
expected:=[0,0,36,0,180,792];
assert counts eq expected;
printf "GL2_COUNT_35 p=3 noncuspidal_rational_point_counts=%o agree=true\n",counts;

AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
X:=MDX1(35,GF(3));
xs:=NoncuspidalPlaces(X,3);
assert #xs eq 12;
for q in [11,13] do
 tf:=HeckeOperatorFast(X,q,xs[1]);
 tr:=HeckeOperator(X,q,xs[1]);
 assert tf eq tr;
 printf "HECKE_REFERENCE_35 p=3 place_degree=3 q=%o exact_divisors_agree=true\n",q;
end for;
print "VERIFY_35_DONE";
quit;
