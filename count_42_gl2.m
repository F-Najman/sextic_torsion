SetColumns(0);
SetMemoryLimit(10^9);
AttachSpec("gl2/magma.spec");
n:=42;p:=5;G:=GL(2,Integers(n));
H:=sub<G|[G![1,0,1,1]] cat [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1]>;
assert GL2Level(H) eq n and #H eq n*EulerPhi(n) and GL2Genus(H) eq 25;
ycnt:=[&+GL2jCounts(H,p^d):d in [1..6]];
a:=[];
for d in [1..6] do
 z:=ycnt[d]-&+([0] cat [e*a[e]:e in [1..d-1]|d mod e eq 0]);
 assert z mod d eq 0;Append(~a,z div d);
end for;
print "GL2_COUNT_42",ycnt,a;
print "GL2_COUNT_42_DONE";
quit;
