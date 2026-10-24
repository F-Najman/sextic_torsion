SetColumns(0);
SetMemoryLimit(10^9);
AttachSpec("gl2/magma.spec");
t0:=Cputime();m:=2;n:=24;p:=5;
G:=GL(2,Integers(n));
H:=sub<G|[G![1,0,m,1]] cat [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 and (a-1) mod m eq 0]>;
ycnt:=[&+GL2jCounts(H,p^d):d in [1..6]];
a:=[];
for d in [1..6] do
 z:=ycnt[d]-&+([0] cat [e*a[e]:e in [1..d-1]|d mod e eq 0]);
 assert z mod d eq 0;Append(~a,z div d);
end for;
assert GL2Genus(H) eq 13;
assert a eq [0,0,16,96,672,2712];
print "INDEPENDENT_M2N24_COUNTS",ycnt,a,Cputime(t0);
print "INDEPENDENT_M2N24_COUNTS_DONE complete=true";
quit;
