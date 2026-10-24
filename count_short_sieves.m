SetColumns(0);
SetMemoryLimit(10^9);
AttachSpec("gl2/magma.spec");
for t in [<2,26,3>,<2,28,3>,<1,36,5>,<1,40,3>,<1,54,5>] do
 m,n,p:=Explode(t);
 G:=GL(2,Integers(n));
 H:=sub<G | [G![1,0,m,1]] cat
  [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 and (a-1) mod m eq 0]>;
 ycnt:=[&+GL2jCounts(H,p^d):d in [1..6]];
 a:=[];
 for d in [1..6] do
  z:=ycnt[d]-&+([0] cat [e*a[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;
  Append(~a,z div d);
 end for;
 printf "SHORT_SIEVE_COUNT m=%o n=%o p=%o genus=%o rational_point_counts=%o closed_place_counts=%o\n",
  m,n,p,GL2Genus(H),ycnt,a;
end for;
print "SHORT_SIEVE_COUNT_DONE";
quit;
