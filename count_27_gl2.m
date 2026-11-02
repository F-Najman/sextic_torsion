/* Independent count of the non-cuspidal closed places of X_1(27) over F_p of degree at most 6,
   by counting elliptic curves with level structure (Sutherland's GL2 package), without a model
   of the curve.  Used for p = 5 in the proof that Z/54Z is not a sextic torsion group: compare
   the line C54_MOD5_PLACES of logs/classes/X1_27_C54_p5.log.
   Usage: magma -n -b [p:=5] count_27_gl2.m */
SetColumns(0);SetMemoryLimit(10^9);
if not assigned p then p:=5; elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
AttachSpec("gl2/magma.spec");
n:=27;G:=GL(2,Integers(n));
H:=sub<G|[G![1,0,1,1]] cat [G![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1]>;
assert GL2Level(H) eq n and #H eq n*EulerPhi(n) and GL2Genus(H) eq 13;
ycnt:=[&+GL2jCounts(H,p^d):d in [1..6]];a:=[];
for d in [1..6] do
 z:=ycnt[d]-&+([0] cat [e*a[e]:e in [1..d-1]|d mod e eq 0]);
 assert z mod d eq 0;Append(~a,z div d);
end for;
printf "GL2_COUNT_27 p=%o rational_point_counts=%o closed_place_counts=%o\n",p,ycnt,a;
print "GL2_COUNT_27_DONE";quit;
