SetColumns(0);
SetMemoryLimit(2500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
t0:=Cputime();X:=MDX11(4,16,GF(9));
for d in [1..3] do
 xs:=NoncuspidalPlaces(X,d);
 print "M4N16_F9_PREFLIGHT",Genus(X),d,#xs,Cputime(t0);
end for;
left:=Seqset(xs);rs:=[];
while #left gt 0 do
 x:=Rep(left);Append(~rs,x);
 for a in [1,3,5,7] do
  y:=DiamondOperator(X,a,x);assert y in xs;Exclude(~left,y);
 end for;
end while;
print "M4N16_F9_DIAMOND_ORBITS",#xs,#rs,Cputime(t0);
x:=rs[1];E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
fast:=HeckeOperatorFast(X,5,x);ys:=[];ms:=[];
for phi in MDIsogenies(E,5) do
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
 y:=ModuliPoint4(X,Codomain(phi),La);ii:=Index(ys,y);
 if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
end for;
reference:=DivisorGroup(Curve(X))!0;
for i in [1..#ys] do reference+:=(Integers()!(ms[i]*3/Degree(ys[i])))*Divisor(ys[i]);end for;
assert fast eq reference;
print "M4N16_F9_HECKE_REFERENCE_EQUAL",Degree(fast),Cputime(t0);
print "M4N16_F9_PREFLIGHT_DONE";
quit;
