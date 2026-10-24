SetColumns(0);
SetMemoryLimit(2500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
Attach("model_6_12.m");
t0:=Cputime();X:=MDX612(GF(7));x:=NoncuspidalPlaces(X,3)[1];
E:=EllipticCurve(X,x);Ls:=LevelStructure612(X,x);
for q in [5,13] do
 ti:=Cputime();fast:=HeckeOperator612(X,q,x);
 print "M6_HECKE_FAST",q,Degree(fast),Cputime(ti);
 ys:=[];ms:=[];ti:=Cputime();
 for phi in MDIsogenies(E,q) do
  La:=rec<recformat<P:PtEll,Q:PtEll>|P:=InverseMod(q,6)*phi(Ls`P),Q:=phi(Ls`Q)>;
  y:=ModuliPoint612(X,Codomain(phi),La);ii:=Index(ys,y);
  if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
 end for;
 reference:=DivisorGroup(Curve(X))!0;
 for i in [1..#ys] do reference+:=(Integers()!(ms[i]*3/Degree(ys[i])))*Divisor(ys[i]);end for;
 assert fast eq reference;
 print "M6_HECKE_REFERENCE_EQUAL",q,Cputime(ti);
end for;
print "M6_HECKE_TEST_DONE",Cputime(t0);
quit;
