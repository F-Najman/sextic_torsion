SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
t0:=Cputime();X:=MDX11(4,20,GF(9));
assert Genus(X) eq 25 and #NoncuspidalPlaces(X,1) eq 0;
for d in [2,3] do
 xs:=NoncuspidalPlaces(X,d);
 left:=Seqset(xs);rs:=[];
 while #left gt 0 do
  x:=Rep(left);Append(~rs,x);
  for a in [1,3,7,9] do
   y:=DiamondOperator(X,a,x);assert y in xs;Exclude(~left,y);
  end for;
 end while;
 print "M4N20_PREFLIGHT_ORBITS",d,#xs,#rs,Cputime(t0);
 assert <d,#xs,#rs> in [<2,96,24>,<3,320,80>];
 for x in xs[1..4] do
  E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
  assert ModuliPoint4(X,E,Ls) eq x;
 end for;
 if d eq 2 then x2:=rs[1];end if;
end for;
E:=EllipticCurve(X,x2);Ls:=LevelStructure(X,x2);q:=7;
fast:=HeckeOperatorQuotientFast(X,q,x2);ys:=[];ms:=[];ti:=Cputime();
for phi in MDIsogenies(E,q) do
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=InverseMod(q,4)*phi(Ls`P),Q:=phi(Ls`Q)>;
 y:=ModuliPoint4(X,Codomain(phi),La);ii:=Index(ys,y);
 if ii eq 0 then Append(~ys,y);Append(~ms,1);else ms[ii]+:=1;end if;
end for;
reference:=DivisorGroup(Curve(X))!0;
for i in [1..#ys] do reference+:=(Integers()!(ms[i]*2/Degree(ys[i])))*Divisor(ys[i]);end for;
assert fast eq reference;
print "M4N20_Q7_REFERENCE_EQUAL",Degree(fast),Cputime(ti);
ti:=Cputime();T:=HeckeOperatorFast(X,29,x2);
A:=T-Divisor(x2)-29*Divisor(DiamondOperator(X,29,x2));
print "M4N20_T29_BENCHMARK",Degree(T),IsPrincipal(A),Cputime(ti);
print "M4N20_PREFLIGHT_DONE",Cputime(t0);
quit;
