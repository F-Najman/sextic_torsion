SetColumns(0);
SetMemoryLimit(1500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
for nd in [<9,4,3,7>,<15,7,2,13>] do
 N,p,d,q:=Explode(nd);
 X:=MDX11(3,N,GF(p));
 xs:=NoncuspidalPlaces(X,d); assert #xs gt 0;
 x:=xs[1]; E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
 time fast:=HeckeOperatorFast(X,q,x);
 t0:=Cputime();
 isog:=MDIsogenies(E,q);ys:=[];ms:=[];
 for phi in isog do
  EG:=Codomain(phi);
  LsG:=rec<recformat<P:PtEll,Q:PtEll>|P:=phi(Ls`P),Q:=phi(Ls`Q)>;
  y:=ModuliPoint3(X,EG,LsG);
  ii:=Index(ys,y);
  if ii eq 0 then Append(~ys,y);Append(~ms,1); else ms[ii]+:=1; end if;
 end for;
 reference:=DivisorGroup(Curve(X))!0;
 for i in [1..#ys] do reference+:=(Integers()!(ms[i]*Degree(x)/Degree(ys[i])))*Divisor(ys[i]); end for;
 print "REFERENCE_TIME",Cputime(t0);
 assert fast eq reference;
 print "HECKE_FAST_REFERENCE_EQUAL",3,N,p,d,q,Degree(fast);
end for;
print "FAST_HECKE_3_TEST_DONE";
quit;
