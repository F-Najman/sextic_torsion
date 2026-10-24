SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
N:=26;t0:=Cputime();X:=MDX11(2,N,GF(3));C:=Curve(X);
function SigFunctions(Y)
 FF:=FunctionField(Curve(Y));u:=FF.1;v:=FF.2;
 b,c:=Explode([FF!f:f in Y`_coordinates]);
 E:=EllipticCurve([0,c,0,(1-b-c)*b,0]);Q:=E![b,b];
 fs:=[u,v,b,c,1-b-c,4*b^2+4*b*c-4*b+c^2];
 for k in [2..(N div 2)] do
  kQ:=k*Q;if kQ[3] ne 0 then Append(~fs,kQ[1]);end if;
 end for;
 return fs;
end function;
cs:=Cusps(X);fs:=SigFunctions(X);groups:=AssociativeArray();
for c in cs do
 sg:=[Valuation(f,c):f in fs];
 if IsDefined(groups,sg) then Append(~groups[sg],c);else groups[sg]:=[c];end if;
end for;
orbits:=[&+[Divisor(c):c in groups[sg]]:sg in Keys(groups)];
print "M2N26_CUSP_SIGNATURE_ORBITS",#groups,Sort([Degree(o):o in orbits]),Cputime(t0);
assert Sort([Degree(o):o in orbits]) eq [1:i in [1..18]] cat [6,6,6];
c0:=Support([o:o in orbits|Degree(o) eq 1][1])[1];
print "M2N26_CLASSGROUP_START",Cputime(t0);
G,m1,m2:=ClassGroup(X);
H:=sub<G|[m2(o-Degree(o)*Divisor(c0)):o in orbits]>;
H1:=sub<G|[m2(Divisor(c)-Divisor(c0)):c in cs|Degree(c) eq 1]>;
Hall:=sub<G|[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs]>;
print "M2N26_CUSPIDAL_GROUPS",Invariants(G),Invariants(H),Invariants(Hall),Cputime(t0);
print "M2N26_RATIONAL_CUSPS_ONLY",Invariants(H1),#H1,#H,Cputime(t0);
xs:=NoncuspidalPlaces(X,6);rs:=PlacesUpToDiamondFast(X,xs);
assert #xs eq 122 and #rs eq 26;
for i in [10,16] do
 x:=rs[i];cl:=m2(Divisor(x)-6*Divisor(c0));
 assert jInvariant(X,x) eq 0;
 print "M2N26_SURVIVOR_CLASS",i,cl in H,cl in Hall,Order(cl),Cputime(t0);
 print "M2N26_SURVIVOR_RATIONAL_CUSPS_ONLY",i,cl in H1;
end for;
print "M2N26_CLASS_INSPECT_DONE",Cputime(t0);
quit;
