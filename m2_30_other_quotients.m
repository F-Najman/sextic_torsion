/* Necessary second-stage filters using the other two Q-defined maps
   X1(2,30) -> X1(30). Input must be the completed first-stage export.
   Q=15Q+16Q splits its 2- and 15-primary parts. The three selected
   nonzero 2-points are 15Q, P2, and P2+15Q. */
SetColumns(0);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("fast_diamond_hash_m2.m");
load "data/m2_30_quotient_survivors.m";
X:=MDX11(2,30,GF(7));CX:=Curve(X);assert Genus(X) eq 25;
assert X`_E[1] eq 0 and X`_E[3] eq 0;
Y:=MDX1(30,GF(7));CY:=Curve(Y);assert Genus(Y) eq 9;
osY:=CuspOrbitsQ(Y);c0Y:=Support([o:o in osY|Degree(o) eq 1][1])[1];
GY,m1Y,m2Y:=ClassGroupRationalReduction(CY,c0Y:Proof:=true);
assert Invariants(GY) eq [4,4,4,136,16320,0];
HQY:=sub<GY|[m2Y(o-Degree(o)*Divisor(c0Y)):o in osY]>;
assert Invariants(HQY) eq [4,8160];QY,qY:=quo<GY|HQY>;
csY:=Cusps(Y);dsY:=[Degree(x):x in csY];
ccY:=[qY(m2Y(Divisor(x)-Degree(x)*Divisor(c0Y))):x in csY];
procedure AddCuspCoset(remaining,start,current,~out)
 if remaining eq 0 then Include(~out,current);return;end if;
 for i in [start..#dsY] do
  if dsY[i] le remaining then AddCuspCoset(remaining-dsY[i],i,current+ccY[i],~out);end if;
 end for;
end procedure;
complementsY:=AssociativeArray();
for k in [0..4] do
 out:={QY|};AddCuspCoset(k,1,QY!0,~out);complementsY[k]:={-z:z in out};
 assert #out eq [1,9,41,67,68][k+1];
end for;
function ReconstructPart(data)
 D:=DivisorGroup(CX)!0;
 for row in data do
  dd,mul,def,coords:=Explode(row);K:=ext<GF(7)|Polynomial(GF(7),def)>;
  assert #K eq 7^dd;
  pp:=Places(CX(K)![K!cs:cs in coords]);
  assert #pp eq 1 and Degree(pp[1]) eq dd and not IsCusp(X,pp[1]);
  D+:=mul*Divisor(pp[1]);
 end for;
 return D;
end function;
fields:=[GF(7^d):d in [1..6]];cache:=AssociativeArray();number:=0;
procedure PushOtherClass(x,variant,~cache,~gc,~number)
 assert variant in [2,3];dd:=Degree(x);
 E:=EllipticCurve(X,x);L:=LevelStructure(X,x);P2:=L`P;Q:=L`Q;
 assert aInvariants(E)[1] eq 0 and aInvariants(E)[3] eq 0;
 assert 2*P2 eq E!0 and P2 ne E!0 and 30*Q eq E!0;
 assert P2 ne 15*Q;
 P:=variant eq 2 select P2+16*Q else P2+Q;
 assert 30*P eq E!0 and 15*P ne E!0 and 10*P ne E!0 and 6*P ne E!0;
 xy:=TateNormalForm_xy(E,P);key:=<dd,MDM2FrobeniusKey(xy,fields[dd],7,dd)>;
 if IsDefined(cache,key) then gc:=cache[key];return;end if;
 y:=ModuliPointCyclicCommonZeros(Y,E,rec<recformat<P:PtEll>|P:=P>);
 dy:=Degree(y);assert dd mod dy eq 0 and not IsCusp(Y,y);
 gc:=(dd div dy)*m2Y(Divisor(y)-dy*Divisor(c0Y));
 cache[key]:=gc;number+:=1;
end procedure;
remaining:=[* *];passed2:=0;passed3:=0;
for idx in [1..#M2_30_QuotientSurvivors] do
 row:=M2_30_QuotientSurvivors[idx];lambda,data:=Explode(row);D:=ReconstructPart(data);
 dd:=Degree(D);assert dd eq &+lambda and dd ge 2 and dd le 6;
 ss,mm:=Support(D);expanded:=&cat[[Degree(ss[j]):k in [1..mm[j]]]:j in [1..#ss]];
 assert Sort(expanded) eq Sort(lambda);
 killed:=false;
 for variant in [2,3] do
  a:=GY!0;
  for j in [1..#ss] do
   gc:=GY!0;PushOtherClass(ss[j],variant,~cache,~gc,~number);a+:=mm[j]*gc;
  end for;
  passes:=qY(a) in complementsY[6-dd];
  print "M2_30_OTHER_QUOTIENT_TEST",idx,lambda,variant,passes,Cputime(t0);
  if not passes then killed:=true;break;end if;
  if variant eq 2 then passed2+:=1;else passed3+:=1;end if;
 end for;
 if not killed then
  Append(~remaining,row);print "M2_30_THREE_QUOTIENT_SURVIVOR",idx,lambda,Cputime(t0);
 end if;
end for;
assert #remaining eq passed3;
Write("data/m2_30_three_quotient_survivors.m",
 "M2_30_ThreeQuotientSurvivors := " cat Sprint(remaining) cat ";" : Overwrite:=true);
printf "M2_30_OTHER_QUOTIENTS_DONE input=%o pass2=%o pass3=%o remaining=%o class_maps=%o complete=true time=%o\n",#M2_30_QuotientSurvivors,passed2,passed3,#remaining,number,Cputime(t0);
quit;
