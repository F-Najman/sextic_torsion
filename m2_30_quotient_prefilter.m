/* Necessary genus-nine quotient obstruction for degree-six divisors on
   X1(2,30)/F7. All cusp complements on Y=X1(30) are allowed, giving a
   deliberate superset of the possible pushed-forward complements.
   No class is discarded on grounds that Y has gonality six. */
SetColumns(0);SetMemoryLimit(16*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
Attach("fast_diamond_hash_m2.m");
X:=MDX11(2,30,GF(7));CX:=Curve(X);assert Genus(X) eq 25;
assert X`_E[1] eq 0 and X`_E[3] eq 0;
places:=[];reps:=[];
for d in [1..6] do
 xx:=NoncuspidalPlaces(X,d);assert #xx eq [0,24,192,732,2736,20000][d];
 Append(~places,xx);Append(~reps,PlacesUpToDiamondM2HashDirect(X,xx));
 assert #reps[d] eq [0,6,48,186,684,5006][d];
 print "M2_30_QUOTIENT_PLACES",d,#xx,#reps[d],Cputime(t0);
end for;
Y:=MDX1(30,GF(7));CY:=Curve(Y);assert Genus(Y) eq 9;
osY:=CuspOrbitsQ(Y);c0Y:=Support([o:o in osY|Degree(o) eq 1][1])[1];
GY,m1Y,m2Y:=ClassGroupRationalReduction(CY,c0Y:Proof:=true);
assert Invariants(GY) eq [4,4,4,136,16320,0];
HQY:=sub<GY|[m2Y(o-Degree(o)*Divisor(c0Y)):o in osY]>;
assert Invariants(HQY) eq [4,8160];QY,qY:=quo<GY|HQY>;
csY:=Cusps(Y);dsY:=[Degree(x):x in csY];
ccY:=[qY(m2Y(Divisor(x)-Degree(x)*Divisor(c0Y))):x in csY];
print "M2_30_QUOTIENT_GROUP",Invariants(GY),Invariants(HQY),Invariants(QY),Cputime(t0);
procedure AddCuspCoset(remaining,start,current,~out)
 if remaining eq 0 then Include(~out,current);return;end if;
 for i in [start..#dsY] do
  if dsY[i] le remaining then AddCuspCoset(remaining-dsY[i],i,current+ccY[i],~out);end if;
 end for;
end procedure;
complementsY:=AssociativeArray();
for k in [0..4] do
 out:={QY|};AddCuspCoset(k,1,QY!0,~out);
 complementsY[k]:={-z:z in out};
 print "M2_30_QUOTIENT_COMPLEMENTS",k,#out;
end for;
fields:=[GF(7^d):d in [1..6]];classcache:=AssociativeArray();
procedure PushClass(x,~cache,~gc,~number)
 dd:=Degree(x);E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`Q;
 assert aInvariants(E)[1] eq 0 and aInvariants(E)[3] eq 0;
 xy:=TateNormalForm_xy(E,P);key:=<dd,MDM2FrobeniusKey(xy,fields[dd],7,dd)>;
 if IsDefined(cache,key) then gc:=cache[key];return;end if;
 y:=ModuliPointCyclicCommonZeros(Y,E,rec<recformat<P:PtEll>|P:=P>);
 dy:=Degree(y);assert dd mod dy eq 0 and not IsCusp(Y,y);
 gc:=(dd div dy)*m2Y(Divisor(y)-dy*Divisor(c0Y));
 cache[key]:=gc;number+:=1;
end procedure;
smallY:=AssociativeArray();number:=0;
for d in [2..4] do
 for i in [1..#places[d]] do
  x:=places[d][i];a:=GY!0;PushClass(x,~classcache,~a,~number);smallY[x]:=a;
 end for;
 print "M2_30_QUOTIENT_SMALL",d,#places[d],number,Cputime(t0);
end for;
function DescribePart(D)
 sup,mul:=Support(D);data:=[* *];
 for i in [1..#sup] do
  dd:=Degree(sup[i]);pt:=RepresentativePoint(sup[i]);K:=GF(7^dd);
  coords:=MDM2ToCommonField(Eltseq(pt),K);
  def:=[Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(7)))];
  assert #def eq dd+1;
  coeffs:=[[Integers()!a:a in Eltseq(z,GF(7))]:z in coords];
  Append(~data,[*dd,mul[i],def,coeffs*]);
 end for;
 return data;
end function;
survivors:=[* *];
for d in [2..6] do
 passed:=0;
 for i in [1..#reps[d]] do
  x:=reps[d][i];a:=GY!0;
  if d le 4 then a:=smallY[x];else PushClass(x,~classcache,~a,~number);end if;
  if qY(a) in complementsY[6-d] then
   passed+:=1;Append(~survivors,[*[d],DescribePart(Divisor(x))*]);
   print "M2_30_QUOTIENT_SURVIVOR",[d],i,Eltseq(a),Cputime(t0);
  end if;
  if i mod 100 eq 0 then print "M2_30_QUOTIENT_PROGRESS",[d],i,#reps[d],passed,number,Cputime(t0);end if;
 end for;
 print "M2_30_QUOTIENT_TYPE",[d],#reps[d],passed,Cputime(t0);
end for;
for ds in [[2,2],[3,2],[4,2],[3,3]] do
 da,db:=Explode(ds);passed:=0;total:=0;
 for x in reps[da] do
  for y in places[db] do
   a:=smallY[x]+smallY[y];total+:=1;
   if qY(a) in complementsY[6-da-db] then
    passed+:=1;Append(~survivors,[*ds,DescribePart(Divisor(x)+Divisor(y))*]);
    print "M2_30_QUOTIENT_SURVIVOR",ds,total,Eltseq(a),Cputime(t0);
   end if;
  end for;
 end for;
 assert total eq #reps[da]*#places[db];
 print "M2_30_QUOTIENT_TYPE",ds,total,passed,Cputime(t0);
end for;
passed:=0;total:=0;
for x in reps[2] do
 for i in [1..#places[2]] do
  for j in [i..#places[2]] do
   y:=places[2][i];z:=places[2][j];a:=smallY[x]+smallY[y]+smallY[z];total+:=1;
   if qY(a) eq QY!0 then
    passed+:=1;Append(~survivors,[*[2,2,2],DescribePart(Divisor(x)+Divisor(y)+Divisor(z))*]);
    print "M2_30_QUOTIENT_SURVIVOR",[2,2,2],total,Eltseq(a),Cputime(t0);
   end if;
  end for;
 end for;
end for;
assert total eq #reps[2]*Binomial(#places[2]+1,2);
print "M2_30_QUOTIENT_TYPE",[2,2,2],total,passed,Cputime(t0);
Write("data/m2_30_quotient_survivors.m","M2_30_QuotientSurvivors := " cat Sprint(survivors) cat ";" : Overwrite:=true);
print "M2_30_QUOTIENT_DATA_BEGIN";print survivors;print "M2_30_QUOTIENT_DATA_END";
printf "M2_30_QUOTIENT_DONE survivors=%o class_maps=%o complete=true time=%o\n",#survivors,number,Cputime(t0);
quit;
