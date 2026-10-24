/* Complementary sign-representation filter for the completed three
   X1(30) quotient survivor export. No assertion that the whole target
   Jacobian overQ is cuspidal is needed: DEHMZ Prop4.16 makes every
   source rational Jacobian class cuspidal, so its pushforward lies in
   the target rational cusp-divisor subgroup used below. */
SetColumns(0);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11_m2_guarded.m");
Attach("fast_hecke_cyclic_commonzero.m");
load "a3_m2_30_model.m";
load "data/m2_30_three_quotient_survivors.m";
X:=MDX11(2,30,GF(7));CX:=Curve(X);assert Genus(X) eq 25;
assert X`_E[1] eq 0 and X`_E[3] eq 0;
printf "A3_M2_30_FILTER_INPUT count=%o time=%o\n",#M2_30_ThreeQuotientSurvivors,Cputime(t0);
A3G,_,A3phi:=ClassGroup(A3K:ReductionDivisor:=Divisor(A3c0),Proof:=true);
A3hs:=[A3phi(OO-Degree(OO)*Divisor(A3c0)):OO in A3orbits];
A3H:=sub<A3G|A3hs>;A3Q,A3qbar:=quo<A3G|A3H>;
A3cc:=[A3phi(Divisor(c)-Degree(c)*Divisor(A3c0)):c in A3cs];
A3ds:=[Degree(c):c in A3cs];
print "A3_M2_30_FILTER_GROUP",Invariants(A3G),Invariants(A3H),Invariants(A3Q),Cputime(t0);
procedure A3AddComplement(remaining,start,current,~out)
 if remaining eq 0 then Include(~out,current);return;end if;
 for i in [start..#A3cs] do
  if A3ds[i] le remaining then
   A3AddComplement(remaining-A3ds[i],i,current+A3qbar(A3cc[i]),~out);
  end if;
 end for;
end procedure;
A3complements:=AssociativeArray();
for k in [0..4] do
 out:={A3Q|};A3AddComplement(k,1,A3Q!0,~out);
 A3complements[k]:={-a:a in out};
 print "A3_M2_30_COMPLEMENTS",k,#out;
end for;
function A3ReconstructPart(data)
 D:=DivisorGroup(CX)!0;
 for row in data do
  dd,mul,def,coords:=Explode(row);L:=ext<GF(7)|Polynomial(GF(7),def)>;
  assert #L eq 7^dd;
  qs:=Places(CX(L)![L!cs:cs in coords]);
  assert #qs eq 1 and Degree(qs[1]) eq dd and not IsCusp(X,qs[1]);
  D+:=mul*Divisor(qs[1]);
 end for;
 return D;
end function;
A3cache:=AssociativeArray();A3remaining:=[* *];A3records:=[* *];A3number:=0;
for idx in [1..#M2_30_ThreeQuotientSurvivors] do
 row:=M2_30_ThreeQuotientSurvivors[idx];lambda,data:=Explode(row);
 D:=A3ReconstructPart(data);dd:=Degree(D);
 assert dd eq &+lambda and dd ge 2 and dd le 6;
 ss,mm:=Support(D);expanded:=&cat[[Degree(ss[j]):k in [1..mm[j]]]:j in [1..#ss]];
 assert Sort(expanded) eq Sort(lambda);
 aa:=A3G!0;parts:=[* *];
 for j in [1..#ss] do
  source:=ss[j];de:=Degree(source);E:=EllipticCurve(X,source);LL:=LevelStructure(X,source);
  values:=A3OrientedValues(E,LL`P,LL`Q);
  key:=<de,A3FrobeniusTupleKey(values,A3fields[de],de)>;
  if IsDefined(A3cache,key) then
   gc,targetdegree,residuemult:=Explode(A3cache[key]);
  else
   target,residuemult,checkkey:=A3PlaceFromLevel(E,LL`P,LL`Q);
   assert checkkey eq key[2];targetdegree:=Degree(target);
   assert targetdegree*residuemult eq de;
   gc:=residuemult*A3phi(Divisor(target)-targetdegree*Divisor(A3c0));
   A3cache[key]:=<gc,targetdegree,residuemult>;A3number+:=1;
  end if;
  aa+:=mm[j]*gc;
  Append(~parts,[*de,mm[j],targetdegree,residuemult,Eltseq(gc)*]);
 end for;
 passes:=A3qbar(aa) in A3complements[6-dd];
 Append(~A3records,[*idx,lambda,Eltseq(aa),passes,parts*]);
 print "A3_M2_30_FILTER_TEST",idx,lambda,passes,Eltseq(aa),Cputime(t0);
 if passes then Append(~A3remaining,row);end if;
end for;
Write("data/m2_30_a3_survivors.m",
 "M2_30_A3Survivors := " cat Sprint(A3remaining) cat ";" : Overwrite:=true);
A3arithmetic:=[*Invariants(A3G),[Eltseq(r):r in Rows(RelationMatrix(A3G))],
 [Eltseq(a):a in A3hs],A3ds,[Eltseq(a):a in A3cc],A3records*];
Write("data/A3_m2_30_arithmetic.txt",Sprint(A3arithmetic) : Overwrite:=true);
print "A3_M2_30_ARITHMETIC_BEGIN";print A3arithmetic;print "A3_M2_30_ARITHMETIC_END";
printf "A3_M2_30_FILTER_DONE input=%o remaining=%o class_maps=%o complete=true time=%o\n",
 #M2_30_ThreeQuotientSurvivors,#A3remaining,A3number,Cputime(t0);
quit;
