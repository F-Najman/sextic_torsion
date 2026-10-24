/* Finite-class sieve for degree-six points on rank-zero X1(36),X1(42).
   Published inputs: J(Q)=Q-rational cuspidal divisor classes; gon_Q>6.
   Odd p gives injective reduction on rational torsion.  All-cuspidal
   reductions are excluded separately by Derickx--Najman Prop5.2(a).

   Every possible noncuspidal part D' and every cuspidal complement C
   is checked by exact finite-group arithmetic, without diamond orbits.
   A class in G_Q represented by an effective Q-cuspidal degree-six
   divisor is already impossible for a noncuspidal rational divisor,
   by the gonality hypothesis.  Remaining classes are printed as
   integer combinations of Q-cusp orbits for characteristic-zero RR.
   Only a completed run with zero remaining classes proves exclusion. */
SetColumns(0);
if not assigned N then N:=36; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned p then p:=5; elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned MemGB then MemGB:=4; elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
assert N in [36,42] and p eq 5;
SetMemoryLimit(MemGB*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
t0:=Cputime();
X:=MDX1(N,GF(p));C:=Curve(X);
cusps:=Cusps(X);orbits:=CuspOrbitsQ(X);
c0:=Support([O:O in orbits|Degree(O) eq 1][1])[1];
printf "CLASS_CUSPIDAL_SETUP N=%o p=%o genus=%o cusp_degrees=%o Qorbit_degrees=%o time=%o\n",
 N,p,Genus(X),Sort([Degree(c):c in cusps]),[Degree(O):O in orbits],Cputime(t0);
Cl,m1,m2:=ClassGroup(C);
printf "CLASS_CUSPIDAL_CLASSGROUP invariants=%o time=%o\n",Invariants(Cl),Cputime(t0);
assert #[e:e in Invariants(Cl)|e eq 0] eq 1;
A:=JOne(N);f:=FrobeniusPolynomial(A,p);
assert Degree(f) eq 2*Genus(X);
jorder:=Integers()!Evaluate(f,1);
assert &*[e:e in Invariants(Cl)|e ne 0] eq jorder;
cc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cusps];
qc:=[m2(O-Degree(O)*Divisor(c0)):O in orbits];
Hall:=sub<Cl|cc>;
Free:=FreeAbelianGroup(#orbits);
hQ:=hom<Free->Cl|qc>;
HQ:=Image(hQ);
assert HQ subset Hall;
expected:=N eq 36 select [12,252,7812] else [182,1092,131040];
assert Invariants(HQ) eq expected;
printf "CLASS_CUSPIDAL_GROUPS N=%o p=%o Jorder=%o Hall=%o Hallorder=%o HQ=%o HQorder=%o quotient_index=%o time=%o\n",
 N,p,jorder,Invariants(Hall),#Hall,Invariants(HQ),#HQ,jorder div #Hall,Cputime(t0);
Qbar,qbar:=quo<Cl|HQ>;

/* Enumerate all effective divisors of fixed degree on a specified set
   of closed places, retaining only their centered divisor classes. */
procedure AddEffectiveClasses(weights,classes,start,remaining,current,~out,~number)
 if remaining eq 0 then Include(~out,current);number+:=1;return;end if;
 for i in [start..#weights] do
  if weights[i] le remaining then
   AddEffectiveClasses(weights,classes,i,remaining-weights[i],current+classes[i],~out,~number);
  end if;
 end for;
end procedure;
complements:=AssociativeArray();
cdegs:=[Degree(c):c in cusps];
for k in [0..3] do
 cs:={Cl|};number:=0;
 AddEffectiveClasses(cdegs,cc,1,k,Cl!0,~cs,~number);
 bycoset:=AssociativeArray();
 for a in cs do
  z:=qbar(a);
  if IsDefined(bycoset,z) then Include(~bycoset[z],a);
  else bycoset[z]:={a};end if;
 end for;
 complements[k]:=bycoset;
 printf "CLASS_CUSPIDAL_COMPLEMENTS degree=%o effective_divisors=%o classes=%o cosets=%o time=%o\n",
  k,number,#cs,#Keys(bycoset),Cputime(t0);
end for;
qcuspidal:={Cl|};qnumber:=0;
AddEffectiveClasses([Degree(O):O in orbits],qc,1,6,Cl!0,~qcuspidal,~qnumber);
assert &and[a in HQ:a in qcuspidal];
printf "CLASS_CUSPIDAL_Q_EFFECTIVE degree=6 divisors=%o classes=%o time=%o\n",qnumber,#qcuspidal,Cputime(t0);

/* These curves have no noncuspidal places of degree1 or2 at p5. */
places:=[NoncuspidalPlaces(X,d):d in [1..6]];
assert #places[1] eq 0 and #places[2] eq 0;
if N eq 36 then assert [#s:s in places] eq [0,0,48,144,660,2856];end if;
printf "CLASS_CUSPIDAL_NONCUSP counts=%o time=%o\n",[#s:s in places],Cputime(t0);
pointclasses:=[[m2(Divisor(x)-d*Divisor(c0)):x in places[d]]:d in [1..6]];
candidates:={Cl|};
procedure TestPart(a,d,Hall,HQ,qbar,complements,qcuspidal,~candidates,~coarse,~fine,~qeffective)
 if a notin Hall then return;end if;
 coarse+:=1;
 tbl:=complements[6-d];key:=-qbar(a);
 if not IsDefined(tbl,key) then return;end if;
 fine+:=1;
 for b in tbl[key] do
  z:=a+b;assert z in HQ;
  if z in qcuspidal then qeffective+:=1;
  else Include(~candidates,z);end if;
 end for;
end procedure;
for d in [3..6] do
 coarse:=0;fine:=0;qeffective:=0;
 for a in pointclasses[d] do
  TestPart(a,d,Hall,HQ,qbar,complements,qcuspidal,~candidates,~coarse,~fine,~qeffective);
 end for;
 printf "CLASS_CUSPIDAL_TYPE type=[%o] total=%o Hall_pass=%o HQ_complement_pass=%o Qeffective_matches=%o cumulative_candidates=%o time=%o\n",
  d,#places[d],coarse,fine,qeffective,#candidates,Cputime(t0);
end for;
coarse:=0;fine:=0;qeffective:=0;total:=0;
for i in [1..#places[3]] do
 for j in [i..#places[3]] do
  a:=pointclasses[3][i]+pointclasses[3][j];total+:=1;
  TestPart(a,6,Hall,HQ,qbar,complements,qcuspidal,~candidates,~coarse,~fine,~qeffective);
 end for;
end for;
assert total eq Binomial(#places[3]+1,2);
printf "CLASS_CUSPIDAL_TYPE type=[3,3] total=%o Hall_pass=%o HQ_complement_pass=%o Qeffective_matches=%o cumulative_candidates=%o time=%o\n",
 total,coarse,fine,qeffective,#candidates,Cputime(t0);

if #candidates gt 0 then
 K:=Kernel(hQ);KB:=Matrix(Integers(),[Eltseq(Free!k):k in Generators(K)]);
 L:=LatticeWithBasis(LLL(KB));
 for z in candidates do
  v:=Eltseq(z@@hQ);w:=ClosestVector(L,Vector(Rationals(),v));
  short:=[v[i]-Integers()!w[i]:i in [1..#v]];
  assert hQ(Free!short) eq z;
  printf "CLASS_CUSPIDAL_CANDIDATE class=%o Qorbit_signatures=%o coefficients=%o\n",
   z,[CuspSignature(X,Support(O)[1]):O in orbits],short;
 end for;
end if;
printf "CLASS_CUSPIDAL_SIEVE_DONE N=%o p=%o D=6 remaining_classes=%o complete=true time=%o\n",N,p,#candidates,Cputime(t0);
quit;
