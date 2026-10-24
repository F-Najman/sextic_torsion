/* Matched integral-unit cusp blocks and candidate-specific multiples.
   Any pair retained is only a necessary local condition, not a point. */
SetColumns(0);SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
Attach("class_rational_reduction.m");
t0:=Cputime();Gs:=[* *];ms:=[* *];Hs:=[* *];hs:=[* *];
signatures:=[];degrees:=[];states:=[* *];
Free:=FreeAbelianGroup(24);
function ExactReps(Y,xs)
 left:=Seqset(xs);rs:=[];
 while #left gt 0 do
  x:=Rep(left);Append(~rs,x);
  for a in [1,5,7] do
   y:=DiamondOperator(Y,a,x);assert y in xs;Exclude(~left,y);
  end for;
 end while;
 return rs;
end function;
for p in [7,13] do
 X:=MDX11(3,18,GF(p));FF:=FunctionField(Curve(X));
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0 and P eq E![0,0];
 fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..9] do
  for j in [0..2] do
   A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
   Append(~fs,A[1]);Append(~fs,A[2]);
  end for;
 end for;
 for k in [1..9] do
  A:=k*Q;B:=k*Q+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
 end for;
 for k in [1..8] do
  A:=k*Q+P;B:=k*Q+2*P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
 end for;
 groups:=AssociativeArray();cs:=Cusps(X);
 for cusp in cs do
  sg:=[Valuation(f,cusp):f in fs];
  if IsDefined(groups,sg) then Append(~groups[sg],cusp);else groups[sg]:=[cusp];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 blocks:=[&+[Divisor(cusp):cusp in groups[sg]]:sg in sigs];
 ds:=[Degree(o):o in blocks];assert Sort(ds) eq [1:i in [1..18]] cat [3:i in [1..6]];
 assert &and[#groups[sg] eq 1:sg in sigs];
 if p eq 7 then signatures:=sigs;degrees:=ds;
 else assert sigs eq signatures and ds eq degrees;end if;
 baseBlock:=Index(ds,1);c0:=Support(blocks[baseBlock])[1];
 G,m1,m2:=ClassGroupRationalReduction(Curve(X),c0);
 gs:=[m2(blocks[j]-ds[j]*Divisor(c0)):j in [1..#blocks]];
 h:=hom<Free->G|gs>;H:=Image(h);
 assert Invariants(H) eq [3,3,3,3,63,63,63];
 Append(~Gs,G);Append(~ms,m2);Append(~Hs,H);Append(~hs,h);
 print "M3N18_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Cputime(t0);
 qp:=p eq 7 select 13 else 7;survplaces:={};surv2:=0;
 for dd in [2,3] do
  xs:=NoncuspidalPlaces(X,dd);rs:=ExactReps(X,xs);passed:=0;
  if p eq 7 then assert #rs eq (dd eq 2 select 9 else 48);
  else assert #rs eq (dd eq 2 select 18 else 282);end if;
  for i in [1..#rs] do
   x:=rs[i];D:=Divisor(x);
   T:=HeckeOperatorFast(X,qp,x);
   principal:=IsPrincipal(T-D-qp*Divisor(DiamondOperator(X,qp,x)));
   if principal then
    passed+:=1;
    print "M3N18_DIAGONAL_HECKE_SURVIVOR",p,dd,i,Cputime(t0);
    if dd eq 2 then surv2+:=1;
    else
     for a in [1,5,7] do Include(~survplaces,DiamondOperator(X,a,x));end for;
    end if;
   end if;
   if i mod 25 eq 0 then print "M3N18_DIAGONAL_HECKE_PROGRESS",p,dd,i,#rs,passed,Cputime(t0);end if;
  end for;
  assert passed eq (dd eq 2 select 0 else 12);
  print "M3N18_DIAGONAL_HECKE_TYPE",p,dd,#rs,passed,Cputime(t0);
 end for;
 assert surv2 eq 0 and #survplaces eq 36;
 cls:=[m2(Divisor(x)-3*Divisor(c0)):x in survplaces];Append(~states,cls);
 print "M3N18_DIAGONAL_STATES",p,#cls,Cputime(t0);
end for;
G7:=Gs[1];G13:=Gs[2];h7:=hs[1];h13:=hs[2];H7:=Hs[1];H13:=Hs[2];
assert Invariants(quo<G7|H7>) eq [3,3,9,0];
for k in Generators(Kernel(h7)) do assert h13(Free!k) eq G13!0;end for;
for k in Generators(Kernel(h13)) do assert h7(Free!k) eq G7!0;end for;
print "M3N18_DIAGONAL_CUSP_KERNELS_MATCH",#H7,#H13,Cputime(t0);
QQ7,q7:=quo<G7|H7>;matches9:=0;matchesminimal:=0;orders:=[];
for i in [1..#states[1]] do
 a:=states[1][i];assert 9*a in H7;t:=Order(q7(a));assert t in [3,9];Append(~orders,t);
 predict9:=h13((9*a)@@h7);predictmin:=h13((t*a)@@h7);
 for j in [1..#states[2]] do
  b:=states[2][j];
  if predict9 eq 9*b then
   matches9+:=1;
   if predictmin eq t*b then
    matchesminimal+:=1;
    print "M3N18_DIAGONAL_PAIR_SURVIVOR",i,j,t,Eltseq(a),Eltseq(b);
   end if;
  end if;
 end for;
end for;
print "M3N18_DIAGONAL_QUOTIENT_ORDERS",orders;
print "M3N18_ARITHMETIC_CODE_BEGIN";
print [* [7,13],
 [Integers()!RootOfUnity(3,GF(7)),Integers()!RootOfUnity(3,GF(13))],
 signatures,degrees,Index(degrees,1),
 Invariants(G7),[Eltseq(r):r in Rows(RelationMatrix(G7))],
 Invariants(G13),[Eltseq(r):r in Rows(RelationMatrix(G13))],
 [Eltseq(h7(Free.i)):i in [1..24]],
 [Eltseq(h13(Free.i)):i in [1..24]],
 [Eltseq(a):a in states[1]],[Eltseq(a):a in states[2]] *];
print "M3N18_ARITHMETIC_CODE_END";
print "M3N18_DIAGONAL_DONE",#states[1],#states[2],matches9,matchesminimal,"complete=true",Cputime(t0);
quit;
