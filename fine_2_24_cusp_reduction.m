/* Exact fine rational-cuspidal-class sieve. Integral unit signatures
   give rational cusp blocks. All cusp complements and repeated [3,3]
   pairs are covered. The chosen degree-one reduction divisor affects
   only performance, not the class group or its exact maps. */
SetColumns(0);SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("class_rational_reduction.m");
Attach("fast_diamond_hash_m2.m");
t0:=Cputime();X:=MDX11(2,24,GF(5));C:=Curve(X);FF:=FunctionField(C);
b,c:=Explode([FF!f:f in X`_coordinates]);d:=b*(1-b-c);
E:=EllipticCurve([0,c,0,d,0]);Q:=E![b,b];fs:=[b,d,c^2-4*d];
for k in [1..12] do A:=k*Q;assert A[3] ne 0 and A[1] ne 0;Append(~fs,A[1]);end for;
cs:=Cusps(X);groups:=AssociativeArray();
for cusp in cs do
 sg:=[Valuation(f,cusp):f in fs];
 if IsDefined(groups,sg) then Append(~groups[sg],cusp);else groups[sg]:=[cusp];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));
blocks:=[&+[Divisor(cusp):cusp in groups[sg]]:sg in sigs];
degrees:=[Degree(o):o in blocks];
assert Sort(degrees) eq [1:i in [1..8]] cat [2:i in [1..6]] cat [4:i in [1..5]];
baseBlock:=Index(degrees,1);c0:=Support(blocks[baseBlock])[1];
print "N24_CUSP_REDUCTION_SIGNATURES",#blocks,Sort(degrees),Cputime(t0);
G,m1,m2:=ClassGroupRationalReduction(C,c0);
A:=JH(48,[25]);assert Dimension(A) eq Genus(X);
assert &*[a:a in Invariants(G)|a ne 0] eq Integers()!Evaluate(FrobeniusPolynomial(A,5),1);
gens:=[m2(o-Degree(o)*Divisor(c0)):o in blocks];
Free:=FreeAbelianGroup(#blocks);h:=hom<Free->G|gens>;H:=Image(h);
assert Invariants(H) eq [2,4,4,120,240];
QQ,qmap:=quo<G|H>;
print "N24_CUSP_REDUCTION_GROUP",Invariants(G),Invariants(H),Invariants(QQ),Cputime(t0);
cuspcosets:=[qmap(m2(Divisor(cusp)-Degree(cusp)*Divisor(c0))):cusp in cs];
ds:=[Degree(cusp):cusp in cs];
procedure AddCoset(remaining,start,current,~S)
 if remaining eq 0 then Include(~S,current);return;end if;
 for i in [start..#ds] do
  if ds[i] le remaining then AddCoset(remaining-ds[i],i,current+cuspcosets[i],~S);end if;
 end for;
end procedure;
complements:=[];
for k in [0..3] do
 S:={QQ!0};Exclude(~S,QQ!0);AddCoset(k,1,QQ!0,~S);
 Append(~complements,{-s:s in S});
 print "N24_CUSP_REDUCTION_COMPLEMENT_COSETS",k,#S;
end for;
assert [#s:s in complements] eq [1,5,23,75];
places:=[NoncuspidalPlaces(X,k):k in [1..6]];
assert [#xs:xs in places] eq [0,0,16,96,672,2712];
liftclasses:={G!0};Exclude(~liftclasses,G!0);allpass:=0;
for dd in [3..6] do
 rs:=PlacesUpToDiamondM2HashDirect(X,places[dd]);passed:=0;
 assert #rs eq [4,32,168,684][dd-2];
 print "N24_CUSP_REDUCTION_REPS",dd,#rs,Cputime(t0);
 for i in [1..#rs] do
  gc:=m2(Divisor(rs[i])-dd*Divisor(c0));cl:=qmap(gc);
  if cl in complements[7-dd] then
   passed+:=1;
   print "N24_CUSP_REDUCTION_SURVIVOR",dd,i,Eltseq(gc),Cputime(t0);
   assert dd eq 6;Include(~liftclasses,gc);
  end if;
  if i mod 25 eq 0 then print "N24_CUSP_REDUCTION_PROGRESS",dd,i,#rs,passed,Cputime(t0);end if;
 end for;
 allpass+:=passed;
 print "N24_CUSP_REDUCTION_TYPE",dd,#rs,passed,Cputime(t0);
end for;
rs3:=PlacesUpToDiamondM2HashDirect(X,places[3]);
classes3:=[m2(Divisor(x)-3*Divisor(c0)):x in places[3]];
passpair:=0;totalpair:=0;
for i in [1..#rs3] do
 gc:=m2(Divisor(rs3[i])-3*Divisor(c0));
 for j in [1..#places[3]] do
  totalpair+:=1;gpair:=gc+classes3[j];
  if qmap(gpair) eq QQ!0 then
   passpair+:=1;Include(~liftclasses,gpair);
   print "N24_CUSP_REDUCTION_PAIR_SURVIVOR",i,j,Eltseq(gpair),Cputime(t0);
  end if;
 end for;
end for;
assert totalpair eq 64;
print "N24_CUSP_REDUCTION_PAIR_TYPE",totalpair,passpair,Cputime(t0);
Kern:=Kernel(h);KB:=Matrix(Integers(),[Eltseq(Free!k):k in Generators(Kern)]);
L:=LatticeWithBasis(LLL(KB));vectors:=[];
for gc in liftclasses do
 v:=Eltseq(gc@@h);w:=ClosestVector(L,Vector(Rationals(),v));
 short:=[v[i]-Integers()!w[i]:i in [1..#v]];
 assert h(Free!short) eq gc;Append(~vectors,short);
end for;
print "N24_LIFT_CODE_BEGIN";
print [*24,6,sigs,degrees,baseBlock,vectors*];
print "N24_LIFT_CODE_END";
print "N24_CUSP_REDUCTION_DONE",allpass+passpair,#liftclasses,"complete=true",Cputime(t0);
quit;
