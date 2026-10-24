/* Canonical normalization for (E,P2,Q_N), N>=6:
   E:y^2=x^3+c*x^2+d*x, P2=(0,0), Q=(b,b).
   For a unit a and (xa,ya)=aQ, normalize using u=ya/xa.
   The resulting coordinates are b_a=xa^3/ya^2,
   c_a=c*xa^2/ya^2. The pair (b_a,c_a) uniquely determines
   the isomorphism class of the ordered level structure.
   Frobenius-canonical keys therefore classify exact closed places.
   This file uses no moduli-place reconstruction or ideal factoring.
   Explicit root-defined field maps avoid Magma's very slow generic
   Embed(L,K) on some equal-size residue-field presentations. */

intrinsic MDM2ToCommonField(bc::SeqEnum,K::FldFin) -> SeqEnum
{ Convert an absolute finite-field pair by an explicit root-defined map }
 L:=Parent(bc[1]);if L cmpeq K then return bc;end if;
 assert Characteristic(L) eq Characteristic(K) and #L eq #K;
 if Degree(L) eq 1 then return [K!(Integers()!z):z in bc];end if;
 f:=DefiningPolynomial(L,PrimeField(L));assert Degree(f) eq Degree(L);
 fk:=PolynomialRing(K)![K!(Integers()!a):a in Coefficients(f)];
 rr:=Roots(fk);assert #rr eq Degree(L);
 mp:=hom<L->K|rr[1][1]>;
 return [mp(z):z in bc];
end intrinsic;

intrinsic MDM2FrobeniusKey(bc::SeqEnum,K::FldFin,p::RngIntElt,d::RngIntElt) -> SeqEnum
{ Canonical prime-field coefficient key of the Frobenius orbit of (b,c) }
 assert #bc eq 2 and #K eq p^d;
 v:=MDM2ToCommonField(bc,K);
 keys:=[];
 for i in [1..d] do
  key:=&cat[[Integers()!c:c in Eltseq(z,PrimeField(K))]:z in v];
  assert #key eq 2*d;Append(~keys,key);v:=[z^p:z in v];
 end for;
 return Sort(keys)[1];
end intrinsic;

intrinsic MDM2CanonicalPair(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ Evaluate only the two canonical parameters, into one common residue field }
 require X`M eq 2 : "Require M=2";
 bc:=[Evaluate(f,x):f in X`_coordinates];
 return MDM2ToCommonField(bc,K);
end intrinsic;

intrinsic MDM2DiamondKeysDirect(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ Exact diamond keys using only two function evaluations at x }
 p:=#BaseRing(X);dd:=Degree(x);bc:=MDM2CanonicalPair(X,x,K);
 b,c:=Explode(bc);E:=EllipticCurve([K!0,c,K!0,b*(1-b-c),K!0]);Q:=E![b,b];
 keys:=[];n:=Level(X);
 for a in [1..n div 2] do
  if GCD(a,n) ne 1 then continue;end if;
  Qa:=a*Q;assert Qa[3] ne 0 and Qa[1] ne 0 and Qa[2] ne 0;
  Append(~keys,MDM2FrobeniusKey([Qa[1]^3/Qa[2]^2,c*Qa[1]^2/Qa[2]^2],K,p,dd));
 end for;
 assert #keys eq EulerPhi(n) div 2;return keys;
end intrinsic;

intrinsic PlacesUpToDiamondM2HashDirect(X::MDCrvMod11,S::SeqEnum[PlcCrvElt] : Trace:=false) -> SeqEnum[PlcCrvElt]
{ Exact M2 hash, evaluating only(b,c) and sharing one finite-field parent }
 if #S eq 0 then return S;end if;
 require X`M eq 2 : "Require M=2";
 p:=#BaseRing(X);assert IsPrime(p);
 dd:=Degree(S[1]);assert &and[Degree(x) eq dd:x in S];
 K:=Parent(Evaluate(X`_coordinates[1],S[1]));assert #K eq p^dd;
 seen:=AssociativeArray();reps:=[];
 for idx in [1..#S] do
  x:=S[idx];if Trace then print "M2_HASH_INPUT",idx;end if;
  bc:=MDM2CanonicalPair(X,x,K);if Trace then print "M2_HASH_EVALUATED",idx;end if;
  key:=MDM2FrobeniusKey(bc,K,p,dd);if Trace then print "M2_HASH_KEY",idx;end if;
  if IsDefined(seen,key) then continue;end if;
  Append(~reps,x);
  for ky in MDM2DiamondKeysDirect(X,x,K) do
   if IsDefined(seen,ky) then assert seen[ky] eq #reps;
   else seen[ky]:=#reps;end if;
  end for;
  if Trace then print "M2_HASH_DIAMONDS",idx,#reps;end if;
  assert IsDefined(seen,key) and seen[key] eq #reps;
 end for;
 assert #Keys(seen) eq #S;return reps;
end intrinsic;

intrinsic MDM2DiamondKeys(X::MDCrvMod11,x::PlcCrvElt,K::FldFin) -> SeqEnum
{ Frobenius keys of all diamond translates of the noncuspidal M2 place x }
 require X`M eq 2 : "Require M=2";
 p:=#BaseRing(X);assert IsPrime(p);
 dd:=Degree(x);assert #K eq p^dd;
 E:=EllipticCurve(X,x);Ls:=LevelStructure(X,x);
 ai:=aInvariants(E);assert ai[1] eq 0 and ai[3] eq 0 and ai[5] eq 0;
 assert Ls`P eq E![0,0] and Ls`Q[1] eq Ls`Q[2];
 c:=ai[2];n:=Level(X);keys:=[];
 for a in [1..n div 2] do
  if GCD(a,n) ne 1 then continue;end if;
  Q:=a*Ls`Q;assert Q[3] ne 0 and Q[1] ne 0 and Q[2] ne 0;
  bc:=[Q[1]^3/Q[2]^2,c*Q[1]^2/Q[2]^2];
  Append(~keys,MDM2FrobeniusKey(bc,K,p,dd));
 end for;
 assert #keys eq EulerPhi(n) div 2;
 return keys;
end intrinsic;

intrinsic PlacesUpToDiamondM2Hash(X::MDCrvMod11,S::SeqEnum[PlcCrvElt]) -> SeqEnum[PlcCrvElt]
{ Exact M2 diamond representatives from canonical (b,c) Frobenius keys }
 if #S eq 0 then return S;end if;
 require X`M eq 2 : "Require M=2";
 p:=#BaseRing(X);assert IsPrime(p);
 dd:=Degree(S[1]);assert &and[Degree(x) eq dd:x in S];
 K:=BaseRing(EllipticCurve(X,S[1]));assert #K eq p^dd;
 seen:=AssociativeArray();reps:=[];
 for x in S do
  E:=EllipticCurve(X,x);Q:=LevelStructure(X,x)`Q;
  key:=MDM2FrobeniusKey([Q[1],aInvariants(E)[2]],K,p,dd);
  if IsDefined(seen,key) then continue;end if;
  Append(~reps,x);
  for ky in MDM2DiamondKeys(X,x,K) do
   if IsDefined(seen,ky) then assert seen[ky] eq #reps;
   else seen[ky]:=#reps;end if;
  end for;
  assert IsDefined(seen,key) and seen[key] eq #reps;
 end for;
 assert #Keys(seen) eq #S;
 return reps;
end intrinsic;
