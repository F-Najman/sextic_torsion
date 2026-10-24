/* Cyclic Tate-normal-form hash using the explicit root-defined common
   finite-field map from the independently validated M2 hash helpers.
   Attach after mdmagma and fast_diamond_hash_m2.m. */
intrinsic PlacesUpToDiamondTateHashRoots(X::MDCrvMod1,S::SeqEnum[PlcCrvElt]) -> SeqEnum[PlcCrvElt]
{ Exact cyclic diamond representatives using explicit finite-field maps }
 if #S eq 0 then return S;end if;
 p:=#BaseRing(X);assert IsPrime(p);d:=Degree(S[1]);
 assert &and[Degree(x) eq d:x in S];
 K:=BaseRing(EllipticCurve(X,S[1]));assert #K eq p^d;
 n:=Level(X);units:=[a:a in [1..n div 2]|GCD(a,n) eq 1];
 assert #units eq EulerPhi(n) div 2;
 seen:=AssociativeArray();reps:=[];
 for x in S do
  E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
  key:=MDM2FrobeniusKey(TateNormalForm_xy(E,P),K,p,d);
  if IsDefined(seen,key) then continue;end if;
  Append(~reps,x);
  for a in units do
   ky:=MDM2FrobeniusKey(TateNormalForm_xy(E,a*P),K,p,d);
   if IsDefined(seen,ky) then assert seen[ky] eq #reps;else seen[ky]:=#reps;end if;
  end for;
  assert IsDefined(seen,key) and seen[key] eq #reps;
 end for;
 assert #Keys(seen) eq #S;return reps;
end intrinsic;
