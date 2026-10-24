/* Attach after the verified cyclic common-zero inverse.  Enumerate all
   three degree-two isogenies, then descend their multiset of image places. */
intrinsic HeckeTwoCommonZeros(X::MDCrvMod1,x::PlcCrvElt) -> DivCrvElt
{ Compute T2 exactly using the original full isogeny enumeration }
 require IsOdd(Characteristic(BaseRing(X))) and IsOdd(Level(X)) : "Require odd characteristic and level";
 E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
 phis:=MDIsogenies(E,2);assert #phis eq 3;
 ys:=[];ms:=[];
 for phi in phis do
  y:=ModuliPointCyclicCommonZeros(X,Codomain(phi),rec<recformat<P:PtEll>|P:=phi(P)>);
  k:=Index(ys,y);
  if k eq 0 then Append(~ys,y);Append(~ms,1);else ms[k]+:=1;end if;
 end for;
 D:=DivisorGroup(Curve(X))!0;
 for k in [1..#ys] do
  mult:=ms[k]*Degree(x)/Degree(ys[k]);assert Denominator(mult) eq 1;
  D+:=(Integers()!mult)*Divisor(ys[k]);
 end for;
 assert Degree(D) eq 3*Degree(x);
 return D;
end intrinsic;
