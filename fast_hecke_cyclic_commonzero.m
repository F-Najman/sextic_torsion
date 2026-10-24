/* Attach after mdmagma, fast_hecke.m and
   fast_hecke_11.m.  The latter supplies
   MDPlaceCommonZeros; the cyclic kernel routine is unchanged. */

intrinsic ModuliPointCyclicCommonZeros(X::MDCrvMod1,E::CrvEll,Ls::Rec) -> PlcCrvElt
{ Recover a cyclic moduli place by factoring the common coordinate-zero ideal }
    xy:=TateNormalForm_xy(E,Ls`P);
    F:=FunctionField(Curve(X));
    return MDPlaceCommonZeros([F.1,F.2],xy);
end intrinsic;

intrinsic HeckeOperatorCommonZeros(X::MDCrvMod1,q::RngIntElt,x::PlcCrvElt) -> DivCrvElt
{ The verified cyclic fast kernel-orbit algorithm with common-zero place reconstruction }
    require IsPrime(q) and IsOdd(q) and GCD(q,Level(X)*Characteristic(BaseRing(X))) eq 1 : "q must be odd prime coprime to level and characteristic";
    E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
    isog:=MDIsogeniesFast(E,P,q);
    D:=DivisorGroup(Curve(X))!0;
    for t in isog do
        EG,PG,e:=Explode(t);
        y:=ModuliPointCyclicCommonZeros(X,EG,rec<recformat<P:PtEll>|P:=PG>);
        mult:=e*Degree(x)/Degree(y);assert Denominator(mult) eq 1;
        D+:=(Integers()!mult)*Divisor(y);
    end for;
    assert Degree(D) eq (q+1)*Degree(x);
    return D;
end intrinsic;

intrinsic HeckeOperatorCommonZeros(X::MDCrvMod1,q::RngIntElt,D::DivCrvElt) -> DivCrvElt
{ Common-zero Hecke operator on a divisor supported on noncuspidal places }
    supp,mults:=Support(D);answer:=DivisorGroup(Curve(X))!0;
    for i in [1..#supp] do answer+:=mults[i]*HeckeOperatorCommonZeros(X,q,supp[i]);end for;
    return answer;
end intrinsic;
