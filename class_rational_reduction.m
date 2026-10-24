/* Identical map construction to Magma's ClassGroup(C::Crv[FldFin]),
   but with a chosen degree-one divisor instead of (x)_infinity.
   Validated on X1(21)/F5 against Frobenius order and independent
   principal-divisor checks; test_class_rational_reduction_v3.log. */
intrinsic ClassGroupRationalReduction(C::Crv,c0::PlcCrvElt : Proof:=true) -> GrpAb,Map,Map
{ Finite-curve divisor class group using the rational place c0 for reduction }
 require Type(BaseRing(C)) eq FldFin and Degree(c0) eq 1 : "Require a finite base field and a degree-one place";
 F:=AlgorithmicFunctionField(FunctionField(C));
 G,_,phi:=ClassGroup(F:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=Proof);
 m2:=map<DivisorGroup(C)->G|D:->phi(FunctionFieldDivisor(D)),g:->CurveDivisor(C,g@@phi)>;
 return G,Inverse(m2),m2;
end intrinsic;
