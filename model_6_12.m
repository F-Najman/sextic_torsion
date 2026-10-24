/* Restore the archived X1(6,12) model, whose extra level point is
   supplied only by its x-coordinate.  The finite-field level structure
   is checked before use.  Attach after fast_hecke_11.m. */

intrinsic MDX612(F::FldFin) -> MDCrvMod11
{ Construct the archived X1(6,12) curve over a finite field containing sixth roots of unity }
 require Characteristic(F) notin [2,3] and (#F-1) mod 6 eq 0 : "Require characteristic prime to6 and sixth roots of unity";
 A<u,v>:=AffineSpace(F,2);
 data:=Split(Read("mdmagma/models_X1_m_n/X1_6_12.txt"));
 for line in data do
  val:=Split(Split(line,"=")[2],";")[1];
  if line[1] eq "X" then equation:=eval(val);end if;
  if line[1] eq "r" then r:=eval(val);end if;
  if line[1] eq "s" then s:=eval(val);end if;
  if line[1] eq "q" then q:=eval(val);end if;
  if line[1] eq "E" then es:=eval(val);end if;
 end for;
 C:=ProjectiveClosure(Curve(A,equation));FF:=FunctionField(C);
 X:=New(MDCrvMod11);X`M:=6;X`N:=12;
 X`_P:=[FF!q];X`_Q:=[FF!0,FF!0];X`_coordinates:=[FF.1,FF.2];
 _InitMDCrvMod(X,12,C,[FF!a:a in es]);
 return X;
end intrinsic;

intrinsic LevelStructure612(X::MDCrvMod11,x::PlcCrvElt) -> Rec
{ Recover and verify the independent order-six point from its archived x-coordinate }
 require X`M eq 6 and X`N eq 12 : "Require X1(6,12)";
 E:=EllipticCurve(X,x);K:=BaseRing(E);Q:=E![0,0];
 assert Order(Q) eq 12;
 xx:=Evaluate(X`_P[1],x);
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));
 R<Y>:=PolynomialRing(K);
 ys:=Roots(Y^2+(a1*xx+a3)*Y-xx^3-a2*xx^2-a4*xx-a6);
 assert #ys eq 2;
 P:=E![xx,ys[1][1]];
 assert Order(P) eq 6;
 z:=WeilPairing(P,2*Q,6);
 assert Order(z) eq 6;
 z0:=K!(PrimitiveElement(BaseRing(X))^((#BaseRing(X)-1) div 6));
 if z ne z0 then P:=-P;end if;
 assert WeilPairing(P,2*Q,6) eq z0;
 return rec<recformat<P:PtEll,Q:PtEll>|P:=P,Q:=Q>;
end intrinsic;

intrinsic Diamond612(X::MDCrvMod11,a::RngIntElt,x::PlcCrvElt) -> PlcCrvElt
{ Diamond operator on the restored X1(6,12) model }
 require GCD(a,12) eq 1 : "Require a coprime to12";
 Ls:=LevelStructure612(X,x);E:=Curve(Ls`P);
 La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*Ls`P,Q:=a*Ls`Q>;
 return ModuliPoint612(X,E,La);
end intrinsic;

intrinsic Representatives612(X::MDCrvMod11,xs::SeqEnum[PlcCrvElt]) -> SeqEnum[PlcCrvElt]
{ Exact diamond orbit representatives; units modulo12 up to simultaneous sign are1,5 }
 left:=Seqset(xs);rs:=[];
 while #left gt 0 do
  x:=Rep(left);Append(~rs,x);Exclude(~left,x);
  y:=Diamond612(X,5,x);assert y in xs;
  Exclude(~left,y);
 end while;
 return rs;
end intrinsic;

intrinsic HeckeOperator612(X::MDCrvMod11,q::RngIntElt,x::PlcCrvElt) -> DivCrvElt
{ Kernel-orbit Hecke operator for the restored X1(6,12) model }
 require IsPrime(q) and IsOdd(q) and GCD(q,12*Characteristic(BaseRing(X))) eq 1 : "Require odd prime q coprime to level and characteristic";
 E:=EllipticCurve(X,x);Ls:=LevelStructure612(X,x);
 isog:=MDIsogeniesFastPoints(E,[Ls`P,Ls`Q],q);
 D:=DivisorGroup(Curve(X))!0;
 for t in isog do
  EG,points,e:=Explode(t);
  // Multiplication by q^-1 on P fixes the chosen pairing component.
  a:=InverseMod(q,6);
  La:=rec<recformat<P:PtEll,Q:PtEll>|P:=a*points[1],Q:=points[2]>;
  y:=ModuliPoint612(X,EG,La);
  mult:=e*Degree(x)/Degree(y);assert Denominator(mult) eq 1;
  D+:=(Integers()!mult)*Divisor(y);
 end for;
 assert Degree(D) eq (q+1)*Degree(x);
 return D;
end intrinsic;

intrinsic ModuliPoint612(X::MDCrvMod11,E::CrvEll,Ls::Rec) -> PlcCrvElt
{ Rational inverse on the archived X1(6,12) model }
 require X`M eq 6 and X`N eq 12 : "Require X1(6,12)";
 P:=Ls`P;Q:=Ls`Q;assert 6*P eq E!0 and 12*Q eq E!0;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));
 x0:=Q[1];y0:=Q[2];
 aa3:=2*y0+a3+a1*x0;aa2:=3*x0+a2;
 aa4:=3*x0^2+2*x0*a2+a4-a1*y0;
 slope:=aa4/aa3;
 aaa1:=2*slope+a1;aaa2:=aa2-slope^2-a1*slope;
 scale:=aaa2/aa3;
 b:=-aaa2^3/aa3^2;c:=1-aaa1*scale;
 r:=b/c;s:=c^2/(b-c);
 u:=(r-1)/(r-s);
 q:=scale^2*(P[1]-x0);
 v:=q*u^2*(u-1)^2/(u+1)-u^2-1;
 assert r eq (u^2+1)/(u^2-u) and s eq (u^2+u+1)/u^2;
 return MDPlaceCommonZeros(X`_coordinates,[u,v]);
end intrinsic;
