/* Independent full-torsion upper bounds from good finite reductions.
   No number-field TorsionSubgroup call is used in this certificate. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
R<V>:=PolynomialRing(Rationals());
fs:=[V^6-3*V^5+V^4/2+4*V^3-V^2-3*V/2+1/2,
 V^6+4*V^5+41*V^4/6+31*V^3/6+4*V^2/3-V/3-1/6];
hs:=[3-3*V-13*V^2+2*V^3+10*V^4-4*V^5,
 (-92-15*V+542*V^2+1074*V^3+744*V^4+216*V^5)/17];
for example in [1,2] do
 f:=fs[example];assert IsIrreducible(f);K<a>:=NumberField(f);
 u:=Evaluate(hs[example],a);v:=a;
 q:=(u+v+2)/(u-v);
 t:=(-u^3-u^2*v-2*u^2+u*v^2+v^3+2*v^2)/(u^3+u^2*v+2*u^2+u*v^2+2*u*v+v^3+2*v^2-2);
 E:=EllipticCurve([K!0,t^2-2*q*t-2,K!0,-(t^2-1)*(q*t+1)^2,K!0]);
 P:=E![0,0];Q:=E![(t+1)*(q*t+1),t*(q*t+1)*(t+1)];
 assert 2*P eq E!0 and P ne E!0;
 assert 24*Q eq E!0 and 12*Q ne E!0 and 8*Q ne E!0 and P ne 12*Q;
 RK<XK>:=PolynomialRing(K);curves:=[*E*];images:=[* *];
 for j in [1..3] do
  B:=[P,12*Q,P+12*Q][j];E2,phi:=IsogenyFromKernel(E,XK-B[1]);
  assert Degree(phi) eq 2 and phi(B) eq E2!0;
  Q2:=phi(Q);P2:=phi(P);
  if j eq 2 then
   assert 12*Q2 eq E2!0 and 6*Q2 ne E2!0 and 4*Q2 ne E2!0;
   assert 2*P2 eq E2!0 and P2 ne E2!0 and P2 ne 6*Q2;
  else assert 24*Q2 eq E2!0 and 12*Q2 ne E2!0 and 8*Q2 ne E2!0;end if;
  Append(~curves,E2);
  Append(~images,[P2,Q2]);
 end for;
 // Isogenous finite curves have equal orders. The reduction bound is48
 // for every target; exact halving below removes the residual factor2.
 expected:=[48,48,48,48];bounds:=[0,0,0,0];disc:=Discriminant(f);
 for p in PrimesInInterval(5,97) do
  if Numerator(disc) mod p eq 0 or Denominator(disc) mod p eq 0 then continue;end if;
  Fp:=GF(p);S<T>:=PolynomialRing(Fp);fp:=S!f;
  for fac in Factorization(fp) do
   ff:=fac[1];dd:=Degree(ff);if dd gt 3 then continue;end if;
   if dd eq 1 then L:=Fp;aa:=-Coefficient(ff,0)/Coefficient(ff,1);
   else L<aa>:=ext<Fp|ff>;end if;
   for k in [1..4] do
    if bounds[k] eq expected[k] then continue;end if;
    coeffs:=[Eltseq(z):z in aInvariants(curves[k])];
    if &or[Denominator(z) mod p eq 0:z in &cat coeffs] then continue;end if;
    ainvs:=[&+[L!Coefficient(S![Fp!z:z in co],j)*aa^j:j in [0..#co-1]]:co in coeffs];
    Ep:=EllipticCurve(ainvs);if Discriminant(Ep) eq 0 then continue;end if;
    n:=#Ep;assert n mod expected[k] eq 0;bounds[k]:=GCD(bounds[k],n);
    print "N24_GOOD_REDUCTION_BOUND",example,k,p,dd,n,bounds[k];
   end for;
  end for;
  if bounds eq expected then break;end if;
 end for;
 assert bounds eq expected;
 for j in [1..3] do
  E2:=curves[j+1];P2,Q2:=Explode(images[j]);
  if j eq 2 then
   assert #Roots(DivisionPolynomial(E2,2)) eq 3;
   halfP:=DivisionPoints(P2,2);
   halfQ:=DivisionPoints(3*Q2,2);
   halfPQ:=DivisionPoints(3*Q2+P2,2);
   assert #halfP eq 0 and #halfQ eq 0 and #halfPQ eq 0;
   print "N24_EXACT_HALVING_COUNTS",example,j,[#halfP,#halfQ,#halfPQ];
  else
   assert #Roots(DivisionPolynomial(E2,2)) eq 1;
   halves:=DivisionPoints(Q2,2);assert #halves eq 0;
   print "N24_EXACT_HALVING_COUNTS",example,j,[#halves];
  end if;
 end for;
 print "N24_REDUCTION_BOUND_EXAMPLE_DONE",example,bounds,Cputime(t0);
end for;
print "N24_REDUCTION_BOUNDS_ALL_DONE exact=true",Cputime(t0);quit;
