/* S3 decomposition: J(source) is isogenous to J(A3-quotient) times
   J1(30)^2 divided by J1(15)^2. This gives independent point orders. */
SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
S:=JH(60,[31]);A:=JOne(15);B:=JOne(30);
assert Dimension(S) eq 25 and Dimension(A) eq 1 and Dimension(B) eq 9;
gg:=0;
for p in [7,11,13,17,19] do
 fs:=FrobeniusPolynomial(S,p);fa:=FrobeniusPolynomial(A,p);fb:=FrobeniusPolynomial(B,p);
 fy:=(fs*fa^2) div fb^2;assert fy*fb^2 eq fs*fa^2 and Degree(fy) eq 18;
 nn:=Integers()!Evaluate(fy,1);assert nn mod 256 eq 0;gg:=GCD(gg,nn);
 if p eq 7 then assert nn eq 69249024;end if;
 print "A3_M2_30_FROBENIUS",p,fy,nn,Factorization(nn),"gcd",gg,Cputime(t0);
end for;
print "A3_M2_30_FROBENIUS_DONE exact=true gcd_order",gg,Cputime(t0);
quit;
