/* Direct exact certificate: evaluate the published model in a sextic
   field; no divisor normalization or finite-field reconstruction is needed. */
SetColumns(0);SetMemoryLimit(4*10^9);t0:=Cputime();
R<V>:=PolynomialRing(Rationals());
f:=V^6+4*V^5+41*V^4/6+31*V^3/6+4*V^2/3-V/3-1/6;
h:=(-92-15*V+542*V^2+1074*V^3+744*V^4+216*V^5)/17;
assert IsIrreducible(f);K<a>:=NumberField(f);u:=Evaluate(h,a);v:=a;
model:=Read("mdmagma/models_X1_m_n/X1_2_24.txt");
for line in Split(model) do
 val:=Split(Split(line,"=")[2],";")[1];
 if line[1] eq "X" then xx:=eval(val);assert xx eq 0;end if;
 if line[1] eq "q" then q:=eval(val);end if;
 if line[1] eq "t" then t:=eval(val);end if;
 if line[1] eq "E" then ainvs:=eval(val);end if;
 if line[1] eq "P" then pp:=eval(val);end if;
 if line[1] eq "Q" then qq:=eval(val);end if;
end for;
E:=EllipticCurve([K!x:x in ainvs]);P:=E!qq;Q:=E!pp;
assert Discriminant(E) ne 0;
assert 2*P eq E!0 and P ne E!0;
assert 24*Q eq E!0 and 12*Q ne E!0 and 8*Q ne E!0;
assert P ne 12*Q;
print "M2N24_SECOND_DIRECT_EXAMPLE_CONTAINS_DONE exact=true",[2,24],Cputime(t0);
print "M2N24_SECOND_DIRECT_FIELD",DefiningPolynomial(K);
print "M2N24_SECOND_DIRECT_AINVS",aInvariants(E);
print "M2N24_SECOND_DIRECT_P",P;
print "M2N24_SECOND_DIRECT_Q",Q;
T,tm:=TorsionSubgroup(E);
print "M2N24_SECOND_DIRECT_FULL_TORSION",Invariants(T),Cputime(t0);
assert Invariants(T) eq [2,24];
print "M2N24_SECOND_DIRECT_FULL_DONE exact=true",Cputime(t0);
quit;
