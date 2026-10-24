/* Complete finite-field portion of the C54 proof via X1(27).
   No finite ClassGroup or Q-normalization is used in this certificate. */
SetColumns(0);SetMemoryLimit(2*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
// These are Q-cusp orbit signatures and degrees.  Their modular-unit
// valuations are preserved at2 and5; only rational/quadratic cusps occur.
ss:=[[25,-23],[-19,2],[-1,2],[-5,1],[-11,4],[-23,1],[5,8]];
ds:=[1,1,2,2,1,1,1];vs:=[10,-1,-1,-1,1,-1,1];
for p in [2,5] do
 X:=MDX1(27,GF(p));C:=Curve(X);assert Genus(C) eq 13;
 os:=[Type(O) eq PlcCrvElt select Divisor(O) else O:O in CuspOrbitsQ(X)];
 sigs:=[CuspSignature(X,Support(O)[1]):O in os];
 assert #Seqset(sigs) eq #os and &and[s in sigs:s in ss];
 cs:=[os[Index(sigs,s)]:s in ss];assert [Degree(O):O in cs] eq ds;
 L:=&+[vs[i]*cs[i]:i in [1..7]];assert Degree(L) eq 6;
 V,vm:=RiemannRochSpace(L);assert Dimension(V) eq 1;
 T:=L+Divisor(vm(V.1));assert IsEffective(T) and Degree(T) eq 6;
 sup,mul:=Support(T);
 printf "C54_CERT_RR p=%o dimension=%o degrees=%o multiplicities=%o time=%o\n",p,Dimension(V),[Degree(s):s in sup],mul,Cputime(t0);
 assert &and[not IsCusp(X,s):s in sup];
 if p eq 2 then
  assert #sup eq 1 and Degree(sup[1]) eq 6 and mul eq [1];
  xs:=[NoncuspidalPlaces(X,d):d in [1..6]];
  assert [#s:s in xs] eq [0,0,0,0,0,9];
  // Use the original reference diamond correspondence, independently
  // of the orbit hash that discovered the single orbit.
  orb:=Seqset(DiamondOrbit(X,sup[1]));
  assert #orb eq 9 and orb eq Seqset(xs[6]);
  printf "C54_CERT_F2_PLACES counts=%o reference_diamond_orbit_size=%o time=%o\n",[#s:s in xs],#orb,Cputime(t0);
 else
  assert Sort([Degree(s):s in sup]) eq [2,4] and mul eq [1,1];
  x:=[s:s in sup|Degree(s) eq 2][1];E:=EllipticCurve(X,x);
  assert #E eq 27 and Invariants(AbelianGroup(E)) eq [27];
  assert MinimalPolynomial(jInvariant(E),GF(5)) eq Polynomial(GF(5),[4,2,1]);
  print "C54_CERT_F5_ODD_GROUP",#E,Invariants(AbelianGroup(E));
 end if;
end for;
printf "C54_FINITE_CERTIFICATE_DONE complete=true time=%o\n",Cputime(t0);
quit;
