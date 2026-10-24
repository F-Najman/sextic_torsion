SetColumns(0);
SetMemoryLimit(2500000000);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
Attach("model_6_12.m");
t0:=Cputime();X:=MDX612(GF(7));
assert #NoncuspidalPlaces(X,1) eq 0 and #NoncuspidalPlaces(X,2) eq 0;
xs:=NoncuspidalPlaces(X,3);rs:=Representatives612(X,xs);
assert #xs eq 96 and #rs eq 48;
print "SIEVE_6_12_SETUP",Genus(X),#xs,#rs,Cputime(t0);
survivors:=[];
for i in [1..#rs] do
 x:=rs[i];D:=Divisor(x);killed:=false;
 for q in [13,19] do
  T:=HeckeOperator612(X,q,x);
  A:=T-D-q*Divisor(Diamond612(X,q,x));
  principal:=IsPrincipal(A);
  print "SIEVE_6_12_TEST",i,q,principal,Cputime(t0);
  if not principal then killed:=true;break;end if;
 end for;
 if not killed then Append(~survivors,i);end if;
end for;
print "SIEVE_6_12_DONE",#rs,#survivors,survivors,Cputime(t0);
quit;
