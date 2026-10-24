SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke_11.m");
t0:=Cputime();X:=MDX11(4,20,GF(9));
assert #NoncuspidalPlaces(X,1) eq 0;
survivors:=[];
for d in [2,3] do
 xs:=NoncuspidalPlaces(X,d);left:=Seqset(xs);rs:=[];
 while #left gt 0 do
  x:=Rep(left);Append(~rs,x);
  for a in [1,3,7,9] do
   y:=DiamondOperator(X,a,x);assert y in xs;Exclude(~left,y);
  end for;
 end while;
 assert <d,#xs,#rs> in [<2,96,24>,<3,320,80>];
 qs:=(d eq 2) select [29,41] else [13,17];
 if d eq 2 then assert &and[q mod 20 in [1,9]:q in qs];end if;
 print "SIEVE_4_20_P3_SETUP",d,Genus(X),#xs,#rs,qs,Cputime(t0);
 for i in [1..#rs] do
  x:=rs[i];D:=Divisor(x);killed:=false;
  for q in qs do
   T:=HeckeOperatorFast(X,q,x);
   A:=T-D-q*Divisor(DiamondOperator(X,q,x));
   principal:=IsPrincipal(A);
   print "SIEVE_4_20_P3_TEST",d,i,q,principal,Cputime(t0);
   if not principal then killed:=true;break;end if;
  end for;
  if not killed then Append(~survivors,<d,i>);end if;
 end for;
 print "SIEVE_4_20_P3_TYPE_DONE",d,#rs,#survivors,Cputime(t0);
end for;
print "SIEVE_4_20_P3_DONE",104,#survivors,survivors,Cputime(t0);
quit;
