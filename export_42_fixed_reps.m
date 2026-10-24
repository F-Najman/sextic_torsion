/* Load into the interrupted, completed-prefix C42 sieve session.
   Export one exact indexed list, with a round-trip place identity check. */
assert N eq 42 and M eq 1 and p eq 5 and D eq 6;
assert [#places[d]:d in [1..6]] eq [0,0,4,183,648,2286];
assert [#reps[d]:d in [1..6]] eq [0,0,2,33,108,392];
function C42DescribeAndCheck(x)
 dd:=Degree(x);print "C42_EXPORT_STEP",1,Cputime(t0);
 pt:=RepresentativePoint(x);print "C42_EXPORT_STEP",2,Cputime(t0);K:=GF(5^dd);
 vv:=MDM2ToCommonField(Eltseq(pt),K);print "C42_EXPORT_STEP",3,Cputime(t0);
 def:=[Integers()!a:a in Eltseq(DefiningPolynomial(K,GF(5)))];
 coords:=[[Integers()!a:a in Eltseq(v,GF(5))]:v in vv];
 L:=ext<GF(5)|Polynomial(GF(5),def)>;
 xx:=Places(Curve(X)(L)![L!cs:cs in coords]);print "C42_EXPORT_STEP",4,Cputime(t0);
 assert #xx eq 1 and xx[1] eq x;
 return [*dd,def,coords*];
end function;
C42FixedRepresentatives:=[* *];
for exportidx in [1..#reps[6]] do
 print "C42_EXPORT_INPUT",exportidx;
 Append(~C42FixedRepresentatives,C42DescribeAndCheck(reps[6][exportidx]));
end for;
assert #C42FixedRepresentatives eq 392;
Write("data/c42_fixed_representatives.m",
 "C42FixedRepresentatives := " cat Sprint(C42FixedRepresentatives) cat ";" : Overwrite:=true);
printf "C42_FIXED_REPRESENTATIVES_EXPORTED count=%o roundtrip=true time=%o\n",#C42FixedRepresentatives,Cputime(t0);
/* Only the [3,3] partition remains outside the four pure sextic shards. */
dd:=6;lambda:=[3,3];divs:=DivisorsOfType(lambda);assert #divs eq 8;
surv:=0;qlist:=[11,13,17];
for idx in [1..#divs] do
 Dv:=divs[idx];support,multiplicities:=Support(Dv);killed:=false;
 for q in qlist do
  AqD:=DivisorGroup(Curve(X))!0;
  for si in [1..#support] do
   key:=<q,support[si]>;
   if not IsDefined(HeckeCache,key) then
    PP:=Divisor(support[si]);
    HeckeCache[key]:=HeckeOp(X,q,PP)-q*DiamondOp(X,q,PP)-PP;
   end if;
   AqD+:=multiplicities[si]*HeckeCache[key];
  end for;
  if not IsPrincipal(AqD) then killed:=true;break;end if;
 end for;
 if not killed then surv+:=1;print "C42_CUBIC_PAIR_SURVIVOR",idx,DescribeData(Dv);end if;
 print "C42_CUBIC_PAIR_TEST",idx,killed,Cputime(t0);
end for;
printf "C42_CUBIC_PAIR_DONE total=8 survivors=%o complete=true time=%o\n",surv,Cputime(t0);
quit;
