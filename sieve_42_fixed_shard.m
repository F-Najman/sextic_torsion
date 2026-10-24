/* One quarter of the SAME exact 392 sextic representatives.
   The master list has an exact original-place reconstruction audit. */
SetColumns(0);
if Type(Shard) eq MonStgElt then Shard:=StringToInteger(Shard);end if;
assert Shard in [0..3];SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");
Attach("fast_hecke.m");
Attach("fast_hecke_11.m");
Attach("fast_hecke_cyclic_commonzero.m");
X:=MDX1(42,GF(5));assert Genus(X) eq 25;
load "data/c42_fixed_representatives.m";
assert #C42FixedRepresentatives eq 392;
indices:=[i:i in [1..392]|i mod 4 eq Shard];assert #indices eq 98;
printf "C42_FIXED_SHARD_SETUP shard=%o/4 count=98 total=392 genus=25 p=5 qlist=[11,13,17] time=%o\n",Shard,Cputime(t0);
function Reconstruct(row)
 dd,def,coords:=Explode(row);assert dd eq 6;
 K:=ext<GF(5)|Polynomial(GF(5),def)>;assert #K eq 5^6;
 pp:=Places(Curve(X)(K)![K!cs:cs in coords]);
 assert #pp eq 1 and Degree(pp[1]) eq 6 and not IsCusp(X,pp[1]);
 return pp[1];
end function;
function DiamondPlace(x,q)
 E:=EllipticCurve(X,x);P:=LevelStructure(X,x)`P;
 y:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=q*P>);
 assert Degree(y) eq 6;return y;
end function;
tested:=[];survivors:=[];
for index in indices do
 x:=Reconstruct(C42FixedRepresentatives[index]);killed:=false;
 for q in [11,13,17] do
  A:=HeckeOperatorCommonZeros(X,q,x)-q*Divisor(DiamondPlace(x,q))-Divisor(x);
  assert Degree(A) eq 0;
  if not IsPrincipal(A) then
   printf "C42_FIXED_CASE shard=%o index=%o nonprincipal_q=%o time=%o\n",Shard,index,q,Cputime(t0);
   killed:=true;break;
  end if;
 end for;
 Append(~tested,index);
 if not killed then
  Append(~survivors,index);
  printf "C42_FIXED_SURVIVOR shard=%o index=%o time=%o\n",Shard,index,Cputime(t0);
 end if;
end for;
assert tested eq indices;
printf "C42_FIXED_SHARD_DONE shard=%o/4 tested=%o survivors=%o indices=%o complete=true time=%o\n",Shard,#tested,#survivors,survivors,Cputime(t0);
quit;
