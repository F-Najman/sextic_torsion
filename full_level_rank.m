/* Exact winding-submodule checks for full levels 7 and 9.
   Let V = direct sum of S_2(Gamma_0(m^2),eps) for even characters eps mod m.
   V is stable under twisting by characters modulo m: m^2 divides the level,
   and cond(eps) divides m. Thus nonvanishing for all its simple new factors
   also proves rank zero after base change to Q(zeta_m).
   See the accompanying mathematical note for this argument and the quotient. */
SetColumns(0);
SetMemoryLimit(2*10^9);
for m in [7,9] do
 D := FullDirichletGroup(m);
 good := true;
 totaldim := 0;
 for eps in GaloisConjugacyRepresentatives(D) do
  if IsOdd(eps) then continue; end if;
  ee := Extend(eps,m^2);
  S := CuspidalSubspace(ModularSymbols(ee,2,1));
  printf "CHARACTER m=%o order=%o cuspidal_dimension=%o\n",m,Order(eps),Dimension(S);
  for Si in NewformDecomposition(S) do
   T := AssociatedNewSpace(Si);
   dd := Dimension(T); ww := Dimension(WindingSubmodule(T));
   printf "SIMPLE_FACTOR m=%o level=%o dimension=%o winding_dimension=%o\n",m,Level(T),dd,ww;
   good := good and dd eq ww;
  end for;
 end for;
 printf "ALL_L_VALUES_NONZERO m=%o result=%o\n",m,good;
 assert good;
end for;
print "FULL_LEVEL_RANK_DONE";
