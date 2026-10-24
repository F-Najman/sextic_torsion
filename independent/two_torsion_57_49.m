// Independent check of Remark 6.4 and Lemma 9.2, written separately from rational_two_torsion.m and
// two_torsion_49_T2.m.  Run from the repository root: magma -n -b independent/two_torsion_57_49.m
// Independent check of the 2-torsion lemmas for J1(57) and J1(49),
// using integral Hecke matrices on H_1(J1(N),Z) and only identities that are
// valid mod 2 without division:  q<q> = T_q^2 - T_{q^2}.
SetColumns(0);
SetMemoryLimit(12*10^9);
t0 := Cputime();
function M2(J, n)
  A := IntegralMatrix(HeckeOperator(J, n));
  return MatrixAlgebra(GF(2), NumberOfRows(A)) ! Eltseq(A);
end function;
// ---- J1(57) ----
J := JOne(57);
printf "J1(57) dim %o\n", Dimension(J);
T5 := M2(J,5); T25 := M2(J,25);
A5 := T5 - 1 - (T5^2 - T25);          // T5 - 1 - 5<5>, mod 2
K := Kernel(A5);
printf "57: dim ker A5 mod 2 = %o\n", Dimension(K);
// <2> = <59>;  59<59> = T59^2 - T_{59^2},  59 = 1 mod 2
T59 := M2(J,59); T3481 := M2(J,3481);
D2 := T59^2 - T3481;
printf "57: (<2>-1) vanishes on ker A5: %o   (%o s)\n", forall{v : v in Basis(K) | v*(D2-1) eq 0}, Cputime(t0);
// sanity: <2>^9 = 1 mod 2
printf "57: <2>^9 = 1 mod 2: %o\n", D2^9 eq 1;
Kt := Kernel(Transpose(A5));
printf "57 (transposed convention): dim %o, (<2>-1) vanishes: %o\n", Dimension(Kt), forall{v : v in Basis(Kt) | v*(Transpose(D2)-1) eq 0};
assert Dimension(K) eq 2 and forall{v : v in Basis(K) | v*(D2-1) eq 0};
assert Dimension(Kt) eq 2 and forall{v : v in Basis(Kt) | v*(Transpose(D2)-1) eq 0};
// ---- J1(49) ----
J := JOne(49);
printf "J1(49) dim %o\n", Dimension(J);
T5 := M2(J,5); T25 := M2(J,25); T2 := M2(J,2);
A5 := T5 - 1 - (T5^2 - T25);
K := Kernel(A5);
KT2 := Kernel(T2);
Kt := Kernel(Transpose(A5)); KT2t := Kernel(Transpose(T2));
printf "49 (transposed): dim ker A5 = %o, meet ker T2 = %o\n", Dimension(Kt), Dimension(Kt meet KT2t);
printf "49: dim ker A5 = %o, dim (ker A5 meet ker T2) = %o  (%o s)\n", Dimension(K), Dimension(K meet KT2), Cputime(t0);
assert Dimension(K) eq 2 and Dimension(K meet KT2) eq 0 and Dimension(Kt) eq 2 and Dimension(Kt meet KT2t) eq 0;
print "TWO_TORSION_CHECK_DONE";
quit;
