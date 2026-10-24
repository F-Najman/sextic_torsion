/* Optional exact JOne49 certificate: if ker(A5) meet ker(T2) is zero
   modulo 2, T2 is nonzero on every rational nonzero 2-torsion point.
   Together with the published cyclic order-two 2-primary subgroup,
   this proves T2 acts as identity there without a finite-flat argument. */
SetColumns(0);
SetMemoryLimit(2*10^9);
t0:=Cputime(); N:=49;
printf "TWO49_START exact_JOne=true sign=0\n";
J:=JOne(N); assert Sign(J) eq 0;
H:=IntegralHomology(J); d:=Rank(H); assert d eq 2*Dimension(J);
printf "TWO49_HOMOLOGY dimension=%o rank=%o time=%o\n",Dimension(J),d,Cputime(t0);
function TZ(J,n)
    return ChangeRing(IntegralMatrix(HeckeOperator(J,n)),Integers());
end function;
T2Z:=TZ(J,2);
T5Z:=TZ(J,5); T25Z:=TZ(J,25);
D5Z:=ChangeRing((ChangeRing(T5Z,Rationals())^2-ChangeRing(T25Z,Rationals()))/5,Integers());
assert Abs(Determinant(D5Z)) eq 1;
A5Z:=T5Z-IdentityMatrix(Integers(),d)-5*D5Z;
assert T2Z*A5Z eq A5Z*T2Z;
K:=Nullspace(ChangeRing(A5Z,GF(2)));
Z:=K meet Nullspace(ChangeRing(T2Z,GF(2)));
printf "TWO49_KERNEL kernel_A5=%o kernel_A5_T2=%o time=%o\n",Dimension(K),Dimension(Z),Cputime(t0);
if not assigned MatrixFile then MatrixFile:="";end if;
if #MatrixFile gt 0 then
 Write(MatrixFile,Sprintf("// Exact integral JOne49 homology matrices.\nd:=%o;",d) : Overwrite:=true);
 Write(MatrixFile,Sprintf("T2Z:=Matrix(Integers(),%o,%o,%o);",d,d,Eltseq(T2Z)));
 Write(MatrixFile,Sprintf("A5Z:=Matrix(Integers(),%o,%o,%o);",d,d,Eltseq(A5Z)));
end if;
printf "TWO49_DONE zero=%o time=%o\n",Dimension(Z) eq 0,Cputime(t0);
quit;
