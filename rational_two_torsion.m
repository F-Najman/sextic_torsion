/* Rigorous upper bound for J1(N)(Q)[2] via integral Betti homology.
   Uses JOne, NOT Js (the latter is merely isogenous and is unsafe at2).
   All homology has sign0.  For odd good q, Aq=Tq-1-q<q> kills rational
   2-torsion by Eichler--Shimura and injectivity of reduction away from2.
   Recover diamond matrices by Tq^2-T(q^2)=q<q> over Z before reducing2.
   Shell timeout bounds the intended N57 exploratory job to three minutes.
   A success certifies that diamond2 acts trivially on this upper bound.
*/
SetColumns(0);
SetMemoryLimit(2*10^9);
if not assigned N then N:=57; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned Qs then Qs:=[5,7,11]; elif Type(Qs) eq MonStgElt then Qs:=[StringToInteger(s):s in Split(Qs,",")]; end if;
if not assigned MatrixFile then MatrixFile:=""; end if;
t0:=Cputime();
printf "RATIONAL_TWO_START N=%o Qs=%o sign=0 exact_JOne=true\n",N,Qs;
J:=JOne(N);
assert Sign(J) eq 0;
printf "RATIONAL_TWO_JONE N=%o dimension=%o time=%o\n",N,Dimension(J),Cputime(t0);
H:=IntegralHomology(J);
d:=Rank(H); assert d eq 2*Dimension(J);
printf "RATIONAL_TWO_HOMOLOGY N=%o rank=%o time=%o\n",N,d,Cputime(t0);
function TZ(J,n)
    // IntegralMatrix expresses the map on the SAME integral homology lattice.
    return ChangeRing(IntegralMatrix(HeckeOperator(J,n)),Integers());
end function;
T2:=TZ(J,2); T4:=TZ(J,4);
D2Q:=(ChangeRing(T2,Rationals())^2-ChangeRing(T4,Rationals()))/2;
D2Z:=ChangeRing(D2Q,Integers());
assert Abs(Determinant(D2Z)) eq 1;
if N eq 57 then assert D2Z^9 eq IdentityMatrix(Integers(),d); end if;
D2:=ChangeRing(D2Z,GF(2)); I:=IdentityMatrix(GF(2),d);
K:=VectorSpace(GF(2),d);
printf "RATIONAL_TWO_DIAMOND N=%o dimension=%o time=%o\n",N,d,Cputime(t0);
if #MatrixFile gt 0 then
    Write(MatrixFile,Sprintf("// Exact JOne(%o) integral homology, sign0. Fresh export path.\nN:=%o; d:=%o;",N,N,d) : Overwrite:=true);
    Write(MatrixFile,Sprintf("D2Z:=Matrix(Integers(),%o,%o,%o);",d,d,Eltseq(D2Z)));
end if;
passed:=false;
for q in Qs do
    assert IsPrime(q) and IsOdd(q) and N mod q ne 0;
    printf "RATIONAL_TWO_HECKE_START N=%o q=%o time=%o\n",N,q,Cputime(t0);
    Tq:=TZ(J,q); Tqq:=TZ(J,q^2);
    DqQ:=(ChangeRing(Tq,Rationals())^2-ChangeRing(Tqq,Rationals()))/q;
    DqZ:=ChangeRing(DqQ,Integers());
    assert Tq*D2Z eq D2Z*Tq;
    AqZ:=Tq-IdentityMatrix(Integers(),d)-q*DqZ;
    Aq:=ChangeRing(AqZ,GF(2));
    K:=K meet Nullspace(Aq);
    assert &and[v*D2 in K:v in Basis(K)];
    action_rank:=#Basis(K) eq 0 select 0 else Rank(Matrix(GF(2),[Eltseq(v*(D2-I)):v in Basis(K)]));
    printf "RATIONAL_TWO_KERNEL N=%o q=%o kernel_dimension=%o diamond_difference_rank=%o time=%o\n",N,q,Dimension(K),action_rank,Cputime(t0);
    if #MatrixFile gt 0 then
        Write(MatrixFile,Sprintf("Aq%oZ:=Matrix(Integers(),%o,%o,%o);",q,d,d,Eltseq(AqZ)));
    end if;
    if action_rank eq 0 then passed:=true; break; end if;
end for;
printf "RATIONAL_TWO_DONE N=%o kernel_dimension=%o diamond_trivial=%o time=%o\n",N,Dimension(K),passed,Cputime(t0);
quit;
