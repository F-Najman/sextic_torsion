assert N eq 51 and gx eq 33 and Dimension(U) eq 30;
assert #cs eq 20 and {*Degree(c):c in cs*} eq {*1^^16,4^^4*};
assert #vs eq 54808 and &+[-Valuation(jaf,c)*Degree(c):c in cs] eq 576;
StopAfterExceptions:=0;
/* Exact local restriction kernels from DifferentialSpace, no series. */
JetFile:="data/fi51_retained_rr_kernels.m";
Write(JetFile,Sprintf("// Direct differential-space restriction certificate.\nN:=%o; As:=[]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",N,D,gx,udim,ds));
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
VO:=VectorSpace(GF(2),gx);localjets:=AssociativeArray();
for i in [1..#cs] do
    maxm:=D div ds[i];cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "C51_RR_KERNEL_START cusp=%o/%o multiplicity=%o time=%o\n",i,#cs,m,Cputime(t0);
        VV,vvmap:=DifferentialSpace(m*Divisor(cs[i]));
        W:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
        Write(JetFile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Dimension(W),gx,Eltseq(BasisMatrix(W))));
        ann:=Nullspace(Transpose(BasisMatrix(W)));assert annprev subset ann;
        extra:=[];
        for v in Basis(ann) do
            if v notin annprev then Append(~extra,v);annprev:=sub<VO|Basis(annprev) cat [v]>;end if;
        end for;
        assert #extra le ds[i];cols cat:=extra;
        while #cols lt m*ds[i] do Append(~cols,VO!0);end while;
        assert annprev eq ann;
        A:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));
        assert Nullspace(A) eq W;
        J:=BasisMatrix(U)*A;localjets[<i,m>]:=J;
        Write(JetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Nrows(J),Ncols(J),Eltseq(J)));
        printf "C51_RR_KERNEL_DONE cusp=%o/%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",i,#cs,m,Dimension(W),Rank(J),Cputime(t0);
    end for;
    Write(JetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),%o,%o,%o);",i,Nrows(A),Ncols(A),Eltseq(A)));
end for;
Write(JetFile,"print \"FI_JET_DATA_READY\";");

VO:=VectorSpace(GF(2),gx); bad:=[]; tested:=0;
for v in vs do
    if Limit gt 0 and tested ge Limit then break; end if;
    A:=ZeroMatrix(GF(2),udim,0);
    for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]); end if; end for;
    assert Ncols(A) eq D; r:=Rank(A);
    if tested lt DirectCheck or r ne D then
        Ediv:=&+[v[i]*Divisor(cs[i]):i in [1..#cs]]; assert Degree(Ediv) eq D;
        V,vm:=DifferentialSpace(Ediv);
        V_original:=sub<VO|[VO!Eltseq(vm(b) @@ om):b in Basis(V)]>;
        direct_rank:=udim-Dimension(U meet V_original); assert direct_rank eq r;
        printf "RANK0_CARTIER_DIRECT N=%o index=%o jet_rank=%o direct_rank=%o time=%o\n",N,tested+1,r,direct_rank,Cputime(t0);
    end if;
    tested+:=1;
    if r ne D then Append(~bad,v); printf "RANK0_CARTIER_EXCEPTION N=%o index=%o rank=%o vector=%o time=%o\n",N,tested,r,v,Cputime(t0); end if;
    if tested mod 1000 eq 0 then printf "RANK0_CARTIER_PROGRESS N=%o tested=%o total=%o exceptions=%o time=%o\n",N,tested,#vs,#bad,Cputime(t0); end if;
    if StopAfterExceptions gt 0 and #bad ge StopAfterExceptions then print "RANK0_CARTIER_EXCEPTION_LIMIT_REACHED"; break; end if;
end for;
printf "RANK0_CARTIER_DONE N=%o D=%o tested=%o total=%o exceptions=%o complete=%o time=%o\n",N,D,tested,#vs,#bad,tested eq #vs,Cputime(t0);
quit;
