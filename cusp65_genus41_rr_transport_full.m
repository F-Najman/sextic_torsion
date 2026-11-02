/* Direct genus41 H3 quotient of X1(65), AS cover of genus21 fiber product (Lemma 9.4).
   Exact j13 and twist parameter certified in setup13_twist_parameter_certified.m.
   Run: magma -n -b cusp65_genus41_rr_transport_full.m (see jobs.txt and README.md). */
SetColumns(0);SetMemoryLimit(24*10^9);t0:=Cputime();
F0<z13>:=FunctionField(GF(2));
j13:=(z13^16+z13^4+1)/(z13^28+z13^25+z13^24+z13^21+z13^20+z13^17+z13^16+z13^13);
A13:=(z13^7+z13^6+z13^2+z13)/(z13^4+z13+1);
assert Evaluate(A13,z13+1)+A13 eq z13^2+z13;
R<t>:=PolynomialRing(F0);
f:=t^12+t^8+t^4+1+j13*(t^7+t^6+t^5);
assert IsIrreducible(f);
K<t5>:=FunctionField(f);assert Degree(K) eq 12 and Genus(K) eq 21;
printf "C65_G41_BASE genus=21 degree=12 time=%o\n",Cputime(t0);
RR<v>:=PolynomialRing(K);
A5:=t5^2/(t5+1)^3; A:=A5+K!A13;
aspoly:=v^2+v+A;assert IsIrreducible(aspoly);
AF<w>:=FunctionField(aspoly);assert Degree(AF) eq 2;
gx:=Genus(AF);assert gx eq 41;
printf "C65_G41_MODEL genus=41 degree=24 time=%o\n",Cputime(t0);
mbT:=hom<F0->AF|AF!z13>;
mbU:=hom<F0->AF|AF!z13+1>;
mT:=hom<K->AF|mbT,1/(AF!t5)>;
mU:=hom<K->AF|mbU,AF!t5>;
tauT:=hom<AF->AF|mT,w+(AF!t5)/(AF!t5+1)>;
tauU:=hom<AF->AF|mU,w+AF!z13>;
for a in [AF!z13,AF!t5,w] do
    assert tauT(tauT(tauT(tauT(a)))) eq a;
    assert tauU(tauU(tauU(tauU(a)))) eq a;
    assert tauT(tauT(a)) eq tauU(tauU(a));
    assert tauT(tauU(a)) eq tauU(tauT(a));
end for;
assert tauT(tauT(w)) eq w+1;
printf "C65_G41_AUTOMORPHISMS order=8 time=%o\n",Cputime(t0);
Omega,om:=SpaceOfHolomorphicDifferentials(AF);
B:=[om(v):v in Basis(Omega)];assert #B eq 41;
etaAF:=Differential(AF!z13);
MT:=Matrix(GF(2),[Eltseq((tauT(b/etaAF)*etaAF) @@ om):b in B]);
MU:=Matrix(GF(2),[Eltseq((tauU(b/etaAF)*etaAF) @@ om):b in B]);
I:=IdentityMatrix(GF(2),gx);
assert MT^4 eq I and MU^4 eq I and MT^2 eq MU^2 and MT*MU eq MU*MT;
U0:=RowSpace(MT-I)+RowSpace(MU-I);
CM:=Matrix(GF(2),[Eltseq(Cartier(b) @@ om):b in B]);
assert CM*MT eq MT*CM and CM*MU eq MU*CM;
U:=RowSpace(BasisMatrix(U0)*(CM-I));udim:=Dimension(U);
printf "C65_G41_FORMS diamond_image=%o cartier_image=%o time=%o\n",Dimension(U0),udim,Cputime(t0);
print "C65_G41_SETUP_READY";

/* TTY continuation after exact genus41/diamond/Cartier setup. */
N:=65;D:=6;Limit:=0;DirectCheck:=1;StopAfterExceptions:=0;
assert gx eq 41 and Dimension(U) eq udim;
JetFile:="fi65_genus41_cartier_jets.m";
jaf:=AF!j13;
basecs:=Poles(j13); cs_all:=[];
for p in basecs do
    cs_all cat:=Decomposition(AF,p);
end for;
cs_all:=Setseq(Seqset(cs_all));
assert &and[Valuation(jaf,c) lt 0:c in cs_all];
pole_degree:=&+[-Valuation(jaf,c)*Degree(c):c in cs_all];
assert pole_degree eq 672 and {*Degree(c):c in cs_all*} eq {*1^^8,2^^4,4^^4*};
printf "C65_G41_CUSPS pole_degree=%o signature=%o time=%o\n",pole_degree,{*Degree(c):c in cs_all*},Cputime(t0);
cs:=[c:c in cs_all|Degree(c) le D];ds:=[Degree(c):c in cs];
steps:=[];
for c in cs do
    if Degree(c) eq 1 then st:=1;
    elif Degree(c) eq 2 then st:=3;
    elif Valuation(AF!z13,c) gt 0 or Valuation(AF!z13+1,c) gt 0 then st:=1;
    else st:=3;end if;
    Append(~steps,st);
end for;
assert {*<ds[i],steps[i]>:i in [1..#cs]*} eq {*<1,1>^^8,<2,3>^^4,<4,1>^^2,<4,3>^^2*};
printf "C65_G41_CUSP_STEPS degrees=%o steps=%o time=%o\n",ds,steps,Cputime(t0);

function WeightedVectors(ds,d)
    if #ds eq 0 then return d eq 0 select [[]] else []; end if;
    result:=[];
    for k in [0..d div ds[1]] do
        for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v); end for;
    end for;
    return result;
end function;
rawvs:=WeightedVectors([ds[i]*steps[i]:i in [1..#cs]],D);
vs:=[[steps[i]*v[i]:i in [1..#cs]]:v in rawvs];
printf "RANK0_CARTIER_COUNT N=%o total=%o time=%o\n",N,#vs,Cputime(t0);
assert #vs eq 1792;

/* Exact restriction matrices recovered from differential-space kernels,
   then transported under the verified modular diamond action. */
JetFile:="data/fi65_genus41_rr_transport.m";
function CuspPerm(h,places)
    out:=[];
    for i in [1..#places] do
        a,b:=TwoGenerators(places[i]); aa:=h(AF!a);bb:=h(AF!b);
        js:=[j:j in [1..#places]|Degree(places[j]) eq Degree(places[i]) and Valuation(aa,places[j]) gt 0 and Valuation(bb,places[j]) gt 0];
        assert #js eq 1;Append(~out,js[1]);
    end for;
    return out;
end function;
pT:=CuspPerm(tauT,cs);pU:=CuspPerm(tauU,cs);
perms_transport:=[pT,pU]; invMs:=[MT^-1,MU^-1];
S:=SymmetricGroup(#cs);G:=sub<S|S!pT,S!pU>;
assert #G eq 8;
for p in perms_transport do
    assert &and[ds[i] eq ds[p[i]] and steps[i] eq steps[p[i]]:i in [1..#cs]];
end for;
orbs:=Orbits(G);
assert {*<#o,ds[Min(o)],steps[Min(o)]>:o in orbs*} eq {*<8,1,1>,<4,2,3>,<2,4,1>,<2,4,3>*};
printf "C65_RR_TRANSPORT_PERMS T=%o U=%o orbits=%o time=%o\n",pT,pU,[Sort(Setseq(o)):o in orbs],Cputime(t0);
VO:=VectorSpace(GF(2),gx);
basisjets:=AssociativeArray();kernelspaces:=AssociativeArray();
transforms:=AssociativeArray();origins:=AssociativeArray();
for orb in orbs do
    rep:=Min(orb);maxm:=D div ds[rep]; cols:=[];annprev:=sub<VO|>;
    for m in [1..maxm] do
        printf "C65_RR_KERNEL_START cusp=%o multiplicity=%o time=%o\n",rep,m,Cputime(t0);
        VV,vvmap:=DifferentialSpace(m*Divisor(cs[rep]));
        W:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
        kernelspaces[<rep,m>]:=BasisMatrix(W);
        ann:=Nullspace(Transpose(BasisMatrix(W)));
        assert annprev subset ann;
        extra:=[];
        for v in Basis(ann) do
            if v notin annprev then Append(~extra,v);annprev:=sub<VO|Basis(annprev) cat [v]>;end if;
        end for;
        assert #extra le ds[rep];cols cat:=extra;
        while #cols lt m*ds[rep] do Append(~cols,VO!0);end while;
        assert annprev eq ann;
        printf "C65_RR_KERNEL_DONE cusp=%o multiplicity=%o kernel_dim=%o rank_on_U=%o time=%o\n",rep,m,Dimension(W),udim-Dimension(U meet W),Cputime(t0);
    end for;
    AP:=Transpose(Matrix(GF(2),[Eltseq(v):v in cols]));assert Nrows(AP) eq gx and Ncols(AP) eq maxm*ds[rep];
    for m in [1..maxm] do
        assert Nullspace(Submatrix(AP,1,1,gx,m*ds[rep])) eq RowSpace(kernelspaces[<rep,m>]);
    end for;
    transforms[rep]:=IdentityMatrix(GF(2),gx);origins[rep]:=rep;q:=[rep];
    while #q gt 0 do
        i:=q[1];Remove(~q,1);
        for k in [1,2] do
            j:=perms_transport[k][i];
            if not IsDefined(transforms,j) then
                transforms[j]:=invMs[k]*transforms[i];origins[j]:=rep;Append(~q,j);
            end if;
        end for;
    end while;
    for i in orb do
        assert origins[i] eq rep;
        basisjets[i]:=transforms[i]*AP;
    end for;
    checki:=Min(orb diff {rep});
    VV,vvmap:=DifferentialSpace(maxm*Divisor(cs[checki]));
    Wcheck:=sub<VO|[VO!Eltseq(vvmap(b) @@ om):b in Basis(VV)]>;
    assert Nullspace(basisjets[checki]) eq Wcheck;
    printf "C65_RR_TRANSPORT_DIRECT source=%o target=%o multiplicity=%o verified=true time=%o\n",rep,checki,maxm,Cputime(t0);
end for;
assert #Keys(basisjets) eq #cs;
Write(JetFile,Sprintf("// Differential-space kernel and diamond-transport certificate.\nN:=%o; As:=[\"tauT\",\"tauU\"]; D:=%o; gx:=%o; o:=1; udim:=%o; JetFactor:=1; ds:=%o; perms:=[];\nlocaljets:=AssociativeArray(); basisjets:=AssociativeArray(); kernelspaces:=AssociativeArray(); TransportMatrices:=AssociativeArray();",N,D,gx,udim,ds) : Overwrite:=true);
Write(JetFile,Sprintf("CoefficientSteps:=%o; TransportPermutations:=%o; TransportOrigins:=%o;",steps,perms_transport,[origins[i]:i in [1..#cs]]));
Write(JetFile,Sprintf("DifferenceImageBasis:=Matrix(GF(2),%o,%o,%o);",Dimension(U0),gx,Eltseq(BasisMatrix(U0))));
Write(JetFile,Sprintf("MT:=Matrix(GF(2),%o,%o,%o); MU:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(MT),gx,gx,Eltseq(MU)));
Write(JetFile,Sprintf("CartierMatrix:=Matrix(GF(2),%o,%o,%o); ImageBasis:=Matrix(GF(2),%o,%o,%o);",gx,gx,Eltseq(CM),udim,gx,Eltseq(BasisMatrix(U))));
for key in Sort(Setseq(Keys(kernelspaces))) do
    W:=kernelspaces[key];Write(JetFile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",key[1],key[2],Nrows(W),Ncols(W),Eltseq(W)));
end for;
localjets:=AssociativeArray();
for i in [1..#cs] do
    A:=basisjets[i];L:=transforms[i];
    Write(JetFile,Sprintf("TransportMatrices[%o]:=Matrix(GF(2),%o,%o,%o);",i,gx,gx,Eltseq(L)));
    Write(JetFile,Sprintf("basisjets[%o]:=Matrix(GF(2),%o,%o,%o);",i,Nrows(A),Ncols(A),Eltseq(A)));
    for m in [1..D div ds[i]] do
        J:=BasisMatrix(U)*Submatrix(A,1,1,gx,m*ds[i]);localjets[<i,m>]:=J;
        Write(JetFile,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Nrows(J),Ncols(J),Eltseq(J)));
    end for;
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
