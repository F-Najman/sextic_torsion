for m,n in [(3,21),(3,15),(3,18),(3,24),(2,24),(2,26),(2,28),(2,30),(6,12),(4,16),(4,20)]:
    N=m*n
    H=[d for d in range(1,N) if gcd(d,N)==1 and d%n==1]
    G=GammaH(N,H)
    print((m,n),"cusps",G.ncusps())
for N in [45,51,57,63,65,27,36,49,35,42,40,32,34,38,39,50,33]:
    print(N,"X1 cusps",Gamma1(N).ncusps())
