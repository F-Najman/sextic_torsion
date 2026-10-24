def g_mn(m,n):
    N=m*m*(n//m) if False else m*n  # placeholder
cases=[(6,12),(3,12),(2,24),(2,26),(2,28),(2,30),(3,15),(3,18),(3,21),(3,24),(3,27),(4,16),(4,20),(4,24),(7,7),(9,9),(2,22)]
for m,n in cases:
    N=m*n  # X_1(m,n) ~ X_H(m*n) with H = {d = 1 mod n}, conjugating Gamma(m) cap Gamma1(n) by diag(m,1)
    H=[d for d in range(1,N) if gcd(d,N)==1 and d%n==1]
    G=GammaH(N,H)
    print((m,n), "genus", G.genus(), "index", G.projective_index())
for N in [45,51,57,63,65,27,36,49,54,32,35,42,81]:
    print(N, "X1 genus", Gamma1(N).genus())
