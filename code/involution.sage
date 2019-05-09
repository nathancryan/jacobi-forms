

### Example -- k=5 -- m=13 /// poly of degree 6
# this is a rational eigenform

AS = ModularSymbols(13,8)

MS = AS.minus_submodule().decomposition()[0]
assert MS.dimension() == 1

ms = MS.gen(0)

msg = sum([(val if is_even(s.i) else -val)*AS(s)   for val,s in ms.manin_symbol_rep()])

print [MS.hecke_operator(p)[0,0] for p in prime_range(20)]


def Slift(m, k, coeff, D0, r0, prec):
    # fix leading coefficient ??
    # in general, a special value of dirichlet L-series of kronecker(D0, *)
    # times self.coeff(0,0) (see ...)
    # for D0=1, k=2, also add peterson of self with Tr(r0)
    # (see Skoruppa "Explicit..." p.514)
    coeff = defaultdict(int, coeff)
    return [ 
        sum((kronecker(D0, a) * a**(k-2) * coeff[D0*l*l/a/a, r0*l/a % (2*m)]
            for a in divisors(l)))
    for l in range(1,prec) ]
