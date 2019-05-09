cdef genus_char(m, abc, D0):
    """
    Genus character mod $\Gamma_0(m)$ as in [S90].

    NOTE
        D0 needs to be a fundamental discriminant which is a square mod $4m$.
    """
    #assert (D0==1 or is_fundamental_discriminant(D0)) and  mod(D0, 4*m).is_square(), \
    #       '%d: must be fundamental discriminant and square mod 4*%d' %(D0, m)
    a, b, c = abc
    #D = b**2 - 4*a*c
    a0 = a // m
    #if gcd([a0, b, c, D0]) != 1:
    #    return 0
    #if D % D0 != 0:
    #    return 0
    #if not mod(D/D0, 4*m).is_square():
    #    return 0
    m2 = prime_to_m_part(m, a0)
    m1 = m // m2
    D1 = prime_to_m_part(D0, m1*a0)
    D2 = D0 // D1
    # check that D1 and D2 are discriminants, or change sign
    if (D1%4==3 or D2%4==3): D1, D2 = -D1, -D2
    # See [GKZ ; I.2, Prop 1 (P3)]
    return kronecker_symbol(D1,m1*a0) * kronecker_symbol(D2,m2*c)
    #NOTE: everything except the last part depends only on D0, a

from collections import defaultdict
from sage.arith.misc import prime_to_m_part, gcd, is_square
from sage.functions.other import sqrt
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.number_field.number_field import is_fundamental_discriminant

cdef polynomial_part(m):
    x = polygen(QQ)
    return x**m.i()

cdef make_dict(a):
    # a is a modular symbol, I make a dictionary s->polynomials where a is a formal sum of symbols polynomials*{s,infty}
    d = defaultdict(int)
    aa = a.modular_symbol_rep()# I need to do this cause a modular symbol itself is not iterable
    for coeff, symbol in aa:
        d[symbol.alpha()] += coeff * polynomial_part(symbol)
        d[symbol.beta()] -= coeff * polynomial_part(symbol)
    return d



cdef inner_product(P,Q,k):
    # P is a summand of a modular symbol corresponding to a weight 2k-2 modular form
    # Q is the polynomial of a quadratic form corresponding to a Heegner cycle
    # k is the weight of the resulting Jacobi form
    #print P, Q, k
    x = P.parent().gen(0)
    newQ = Q**(k-2)
    degree = 2*k-4 # w in the paper
    valuation = P.valuation(x)
    res = 0
    for l in range(valuation, degree+1):
        res += (-1)**l*P[l]*newQ[degree-l]/binomial(degree,l)
    return res


from sage.functions.other import ceil, floor, binomial
from sage.arith.srange import xsrange
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.arith.misc import kronecker_symbol
from sage.rings.integer_ring import ZZ

cdef iterate_with_Pk(P, long k, s, long Delta0, r0, long m, hi, lo=0, res=None):
    """
    Return all binary quadratic forms $Q = (a, b, c)$ s.t.
     1.  $a (as^2 + bs + c) < 0$
     2.  $m | a$
     3.  $0 < lo <= disc(a,b,c)/Delta0 < hi$
    sorted accordingly to $( b^2-4ac, b mod 2m)$.
    
    INPUT
        P -- polynomial for inner_product
        k -- weight for inner product
        s -- a rational number
        Delta0 -- a fundamental discriminant which is a square mod $4m$
        r0 -- a solution of $r^2 \equiv Delta0 \bmod 4m$
        m -- a positive integer
        hi -- an integer
        lo -- an integer
    """
    if res is None:
        res = defaultdict(int)
    cdef long bmin, bmax
    cdef long cmin, cmax, cmax1
    cdef long D, Delta1
    cdef long a, b, c
    cdef long b2, bmod
    cdef long gc
    ### used in precalc for genus_char
    #cdef long a0, mc, ma, D1, D2, ks
    cdef long Dmin = ceil(lo)*abs(Delta0)
    cdef long Dmax = (ceil(hi)-1)*abs(Delta0)
    cdef long p = s.numerator(), q = s.denominator()
    cdef long p2=p*p, pq=p*q, q2=q*q
    cdef long sqrtDmaxq = floor(sqrt(Dmax)*q)
    cdef long amax = Dmax * q2 // 4
    if k == 2:
        val = P[0] # may be rational
    else:
        x = polygen(QQ)
    ##assert 4*amax <= Dmax*q*q < 4*(amax+1)
    # a multiple of m and
    # 0 < a <= Dmax * q**2 / 4
    for a in range(m, amax+1, m):
        bmin = -(( 2*a*p + sqrtDmaxq) // q)
        bmax =  ((-2*a*p + sqrtDmaxq) // q)
        ##assert bmin - 1 <  ( (-2*a*p - sqrt(Dmax)*q) / q) <= bmin
        ##assert bmax     <= ( (-2*a*p + sqrt(Dmax)*q) / q) <  bmax+1
        # b integer and
        # -2*a*s - sqrt(Dmax) <= b <= -2*a*s + sqrt(Dmax)
        for b in range(bmin, bmax+1):
            #D0 = -((-(q*b + 2*a*p)**2)  //  (q*q))
            b2 = b*b
            cmin = -((Dmax - b2) //  (4*a))
            cmax = -((a*p2+b*pq) // (q2))-1
            ## this modification for using Dmin
            cmax1 =  (b2-Dmin) // (4*a)
            if cmax > cmax1 : cmax = cmax1
            ##assert cmin - 1 <  (b**2-Dmax)/(4*a) <= cmin
            ##assert cmax     <= (b**2-(b+2*a*s)**2)/(4*a)
            ##assert cmax     <= (b**2-Dmin)/(4*a)
            ##assert cmax+1 > (b**2-Dmin)/(4*a) or cmax+1 > (b**2-(b+2*a*s)**2)/(4*a)
            #
            # c integer and
            # (b**2 - Dmax) / 4a <= c < (-a*p**2 - b*p*q) / q**2
            #                   and c <= (b**2 - Dmin) / 4a
            # ALSO we need
            # 4ac - b2 = D1 * Delta0
            # 4ac - b2 = -b2 (mod 4a)
            #
            # c = (b2 + D1 * Delta0) / 4a   ### cmin = -((Dmax-b2) // (4*a))
            # b2 - 4*a*cmin = (Dmax - (Dmax-b2)%(4*a))
            #      D/Delta0 --> xx/Delta0
            # cmin = b2 - (xx - xx%Delta0)
            if cmin < cmax+1:
                bmod = b%(2*m)
                # PRECALC genus symbol -- it seems not worth it
                # a0 = a // m
                # mc = prime_to_m_part(m, a0)
                # ma = m // mc
                # D1 = prime_to_m_part(Delta0, ma*a0)
                # D2 = Delta0 // D1
                # # check that D1 and D2 are discriminants, or change sign
                # if (D1%4==3 or D2%4==3): D1, D2 = -D1, -D2
                # ks = kronecker_symbol(D1,ma*a0) * kronecker_symbol(D2,mc)
                ##### Now genus_char = ks * kronecker_symbol(D2, c)
                for c in range(cmin, cmax+1):
                    D = b2 - 4*a*c
                    if D % Delta0 != 0:
                        continue
                    Delta1 = D // Delta0
                    if Delta1 % 4 > 1:
                        continue
                    # with PRECALC genus symbol
                    #gc = ks * kronecker_symbol(D2, c)
                    gc = genus_char(m, (a,b,c), Delta0)
                    if k > 2 and gc != 0:
                        val = inner_product(P, a*x*x+b*x+c, k)
                    if gc == 1:
                        res[Delta1, bmod] -= val
                    elif gc == -1:
                        res[Delta1, bmod] += val
    return res

cdef check_D0_r0(m, D0, r0=None):
    D0_mod_4m = mod(D0, 4*m)
    if r0 is not None and D0_mod_4m == r0*r0:
        return r0
    if D0_mod_4m.is_square():
        if r0 is None:
            return D0_mod_4m.sqrt().lift()
        sug_r0 = [x for x in D0_mod_4m.sqrt(all=True) if x < 2*m]
        raise ValueError, "D0 must be r0^2 mod 4m, try r0 in %s" % sug_r0
    else:
        e = D0.sign()
        sug_D0 = [d for d in range(e,e*1000,e) if is_fundamental_discriminant(d)
                    and mod(d,4*m).is_square()]
        raise ValueError, "D0 must be a square mod 4m, try D0 in %s" % sug_D0[:5]

def compute_jacobi_fast(ms, D0, r0=None, Dmax=100, Dmin = 0):
    #if ms.ambient_module().sign() != 0:
    #    raise ValueError, "Need a modular symbol in ambient module with sign=0"
    k_jacobi = ms.parent().weight() / 2 + 1
    m_jacobi = ms.parent().level()
    r0 = check_D0_r0(m_jacobi, D0, r0)
    d = make_dict(ms)
    res = defaultdict(int)
    for s in d:
        if s.is_infinity():
            continue
        P = d[s]
        res = iterate_with_Pk(P, k_jacobi, QQ(s),
                D0, r0, m_jacobi, Dmax, Dmin, res)
    eps = (-1)**k_jacobi * D0.sign()
    m1 = 2*m_jacobi
    res_normal = defaultdict(int)
    for (D1, b), val in (res.iteritems()):
        rs = [r1.lift() for r1 in mod(D1, 4*m_jacobi).sqrt(all=True)
                          if r1 < 2*m_jacobi]
        # TODO / FIXME
        # we are skipping squares for the moment, because we are
        # missing some extra term
        if is_square(D1*D0):
            if D1 != 0:
                for r1 in rs:
                    res_normal[D1,r1] = 'N/A' #float('nan')
            continue
        for r1 in rs:
            if (r0*r1-b)%m1 == 0:
                #print r1,
                Dr1 = D1, r1
                res_normal[Dr1] = res_normal[Dr1] + val
            if (r0*r1+b)%m1 == 0:
                #print r1,
                Dr1 = D1, r1
                res_normal[Dr1] = res_normal[Dr1] - eps*val
    return dict(res_normal)
