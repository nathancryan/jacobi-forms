class BinaryQF_with_index (BinaryQF):

    def __init__(self, m, abc_triple):
        self._m = ZZ(m)
        BinaryQF.__init__(self, abc_triple)
        if self[0] % self._m != 0:
            raise ValueError, "first coefficient must be a multiple of the index"


    def _repr_(self):
        return "Binary Quadratic Form of index %d, with coefficients %s, disc = %d" \
               % (self._m, tuple(self), self.discriminant())


    def index(self):
        return self._m


    def genus_char(self, D0):
        """
        Genus character mod $\Gamma_0(m)$ as in [S90].

        NOTE
            D0 needs to be a fundamental discriminant which is a square mod $4m$.
        """
        D = self.discriminant()
        m = self.index()
        assert (D0==1 or is_fundamental_discriminant(D0)) and  mod(D0, 4*m).is_square(), \
               '%d: must be fundamental discriminant and square mod 4*%d' %(D0, m)
        a, b, c = self
        a0 = a // m
        if gcd([a0, b, c, D0]) != 1:
            return 0
        if D % D0 != 0:
            return 0
        if not mod(D/D0, 4*m).is_square():
            return 0
        m1 = prime_to_m_part(m, self[2])
        m2 = m // m1
        n = QuadraticForm(ZZ, 2, [m1 * a0, b, m2 * c]).basiclemma(D0)
        return kronecker_symbol(D0, n)



def iterate(s, m, Dmax, Dmin = 0):
    """
    Return all binary quadratic forms $Q = (a, b, c)$ s.t.
     1.  $a (as^2 + bs + c) < 0$
     2.  $m | a$
     3.  $0 <= Dmin < disc(a,b,c) <= Dmax$
    sorted lexicographically accordingly to $(b^2-4ac, b mod 2m)$
    """
    res = {}
    p = s.numerator()
    q = s.denominator()
    sqrtDmax = floor(sqrt(Dmax))
    amax = Dmax * q**2 / 4
    for a in xsrange(m, amax, m):
        bmin = floor(-2*a*s - sqrtDmax)+1
        bmax = ceil(-2*a*s + sqrtDmax)
        for b in xsrange(bmin, bmax):
            D0 = max(Dmin, (b+2*a*s)**2)
            cmin = floor((b**2 - Dmax)/(4*a))
            cmax = ceil((b**2 - D0)/(4*a))
            for c in xsrange(cmin, cmax):
                D = b**2 - 4*a*c
                assert a%m == 0, "Bad form: %s (1)" % ((a,b,c),)
                assert a*s**2 + b*s + c < 0, "Bad form: %s" % ((a,b,c),)
                #
                key = (D, b % (2*m))
                #print key, "-->", (a,b,c)
                res[ key ] = res.get(key, 0) - 1
                #
                key = (D, -b % (2*m))
                #print key, "-->", (-a,-b,-c)
                res[ key ] = res.get(key, 0) + 1
    return res


def iterate_with_f( f, s, Delta0, r0, m, hi, lo = 0):
    """
    Return all binary quadratic forms $Q = (a, b, c)$ s.t.
     1.  $a (as^2 + bs + c) < 0$
     2.  $m | a$
     3.  $0 < lo <= disc(a,b,c)/Delta0 < hi$
    sorted accordingly to $( b^2-4ac, b mod 2m)$.

    INPUT
        f -- a function on rational polynomials of degree less or equal to 2.
        s -- a rational number
        Delta0 -- a fundamental discriminant which is a square mod $4m$
        r0 -- a solution of $r^2 \equiv Delta0 \bmod 4m$
        m -- a positive integer
        hi -- an integer
        lo -- an integer
    """
    m1 = 2*m / gcd(2*m, Delta0)
    if r0 == 0: r0 = m1
    Dmin, Dmax = lo*abs(Delta0), (hi-1)*abs(Delta0)
    x = polygen(QQ)
    res = {}
    p = s.numerator()
    q = s.denominator()
    sqrtDmax = floor(sqrt(Dmax))
    amax = Dmax * q**2 / 4
    for a in xsrange(m, amax+1, m):
        bmin = ceil(-2*a*s - sqrt(Dmax))
        bmax = floor(-2*a*s + sqrt(Dmax))
        for b in xsrange(bmin, bmax+1):
            D0 = max(Dmin, (b+2*a*s)**2)
            cmin = ceil((b**2 - Dmax)/(4*a))
            cmax = floor((b**2 - D0)/(4*a))
            for c in range(cmin, cmax+1):
                D = b**2 - 4*a*c
		# TODO / FIXME
		# we are skipping squares for the moment, because we are
		# missing some extra term
		if D.is_square():
			continue
                if D % Delta0 != 0 :
                    continue
                assert a%m == 0, "Bad form: %s (1)" % ((a,b,c),)
                assert a*s**2 + b*s + c < 0, "Bad form: %s" % ((a,b,c),)
                #
                Q = a*x*x + b*x + c
                #
                key = (D/Delta0, (b/r0) % m1)
                #print key, "-->", (a,b,c)
                f1 = f(Q)
                if f1 != 0:
                    res[ key ] = res.get(key, 0) - f1
                #
                key = (D/Delta0, (-b/r0) % m1)
                #print key, "-->", (-a,-b,-c)
                #mQ = a*x*x - b*x + c
                f2 = f(-Q)
                if f2 != 0:
                    res[ key ] = res.get(key, 0) + f2
    return res


def dict_plus_dict(d1, d2):
    return dict( (k, d1.get(k,0) + d2.get(k,0)) for k in union(d1.keys(), d2.keys()) )

def dict_by_scalar(d1, c):
    return dict( (k, d1[k]*c) for k in d1.keys() )

def dict_content(d):
    return QQ(0).content(map(QQ, d.values()))

def dict_content(d):
    return QQ(0).content(map(QQ, d.values()))

def dict_primitive(d):
    c = dict_content(d)
    if c == 0:
        c = 1
    return dict_by_scalar(d, 1/c)

def dict_kill_zeros(d):
    dict( (k, d[k]) for k in d if d[k] != 0)

#######################################


def compute_jacobi(ms, D0, r0=None, Dmax=100, Dmin = 0):
    k_jacobi = ms.parent().weight() / 2 + 1
    m_jacobi = ms.parent().level()
    if r0 == None:
        r0 = mod(D0,4*m_jacobi).sqrt().lift()
    assert (D0 - r0*r0) % (4*m_jacobi) == 0, "Bad parameters"
    d = make_dict(ms)
    res = {}
    for s in d:
        if s.is_infinity():
            continue
        P = d[s]
        def f(Q):
            char = BinaryQF_with_index(m_jacobi, (Q[2], Q[1], Q[0])).genus_char(D0)
            return char * inner_product(P, Q, k_jacobi)
        res1 = iterate_with_f(f, QQ(s), D0, r0, m_jacobi, Dmax, Dmin)
        #print s, P, "-->", res1
        res = dict_plus_dict(res, res1)
    #res = dict_kill_zeros(res)
    return res



######################################

class JacobiForm_fourier_expansion:

    def __init__(self, coeff, Dmax):
        """
        N: Should be initialized later by an instance of class
        JacobiForm_element, for now initilize by a Modular symbol
        """
        self._coeff = coeff
        self._Dmax = Dmax

    def __str__(self):
        return str(self._coeff)

    def __repr__(self):
        return repr(self._coeff)

    def coeff(self):
        return self._coeff

    def coeff_by_r(self):
        pass





