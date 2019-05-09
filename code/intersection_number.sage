## k = 3 #weight
## m = 11 #index
## B = BinaryQF([1,2,3])
## BB = B.polynomial()
## R.<X,Y> = PolynomialRing(ZZ)
## Q = R(BB) # got to be a better way of doing this
## period = ModularSymbols(Gamma0(m),2*k-2).gens()[0]
#P = period.modular_symbol_rep()[0][1].polynomial_part() # got to be a better way of doing this

def inner_product(P,Q,k):
    # P is a summand of a modular symbol corresponding to a weight 2k-2 modular form
    # Q is the polynomial of a quadratic form corresponding to a Heegner cycle
    # k is the weight of the resulting Jacobi form
    #print P, Q, k
    if P.is_zero():
        return 0
    x = P.parent().0
    newQ = Q^(k-2)
    degree = 2*k-4 # w in the paper
    valuation = P.valuation(x)
    res = sum((-1)^l*P[l]*newQ[degree-l]/binomial(degree,l) for l in range(valuation, degree+1))
    return res

from collections import defaultdict

def polynomial_part(m):
    x = polygen(QQ)
    return x**m.i()

def make_dict(a):
    # a is a modular symbol, I make a dictionary s->polynomials where a is a formal sum of symbols polynomials*{s,infty}
    d = defaultdict(int)
    aa = a.modular_symbol_rep()# I need to do this cause a modular symbol itself is not iterable
    for coeff, symbol in aa:
        d[symbol.alpha()] += coeff * polynomial_part(symbol)
        d[symbol.beta()] -= coeff * polynomial_part(symbol)
    return d

