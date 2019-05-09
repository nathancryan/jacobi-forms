###
#

# Start by attaching the .pyx file
#
# attach("compute_jacobi.pyx")
#
# Then run the examples below

####
# WEIGHT 2

# level = 11 , weight = 2  (11a)
#
# use with positive discriminants (skew holomorphic)
#
ms11 = ModularSymbols(Gamma0(11), 2).minus_submodule().cuspidal_submodule().gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms11, 1, 1, 100)
# compute_jacobi_fast(ms11, 5, 7, 200)
# compute_jacobi_fast(ms11, 12, 10, 200)
# compute_jacobi_fast(ms11, 33, 11, 200) ### IS ZERO, WHY?



# level = 37 , weight = 2 , sign -1 (37a)
#
# use with negative discriminants (holomorphic)
#
ms37a = kernel(ModularSymbols(Gamma0(37), 2).plus_submodule().cuspidal_submodule().atkin_lehner_operator() - 1).gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms37a, 1, 1, 100)  ### IS ZERO, WHY?
# compute_jacobi_fast(ms37a, -3, 21, 200)
# compute_jacobi_fast(ms37a, -4, 12, 200)


# level = 37 , weight = 2 , sign +1 (37b)
#
# use with positive discriminants (skew holomorphic)
#
ms37b = kernel(ModularSymbols(Gamma0(37), 2).minus_submodule().cuspidal_submodule().atkin_lehner_operator() + 1).gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms37b, 1, 1, 100)


# level = 277, weight 2, sign +1
#
# use with negative discriminants (holomorphic)
# compare with the paramodular form that is a nonlift
# 
ms277_basis = kernel(ModularSymbols(Gamma0(277), 2).plus_submodule().cuspidal_submodule().new_submodule().atkin_lehner_operator() - 1).gens()
#
# CALCULATIONS:
#
# jmf277_basis = [compute_jacobi_fast(phi,-7,47,500) for phi in ms277_basis]
#
# from Poor Yuen we know the coefficients of the nonlift. After using the Twin operator 
# realizing that they use [2a,b,2c] notation, these are the coefficients we know
#
# known_coeffs_277 = {(-3,233):-3, (-4,120):-2, (-67,-137):-8, (-64,-74):-1, (-88,-110):-1, (-187,-87):4, (-264,-104):-10, (-276,-276):9, (-304,-132):-6, (-340,-136):7, (-448,-178):8, (-479,-181):8, (-491,-127):4}
#
# construct the matrix of coefficients of the basis elements
#
#

# level = 554 = 2*277, weight 2, sign +1
#
# use with negative discriminants (holomorphic)
# compare with the paramodular form that is a nonlift
#
ms554_basis = kernel(ModularSymbols(Gamma0(2*277), 2).plus_submodule().cuspidal_submodule().new_submodule().atkin_lehner_operator() - 1).gens()
#
# CALCULATIONS:
# 
# jmf554_basis = [compute_jacobi_fast(phi,-7,47,400) for phi in ms554_basis]
#
#from Poor Yuen we know the coefficients of the nonlift. After using the Twin operator 
# realizing that they use [2a,b,2c] notation, these are the coefficients we know
#
# known_coeffs_554 = {(-7,601):-1, (-12,466):6, (-16,240):6, (-116,-182):7, (-256,-148):-2, (-268,-274): 16, (-351,-155): 10, (-352,-220): 3, (-399,-175):-4} 


# level = 15, weight=2, sign=+1
# composite level, skew holomorphic
# compare with Ariel and Gonzalo's paper in Exp Math
ms15 = ModularSymbols(Gamma0(15), 2).minus_submodule().cuspidal_submodule().gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms15,61,1,50)




####
## HIGHER WEIGHT

# level = 3 , weight = 16 (jacobi weight = 9), sign -1
#
# use with positive discriminants (skew holomorphic)
#
ms3k16a = (ModularSymbols(Gamma0(3), 16).plus_submodule().new_submodule().cuspidal_subspace().atkin_lehner_operator()-1).kernel().gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms3k16a, 1, 1, 100)
# compute_jacobi_fast(ms3k16a, 13, 1, 300)


# level = 3, weight = 16 (jacobi weight = 9), sign +1
#
# use with negative discriminants (holomorphic)
#
ms3k16b = (ModularSymbols(Gamma0(3), 16).minus_submodule().new_submodule().cuspidal_subspace().atkin_lehner_operator()+1).kernel().gen(0)
#
# EXAMPLES:
#
# compute_jacobi_fast(ms3k16b, -3, 3, 300)   ## SHOULD BE 0


# level = 1, weight = 18 (jacobi weight = 10), sign -1
#
# use with negative discriminants (holomorphic)
# compare to the coefficients of the Maass lift Siegel modular form
#
ms1k18a = (ModularSymbols(1, 18).plus_submodule().new_submodule().cuspidal_subspace().atkin_lehner_operator()-1).kernel().gen(0)
#
# EXAMPLES: 
#
# compute_jacobi_fast(ms1k18a*5/224, -4, 0, 100) # Agrees with LMFDB
# compute_jacobi_fast(ms1k18a*5/448, -3, 1, 300) # Agrees with LMFDB



