#############################################################
#                                                           #
#           Constructing the quadratic form                 #
#                                                           #
#############################################################


def trace_zero_subspace(f):
    r"""
        Return a basis for the subspace of elements of trace zero inside a maximal order of a given number field.

        INPUT:

        - ``f`` a polynomial that defines a number field K

        EXAMPLES::

	sage: f = x^4 - 2*x^3 - 9*x^2 + 5*x + 16
	sage: trace_zero_subspace(f)
	[-16*a^3 + 34*a^2 + 1, -8*a^3 + 17*a^2 + a, -22*a^3 + 47*a^2]

        

    """
    K.<a> = NumberField(f) # construct the number field
    O_K = K.maximal_order()
    bas = O_K.basis() # find a basis
    dim = len(bas)
    M = Matrix(ZZ,dim,dim,[[x.trace() for x in bas] for _ in range(dim)]) # calculate the matrix representation of the trace map in the basis
    kerM = (M.transpose()).kernel() # find a basis for subspace of elements of trace zero

    return [ sum([v[i]*bas[i] for i in range(dim)]) for v in kerM.basis()]

def trace_zero_matrix(f):
    r"""
        Return a matrix for the subspace of elements of trace zero inside a maximal order of a given number field.

        INPUT:

        - ``f`` a polynomial that defines a number field K

        EXAMPLES::



        sage: f = x^4 - 2*x^3 - 9*x^2 + 5*x + 16
        sage: trace_zero_matrix(f)

        [191012  93618 261632]
        [ 93618  45886 128229]
        [261632 128229 358362]
    """
    bas = trace_zero_subspace(f) # find a basis for subspace of elements of trace zero
    ll = len(bas)
    MF = matrix(ZZ,ll,ll,[[(bas[i]*bas[j]).trace() for i in range(ll)] for j in range(ll)])  
    return MF

def compute_qf(f):
    r"""
        Return the quadratic form related to the matrix of the trace zero subspace inside a maximal order of a given number field.

        INPUT:

        - ``f`` a polynomial that defines a number field K

        EXAMPLES::

        sage: f = x^4 - 2*x^3 - 9*x^2 + 5*x + 16
        sage: compute_qf(f)

        Quadratic form in 3 variables over Integer Ring with coefficients: 
        [ 95506 93618 261632 ]
        [ * 22943 128229 ]
        [ * * 179181 ]

        sage: f = x^5 - 15*x^3 - 12*x^2 + 9*x + 4
        sage: compute_qf(f)

        Quadratic form in 4 variables over Integer Ring with coefficients: 
        [ 512674 5355 282974 1484437 ]
        [ * 15 1477 7762 ]
        [ * * 39048 409674 ]
        [ * * * 1074630 ]

        sage: f = x^6 - 10*x^4 - 11*x^3 + 6*x^2 + 8*x + 1
        sage: compute_qf(f)

        Quadratic form in 5 variables over Integer Ring with coefficients: 
        [ 6285671 -15073 16474954 23684498 38105438 ]
        [ * 10 -19753 -28393 -45686 ]
        [ * * 10795354 31038941 49937831 ]
        [ * * * 22310898 71790958 ]
        [ * * * * 57751400 ] 
    """
   
    M = trace_zero_matrix(f)
    return QuadraticForm(M)


########################################################
#                                                      #
#       To process files downloaded from LMFDB         #
#                                                      # 
########################################################

def candidate_fields(file,deg=4):
    r"""
	Return a dictionary, indexed by discriminants, of number fields of a 
        fixed degree for which there is more than one isomorphism class of 
        number fields.  We also require that the discriminant be fundamental 
        and relatively prime to the degree.

        INPUT:

        - ``file``- the path to a .sage files in the format produced by the global fields database on the LMFDB
        - ``deg`` - the degree of the numberfield

        EXAMPLES::
    
        sage: fields = candidate_fields("fields-deg-4-small.sage",4)
        sage: min(fields.keys())
        35537
        sage: max(fields.keys())
        3977033
        sage: len(fields.keys())
        403
    """    

    attach(file)
    good_fields = []
    for l in data:
        if (gcd(l[Integer(1)],Integer(deg))==Integer(1) and is_fundamental_discriminant(l[Integer(1)]) and l[1]>0):
            good_fields.append(l)
    temp = {}
    for tup in good_fields:
    	temp[tup[1]] = []
    for tup in good_fields:
    	temp[tup[1]].append(tup)
    ret = {}
    for dd in temp:
    	if len(temp[dd]) > 1:
	   ret[dd] = temp[dd]
    return ret

def candidate_fields_malle(file,deg=4):
    r"""
	Return a dictionary, indexed by discriminants, of number fields of a 
        fixed degree for which there is more than one isomorphism class of 
        number fields.  We also require that the discriminant be fundamental 
        and relatively prime to the degree.  These files are created in MAGMA by	Malle

        INPUT:

        - ``file``- the path to a .sage file based on files sent to Guillermo by Malle
        - ``deg`` - the degree of the numberfield

        EXAMPLES::
    
        sage: fields = candidate_fields("S4fundmult_1.sage",4)
        sage: min(fields.keys())
        
        sage: max(fields.keys())
        
        sage: len(fields.keys())
        
    """    

    attach(file)
    good_fields = []
    for l in data:
        if (gcd(l[Integer(0)],Integer(deg))==Integer(1) and is_fundamental_discriminant(l[Integer(0)]) and l[0]>0):
            good_fields.append(l)
    ret = {}
    for tup in good_fields:
     	ret[tup[0]] = []
    for dd in good_fields:
        dd[0],dd[1] = dd[1],dd[0]
        ret[dd[1]].append(dd)
    return ret



PREC = 10000 # I set the precision of the q-expansion manually in case I need to compute more coefficients


###############################################################
#                                                             #
#  Find primes represented by forms of the same discriminant  #
#                                                             #
###############################################################

def find_primes_represented(qexp):
    r"""
	Finds the list of all primes represented by a quadratic form Q

        INPUT:
        -  ``qq`` - the q-expansion of the theta series associated to the quadratic form

        EXAMPLES::
	sage: f = x^6 - 10*x^4 - 11*x^3 + 6*x^2 + 8*x + 1
	sage: qf = compute_qf(f)
	sage: theta_series = qf.theta_series(100)
	sage: find_primes_represented(theta_series)
	[37, 41, 43, 61, 67, 73, 79, 83, 89, 97]
    """

    primes = []
    for i in range(qexp.prec()):
        if (is_prime(i) and qexp[i] > 0):
	    primes.append(i)
    return primes

def find_represent_same_prime_file(file, deg=4):
    r""" 
        Searches a list of admissible number fields until it finds the first discriminant in the list
        for which there are two associated theta series that represent the same prime.

        See Section 5 of Barquero-Sanchez, Mantilla-Soler and Ryan.

    """

    ff = candidate_fields(file,deg)
    for disc in ff:
        mm = len(ff[disc])
        forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(mm)]
        for i in range(mm):
            lli = find_primes_represented(forms[i])
            if len(lli)==0:
	        continue
            ssi = set(lli)
	    for j in range(i+1,mm):
               llj = find_primes_represented(forms[j])
               ssj = set(llj)
	       if len(ssi.intersection(ssj))>0:
	           return disc, i, j 


###############################################################
#                                                             #
#  Find forms with same minimum                               #
#                                                             #
###############################################################



def find_min_num_represented(qexp):
    r"""
	Find the minumum positive integer represented  by a quadratic form Q

        INPUT:
        -  ``qq`` - the q-expansion of the theta series associated to the quadratic form

        EXAMPLES::
	sage: f = x^6 - 10*x^4 - 11*x^3 + 6*x^2 + 8*x + 1
	sage: qf = compute_qf(f)
	sage: theta_series = qf.theta_series(100)
        sage: find_min_num_represented(theta_series)
        10
    """
    for i in range(1,qexp.prec()):
        if qexp[i] > 0:
	    return i
    print "increase prec"
    return


def find_represent_same_min_file(file, deg=4):

    r""" 
        Searches a list of admissible number fields until it finds the first discriminant in the list
        for which there are two associated theta series that have the same minimum

        See Section 5 of Barquero-Sanchez, Mantilla-Soler and Ryan.

    """
    ff = candidate_fields(file,deg)
    for disc in ff:
        mm = len(ff[disc])
    	forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(mm)]
        ll = [find_min_num_represented(f) for f in forms]
	ss = set(ll)
        if len(ss) < len(ll):
	   print disc


##############################################################
#                                                            #
#  Checking for linear independence                          #
#                                                            #
##############################################################


def check_lin_ind(forms):
    r"""

        Returns True of False depending on whether or not the theta-series in the list are linearly independent or not.

        If they are not linearly independent and that can be determined within the precision of the forms, it will return True.

        If they appear to be linearly dependent but you don't have enough coefficients to know it for sure (you would have to up to the Sturm bound to really know), it'll tell you to increase the precision.

	INPUT:
	- ``forms`` a list of q-expansions to the same precision


	EXAMPLES::
  	sage: f = x^4 - 2*x^3 - 9*x^2 + 5*x + 16
	sage: qf = compute_qf(f)
	sage: ts1 = qf.theta_series(50)
	sage: g = x^4 - x^3 - 8*x^2 - 3*x + 4
	sage: qg = compute_qf(g)
	sage: ts2 = qg.theta_series(50)
	sage: ll = [ts1,ts2]
	sage: check_lin_ind(ll)
	True
	sage: ll = [ts1,ts1]
	sage: check_lin_ind(ll)
	increase precision perhaps a lot

    """

    dd = len(forms)
    mat = zero_matrix(dd)
    prec = forms[0].prec()
    if dd > prec:
       print "increase precision to at least ",dd
       return
    for row in range(dd):
        mat[row,0] = 1
    rk = mat.rank()
    col = 1
    for nn in range(1,prec):
        if rk == dd:
	    break
	for row in range(dd):
	   mat[row,col] = forms[row][nn] 
        if mat.rank() > rk:
	   col = col + 1
           rk = mat.rank()
    if rk < dd:
        print "increase precision perhaps a lot"
	return 
    return mat.det() != 0


def check_lin_ind_file(file, deg = 4):
    r"""
        Searches a list of admissible number fields and prints those discriminants for which the associated theta series are not linearly independent

        See Section 5 of Barquero-Sanchez, Mantilla-Soler and Ryan

    """

    ff = candidate_fields(file, deg)
    for disc in ff:
        ll = len(ff[disc])
        forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(ll)]
        lin_ind = check_lin_ind(forms)
	if not lin_ind:
	    print disc


def check_lin_ind_file_malle(file, deg = 4):
    r"""
        Searches a list of admissible number fields and prints those discriminants for which the associated theta series are not linearly independent

        See Section 5 of Barquero-Sanchez, Mantilla-Soler and Ryan

    """

    ff = candidate_fields_malle(file, deg)
    for disc in ff:
        ll = len(ff[disc])
        forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(ll)]
        lin_ind = check_lin_ind(forms)
	if not lin_ind:
	    print disc
