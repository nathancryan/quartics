def trace_zero_subspace(f):
    K.<a> = NumberField(f) # construct the number field
    O_K = K.maximal_order()
    bas = O_K.basis() # find a basis
    dim = len(bas)
    M = Matrix(ZZ,dim,dim,[[x.trace() for x in bas] for _ in range(dim)]) # calculate the matrix representation of the trace map in the basis
    kerM = (M.transpose()).kernel() # find a basis for subspace of elements of trace zero
    #print kerM
    return [ sum([v[i]*bas[i] for i in range(dim)]) for v in kerM.basis()]

def trace_zero_matrix(f):
    bas = trace_zero_subspace(f) # find a basis for subspace of elements of trace zero
    #print bas
    ll = len(bas)
    
    MF = matrix(ZZ,ll,ll,[[(bas[i]*bas[j]).trace() for i in range(ll)] for j in range(ll)])  
    return MF

def qf_matrix_to_list(mat):
    ret = []
    for i in range(mat.nrows()):
    	for j in range(i,mat.nrows()):
	    ret.append(mat[i][j])
    return ret

def compute_qf(f):
    M = trace_zero_matrix(f)
    return QuadraticForm(M)

##################################################################################################

def sturm_bound(qf):
    N = qf.level()
    mu = Gamma0(N).index()
    kappa = qf.dim()
    return floor(kappa*mu/24 - (mu-1)/N)


##################################################################################################


def candidate_fields(file,deg=4):
    attach(file)
    good_fields = []
    for l in data:
        if (gcd(l[Integer(1)],Integer(deg))==Integer(1) and l[Integer(1)].is_squarefree() and l[1]>0):
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


#################################################################################################

def find_prime_exponents(qq):
    prime_list = []
    for i in range(qq.prec()):
    	if qq[i] != 0 and is_prime(i):
            prime_list.append(i)
    return prime_list


#################################################################################################

PREC = 10000


def check_lin_ind(forms):
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
    ff = candidate_fields(file, deg)
    for disc in ff:
        ll = len(ff[disc])
        forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(ll)]
	print disc, check_lin_ind(forms)
	

def find_min_num_represented(qexp):
    for i in range(1,qexp.prec()):
        if qexp[i] > 0:
	    return i
    print "increase prec"
    return


def find_represent_same_min_file(file, deg=4):
    ff = candidate_fields(file,deg)
    for disc in ff:
        mm = len(ff[disc])
    	forms = [compute_qf(ff[disc][i][0]).theta_series(PREC) for i in range(mm)]
        ll = [find_min_num_represented(f) for f in forms]
	ss = set(ll)
        if len(ss) < len(ll):
	   print disc

def find_primes_represented(qexp):
    primes = []
    for i in range(qexp.prec()):
        if (is_prime(i) and qexp[i] > 0):
	    primes.append(i)
    return primes

def find_represent_same_prime_file(file, deg=4):
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