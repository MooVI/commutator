class MPO():
    def __init__(self, data):
        if isinstance(data, list):
            self.data = np.empty((len(data), len(data[0])), dtype=np.object)
            for i, row in enumerate(data):
                for j, el in enumerate(row):
                    self.data[i, j] = NCSum(el)
        else:
            self.data = np.empty(data.shape, dtype=np.object)
            for ix,iy in np.ndindex(data.shape):
                self.data[ix, iy] = NCSum(data[ix, iy])


    def __getitem__(self, ind):
        if isinstance(ind, str):
            if ind == 'A':
                return self.data[1:-1, 1:-1]
            elif ind == 'b':
                return self.data[1:-1, -1:]
            elif ind == 'c':
                return self.data[0:1, 1:-1]
            elif ind == 'd':
                return self.data[0,-1]
            else:
                KeyError("Unrecognised key " + ind)
        return self.data[ind]

    def __setitem__(self, ind, value):
        if isinstance(ind, str):
            if ind == 'A':
                self.data[1:-1, 1:-1] = value
            elif ind == 'b':
                self.data[1:-1, -1:] = value
            elif ind == 'c':
                self.data[0:1, 1:-1] = value
            elif ind == 'd':
                self.data[0,-1] = value
            else:
                KeyError("Unrecognised key " + ind)
        else:
            self.data[ind] = value

    def chi(self):
        return self.data.shape[0]-2


    def shape(self):
        return self.data.shape

    def __str__(self):
        return str(self.data)

    def __add__(self, other):
        retchi = self.chi()+other.chi()
        ret = MPOzeros([retchi+2,retchi+2])
        ret[0,0] = self[0,0]
        ret[-1,-1] = self[-1,-1]
        ret['c'] = np.concatenate((self['c'], other['c']), axis = 1)
        ret['b'] = np.concatenate((self['b'], other['b']), axis = 0)
        ret['A'] = np.block([
            [self['A'], MPOzeros([self.chi(), other.chi()]).data],
            [MPOzeros([other.chi(), self.chi()]).data, other['A']]
        ])
        ret['d'] = self['d'] + other['d']
        #else:
        #    TypeError("Cannot add iMPOs with top left or bottom right elements unequal.")
        return ret


    def __mul__(self, other):
        """Disjoint multiply"""
        if isinstance(other, MPO):
            return multiply_MPOS(self, other)
        if isinstance(other, (NCProduct, NCSum, list, str, np.ndarray)):
            raise TypeError("Convert to MPO before multiplying by NCProduct")
        else:
            ret = MPO(self.data)
            ret['b'] *= other
            ret['d'] *= other
            return ret

    def __rmul__(self, other):
        """Disjoint multiply"""
        if isinstance(other, MPO):
            return multiply_MPOS(other, self)
        if isinstance(other, (NCProduct, NCSum, list, str, np.ndarray)):
            raise TypeError("Convert to MPO before multiplying by NCProduct")
        else:
            ret = MPO(self.data)
            ret['b'] *= other
            ret['d'] *= other
            return ret

    def _delrowcol(self,i):
        self.data = np.delete(self.data, i, axis=0)
        self.data = np.delete(self.data, i, axis=1)


    def conjugate(self):
        ret = MPO(self.data)
        for ix,iy in np.ndindex(ret.shape()):
            ret[ix, iy] = conjugate_group(self[ix, iy])
        return ret


    def hermitian_conjugate(self):
        ret = MPO(self.data)
        for ix,iy in np.ndindex(ret.shape()):
            ret[ix, iy] = hermitian_conjugate_group(self[ix, iy])
        return ret


    def _remove_zero_rowcols(self):
        i = 1
        while i < self.shape()[0]-1:
            if np.all(np.equal(self[:, i],0)):
                self._delrowcol(i)
            else:
                i += 1
        j = 1
        while j < self.shape()[1]-1:
            if np.all(np.equal(self[j, :],0)):
                self._delrowcol(j)
            else:
                j += 1

    def left_canonical_form(self, simplify=True):
        if simplify:
            for ix,iy in np.ndindex(self.shape()):
                self[ix,iy] = full_simplify_group(self[ix,iy])
        self._remove_zero_rowcols()
        i = 1
        while i <= self.chi():
            print(i)
            r = []
            c = [sum(inner(self[row,j], self[row,i]) for row in range(j+1)) for j in range(i)]
            r.append(c[0]/(1-inner(self[0,0],self[i,i])))
            j = 1
            while j < i:
                Kj = [- inner(self[col, j], self[i,i]) for col in range(j)]
                rj = (c[j]-sum(Kj[col]*r[col] for col in range(j)))/(1-inner(self[j,j],self[i,i]))
                r.append(rj)
                j += 1



            for row in range(i):
                self[row, i] -= sum(r[col]*self[row,col] for col in range(i))

            for j in range(i):
                self[j, i:] += r[j]*self[i,i:]

            s = sympy.together(sympy.sqrt(sum(inner(self[row,i], self[row,i]) for row in range(i))/(1-inner(self[i,i],self[i,i]))))



            print(s)

            zero = False
            try:
                if s < 1e-12:
                    zero = True
            except:
                n, d = s.as_numer_denom()
                if n.subs({f:0.1}) < 1e-12:
                    zero = True

            if zero:
                print("deleting: " + str(i))
                self._delrowcol(i)
            else:
                self[:i,i] *= (1/s)
                self[i,(i+1):] *= s
                i += 1

            if simplify:
                for ix, iy in np.ndindex(self[:(i+1),i:].shape):
                    self[ix,iy] = self[ix,iy].simplify_scalars()


        self._remove_zero_rowcols()







    def right_canonical_form(self, simplify=True):
        self.data = np.rot90(self.data,2).T
        if simplify:
            for ix,iy in np.ndindex(self.shape()):
                self[ix,iy] = full_simplify_group(self[ix,iy])
        self._remove_zero_rowcols()
        i = 1
        while i <= self.chi():
            print(i)
            r = []
            c = [sum(inner(self[row,j], self[row,i]) for row in range(j+1)) for j in range(i)]
            r.append(c[0]/(1-inner(self[0,0],self[i,i])))
            j = 1
            while j < i:
                Kj = [- inner(self[col, j], self[i,i]) for col in range(j)]
                rj = (c[j]-sum(Kj[col]*r[col] for col in range(j)))/(1-inner(self[j,j],self[i,i]))
                r.append(rj)
                j += 1



            for row in range(i):
                self[row, i] -= sum(r[col]*self[row,col] for col in range(i))

            for j in range(i):
                self[j, i:] += r[j]*self[i,i:]

            s = sympy.cancel(sympy.sqrt(sum(inner(self[row,i], self[row,i]) for row in range(i))/(1-inner(self[i,i],self[i,i]))))


            print(s)

            zero = False
            try:
                if s < 1e-12:
                    zero = True
            except:
                n, d = s.as_numer_denom()
                if sympy.expand(n) == 0:
                    zero = True

            if zero:
                print("deleting: " + str(i))
                self._delrowcol(i)
            else:
                self[:i,i] *= (1/s)
                self[i,(i+1):] *= s
                i += 1

            if simplify:
                for ix, iy in np.ndindex(self[:(i+1),i:].shape):
                    self[ix,iy] = self[ix,iy].simplify_scalars()


        self._remove_zero_rowcols()
        self.data = np.rot90(self.data,2).T


    def expand(self, N):
        L = self.shape()[0]
        left = np.zeros(L, dtype = int)
        left[0] = 1
        right = np.zeros((L,1), dtype = int)
        right[-1,0] = 1
        ret = np.dot(left, self.data)
        print(left)
        print(ret)
        for i in range(1,N):
            trans = np.empty(self.data.shape, dtype=np.object)
            for ix,iy in np.ndindex(trans.shape):
                trans[ix,iy] = self.data[ix,iy].translate(i)
            ret = np.dot(ret, trans)
        return np.dot(ret, right)[0]

def MPOzeros(shape):
    return MPO(np.zeros(shape, dtype=int))

def to_iMPO(group):
    mpo = MPOzeros((2,2)).data
    typ = (type(group[0]))
    mpo[0,0] = typ(1,[])
    mpo[-1,-1] = typ(1,[])
    for el in group:
        supp = el.support_indices()
        if len(supp) == 0:
            ValueError("Cannot iMPOize term with bare identity.")
        if supp[0] != 1:
            ValueError("Cannot iMPOize term with support not on first site.")
        if len(supp) == 1:
           mpo[0,-1] += el
        else:
            old = mpo.shape[0]
            extra = supp[-1]-1
            zerocol = MPOzeros((old, extra)).data
            mpo = np.concatenate((mpo[:,:-1], zerocol, mpo[:,-1:]), axis=1)
            mpo[0, old-1] = el.get_operator_on_site(1)
            zerorow = MPOzeros((extra, mpo.shape[1])).data
            mpo = np.concatenate((mpo[:-1,:], zerorow, mpo[-1:,:]), axis=0)
            for i in range(old, old+extra):
                mpo[i-1, i] = el.get_operator_on_site(i).translate(old-1-i)
            mpo[old+extra-2, old+extra-1] *= el.scalar
    ret = MPO(mpo)
    return ret


def mul_nparray(mul1, mul2):
    if isinstance(mul1, NCSum):
        return lmul_nparray(mul1, mul2)
    else:
        return rmul_nparray(mul1, mul2)

def lmul_nparray(ncsum, nparray):
    ret = np.empty(nparray.shape, dtype=np.object)
    for ix,iy in np.ndindex(nparray.shape):
        ret[ix, iy] = ncsum*nparray[ix, iy]
    return ret

def rmul_nparray(nparray, ncsum):
    ret = np.empty(nparray.shape, dtype=np.object)
    for ix,iy in np.ndindex(nparray.shape):
        ret[ix, iy] = nparray[ix, iy]*ncsum
    return ret


def multiply_MPOS(first, other):
    mul = mul_nparray
    sqchi = first.chi()*other.chi()
    retchi = 2*(other.chi()+first.chi()) + sqchi
    ret = MPOzeros([retchi+2,retchi+2])
    ret[0,0] = first[0,0]*other[0,0]
    ret[-1,-1] = first[-1,-1]*other[-1,-1]
    ret['d'] = first['d']*other['d']
    ret['c'] = np.concatenate(
        (
            mul(first[0,0],other['c']),
            mul(first['c'],other[0,0]),
            np.kron(first['c'], other['c']),
            mul(first['c'],other['d']),
            mul(first['d'],other['c'])
        ),
        axis = 1)
    ret['b'] = np.concatenate(
        (
            mul(first['d'],other['b']),
            mul(first['b'],other['d']),
            np.kron(first['b'], other['b']),
            mul(first['b'],other[-1,-1]),
            mul(first[-1,-1], other['b'])
        ),
        axis = 0)

    z11 = MPOzeros((first.chi(), first.chi())).data
    z12 = MPOzeros((first.chi(), other.chi())).data
    z21 = MPOzeros((other.chi(), first.chi())).data
    z22 = MPOzeros((other.chi(), other.chi())).data
    #zss = MPOzeros((retchi, retchi)).data
    z1s = MPOzeros((first.chi(), sqchi)).data
    zs1 = MPOzeros((sqchi, first.chi())).data
    z2s = MPOzeros((other.chi(), sqchi)).data
    zs2 = MPOzeros((sqchi, other.chi())).data

    ret['A'] = np.block([
        [mul(first[0,0],other['A']),z21,                        np.kron(first['c'], other['A']), np.kron(first['c'], other['b']), z22                            ],
        [z12,                       mul(first['A'],other[-1,1]),np.kron(first['A'], other['c']), z11                            , np.kron(first['b'], other['c'])],
        [zs2,                       zs1,                        np.kron(first['A'], other['A']), np.kron(first['A'], other['b']), np.kron(first['b'], other['A'])],
        [z12,                       z11,                        z1s                            , mul(first['A'], other[-1,-1])  , z12                            ],
        [z22,                       z21,                        z2s                            , z21                            , mul(first[-1,-1],other['A'])   ],
    ])

    return ret

def column_overlaps(impo):
    overlap = np.empty(impo.shape(), dtype=np.object)
    for ix,iy in np.ndindex(overlap.shape):
        overlap[ix, iy] = sum(inner(impo[row, ix], impo[row, iy]) for row in range(overlap.shape[0]))
    return overlap

def row_overlaps(impo):
    overlap = np.empty(impo.shape(), dtype=np.object)
    for ix,iy in np.ndindex(overlap.shape):
        overlap[ix, iy] = sum(inner(impo[ix, col], impo[iy, col]) for col in range(overlap.shape[1]))
    return overlap

def inverse_pauli_vector(vec):
    idpart = inner(Sg(1,""), vec)*Sg(1,"")
    rvec = vec - 2*idpart
    return rvec*(1/inner(rvec,vec))

# def __idfree_gauge(impo):
#     chi = impo.chi()
#     s = [inverse_pauli_vector(impo("A")[-1,-1] - Sg(1, ""))*impo["b"][-1, 0]]
#     j = -2
#     while j > 0:
#         sj = inverse_pauli_vector(impo("A")[j,j] - Sg(1, "")*(impo["b"][j,0]-sum(impo["A"][j, col]*s[col] for col in range(j,0)))
#         s.insert(0, sj)
#         j -=1
#     ret = MPO(impo)


def idfree_gauge(impo):
    Id = [Sg(1,[])]
    chi = impo.chi()
    s = impo.chi()*[0]
    s[-1] = inner(Id, impo["b"][-1, 0])/(inner(Id,impo["A"][-1,-1]) - 1)
    j = impo.chi()-2
    while j >= 0:
        print(s)
        print(j)
        print(impo.chi())
        s[j] = (inner(Id,impo["b"][j,0])-sum(inner(Id, impo["A"][j, col])*s[j] for col in range(j+1,impo.chi())))/(inner(Id, impo["A"][j,j]) - 1)
        j -=1
    ret = MPO(impo.data)
    ret[0, -1] += -sum(s[col]*ret[0, col+1] for col in range(ret.chi()))
    for row in range(1,ret.chi()+1):
        ret[row, -1] += -sum(s[col]*ret[row, col+1] for col in range(row-1, ret.chi())) + s[row-1]*ret[-1,-1]
    ret._remove_zero_rowcols()
    return ret

def nosqrtleft_canonical_form(impo, fvars, simplify=True):
        ret = MPO(impo.data)
        if simplify:
            for ix,iy in np.ndindex(ret.shape()):
                ret[ix,iy] = full_simplify_group(ret[ix,iy])
        ret._remove_zero_rowcols()
        norms = [sum(inner(el, el) for el in ret[0,:1])]
        i = 1
        while i <= ret.chi():
            print(i)
            r = []
            c = [sum(inner(ret[row,j], ret[row,i]) for row in range(j+1)) for j in range(i)]
            try:
                r.append(c[0]/(norms[0]-inner(ret[0,0],ret[i,i])))
            except:
                r.append(0)
            j = 1
            while j < i:
                Kj = [- inner(ret[col, j], ret[i,i]) for col in range(j)]
                rj = (c[j]-sum(Kj[col]*r[col] for col in range(j)))/(norms[j]-inner(ret[j,j],ret[i,i]))
                r.append(rj)
                j += 1

            for row in range(i):
                ret[row, i] -= sum(r[col]*ret[row,col] for col in range(i))

            for j in range(i):
                ret[j, i:] += r[j]*ret[i,i:]

            s =  sum(inner(el, el) for el in ret[:i,i])

            print(s)

            zero = is_zero(s, fvars)

            if zero:
                print("deleting: " + str(i))
                ret._delrowcol(i)
            else:
                norms.append(s+inner(ret[i,i],ret[i,i]))
                i += 1

            if simplify:
                for ix, iy in np.ndindex(ret[:(i+1),i:].shape):
                    ret[ix,iy] = ret[ix,iy].simplify_scalars()


        ret._remove_zero_rowcols()
        return ret

def is_zero(s, fvars):
    try:
        if s < 1e-12:
            return True
    except:
        #n, d = s.as_numer_denom()
        test = s.xreplace(dict(zip(fvars, [random.random() for _ in fvars]))).evalf(chop=True)
        print(test)
        if sympy.re(test) < 1e-12:
            return True
    return False

def remove_irrev_rowcols(impo):
    ret = MPO(impo.data)
    i = 1
    dim = ret.shape()[0]
    while i < dim - 1 :
        if np.all(np.equal(ret[:i, i],0)):
            ret._delrowcol(i)
        else:
            i += 1
    j = dim - 2
    while j > 0:
        if np.all(np.equal(ret[j, (j+1):],0)):
            ret._delrowcol(j)
        else:
            j -= 1

    return ret

def Rfree_gauge(impo):
    Id = NCSum(Sg(1,[]))
    chi = impo.chi()
    s = impo.chi()*[0]
    s[-1] = inverse_pauli_vector(impo["A"][-1,-1] - Id)*impo["b"][-1, 0]
    j = impo.chi()-2
    while j >= 0:
        print(s)
        print(j)
        print(impo.chi())
        s[j] = inverse_pauli_vector(impo["A"][j,j] - Id)*(impo["b"][j,0]-sum(impo["A"][j, col]*s[j] for col in range(j+1,impo.chi())))
        j -=1
    ret = MPO(impo.data)
    ret[0, -1] += -sum(ret[0, col+1]*s[col] for col in range(ret.chi()))
    for row in range(1,ret.chi()+1):
        ret[row, -1] += -sum(ret[row, col+1]*s[col] for col in range(row-1, ret.chi())) + s[row-1]*ret[-1,-1]
    ret._remove_zero_rowcols()
    return ret

Id = NCSum(Sg(1,"Id"))
fvargen = sympy.numbered_symbols("a", real = True)
fvars = []
H= to_iMPO(Sg(1,"z1 z2")+ Sg(f, "x1"))
H.left_canonical_form()
mpos = []
sols = []
for mpo in mpo_gen(3, fvargen, fvars):
    mpos.append(mpo)
    print(mpo)
    mul = H*mpo
    comm = mul + (mul*(-1)).hermitian_conjugate()
    comm = remove_irrev_rowcols(comm)
    Rfreecomm = Rfree_gauge(comm)
    sols.append(sympy.solve((sympy.simplify(el.scalar) for el in Rfreecomm[0,-1]), fvars, dict = True))
    print(sols[-1])


import itertools
def pauli_gen(fvar):
    ps = ["Id", "x1", "y1", "z1"]
    for p in ps:
        yield Sg(fvar, p)

def pauli_noid_gen(fvar):
    ps = ["x1", "y1", "z1"]
    for p in ps:
        yield Sg(fvar, p)

def edge_mpo_gen(dim, fvargen, fvars):
    newfvars = [next(fvargen) for _ in range(dim*(dim+1)//2-3)]
    fvars += newfvars
    ps = [["Id", "x1", "y1", "z1"]]
    ps_noid = [["x1", "y1", "z1"]]
    pscomb = ps_noid*(dim-2)
    for i in range(1, dim-1):
        pscomb += ps_noid + (dim-2-i)*ps
    for pslist in itertools.product(*pscomb):
        ips = iter(pslist)
        ifv = iter(fvars)
        first_row = [[0]+[Sg(next(ifv), next(ips)) if j != (dim-2) else Sg(1, "z1") for j in range(dim-1)]]
        middle = [[0 if j < i else Sg(next(ifv), next(ips)) if j != (dim-1) else Sg(next(ifv), "z1")for j in range(dim)] for i in range(1, dim-1)]
        last_row = [[0]*(dim-1) + [Sg(1,[])]]
        yield MPO(first_row + middle + last_row)

def mpo_gen(dim, fvargen, fvars):
    newfvars = [next(fvargen) for _ in range(dim*(dim+1)//2-3)]
    fvars += newfvars
    ps = [["Id", "x1", "y1", "z1"]]
    ps_noid = [["x1", "y1", "z1"]]
    pscomb = ps_noid*(dim-2)
    for i in range(1, dim-1):
        pscomb += ps_noid + (dim-2-i)*ps
    for pslist in itertools.product(*pscomb):
        ips = iter(pslist)
        ifv = iter(fvars)
        first_row = [[Sg(1, "Id")]+[Sg(next(ifv), next(ips)) if j != (dim-2) else Sg(1, "z1") for j in range(dim-1)]]
        middle = [[0 if j < i else Sg(next(ifv), next(ips)) if j != (dim-1) else Sg(next(ifv), "z1")for j in range(dim)] for i in range(1, dim-1)]
        last_row = [[0]*(dim-1) + [Sg(1,[])]]
        yield MPO(first_row + middle + last_row)


def mma_solve(eqs, fvars):
    syms = set([sym for eq in eqs for sym in eq.free_symbols])
    syms.update(fvars)
    mathematica_parser.vardict = dict(zip([str(sym) for sym in syms], syms))
    streqs = '{' + ','.join([mstr(eq) + '==0' for eq in eqs]) + '}'
    strfvars  = '{' + ','.join([mstr(fvar) for fvar in fvars]) + '}'
    print(streqs)
    solstr = mma.evaluate(wl.ToString(wlexpr(f'Solve[{streqs},{strfvars}, Domain=Reals]'), wl.InputForm))
    print(solstr)
    sols = []
    if solstr == '{}':
        return sols
    for sol in solstr[1:-2].split('}, '):
        sols.append((dict(),[]))
        var = sol.find('->')
        for rule in sol[1:].split(', '):
            var, ans  = rule.split(' -> ')
            cond = 'ConditionalExpression'
            if ans[:len(cond)] == cond:
                ans, cstr = ans[len(cond)+1:
            sols[-1][0][mathematica_parser.vardict[var]] = mmaparser(ans, mathematica_parser.vardict)
    return sols
