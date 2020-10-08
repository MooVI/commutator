import numpy as np
from commutator import NCSum, full_simplify_group, NCProduct, conjugate_group, hermitian_conjugate_group
from commutator import hermitian_trace_inner_product as inner
from commutator import SigmaProduct as Sg
import sympy
import sympy

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

    def canonical_form(self, simplify=True):
        if simplify:
            for ix,iy in np.ndindex(self.shape()):
                self[ix,iy] = full_simplify_group(self[ix,iy])
        self._remove_zero_rowcols()
        i = 1
        while i <= self.chi():
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

            s = sympy.sqrt(sum(inner(self[row,i], self[row,i]) for row in range(i)))

            if sympy.simplify(s) == 0:
                self._delrowcol(i)
            else:
                self[:i,i] *= 1/s
                self[i,i:] *= s
                i += 1
        self._remove_zero_rowcols()

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
    ret.canonical_form()
    return ret


def mul_nparray(ncsum, nparray):
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
    sqchi = first.chi()*other.chi()
    retchi = 2*(other.chi()+first.chi()) + sqchi
    ret = MPOzeros([retchi+2,retchi+2])
    ret[0,0] = first[0,0]*other[0,0]
    ret[-1,-1] = first[-1,-1]*other[-1,-1]
    ret['d'] = first['d']*other['d']
    ret['c'] = np.concatenate(
        (other['c'],
         first['c'],
         np.kron(first['c'], other['c']),
         rmul_nparray(first['c'],other['d']),
         mul_nparray(first['d'],other['c'])),
        axis = 1)
    zerocol = MPOzeros((first.chi()+other.chi(),1)).data
    ret['b'] = np.concatenate(
        (zerocol,
        np.kron(first['b'], other['b']),
         first['b'],
         other['b']),
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
        [other['A'], z21,      np.kron(first['c'], other['A']), np.kron(first['c'], other['b']), z22                           ],
        [z12,        first['A'],np.kron(first['A'], other['c']), z11                         , np.kron(first['b'], other['c'])],
        [zs2,        zs1,      np.kron(first['A'], other['A']), np.kron(first['A'], other['b']), np.kron(first['b'], other['A'])],
        [z12,        z11,      z1s                           , first['A']                   , z12                           ],
        [z22,        z21,      z2s                           , z21                         , other['A']                    ],
    ])

    return ret
