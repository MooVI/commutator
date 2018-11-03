#!/home/kempj/py343ve/bin/python

from sys import argv
from commutator import load_group, print_group, square_to_find_identity_scalar_up_to_order, substitute_group
from commutator import collect_terms
from commutator import Ncproduct as N
from sympy import latex, expand
from mathematica_printer import mstr

#Note squaring this is nonsense, however
def square_to_find_identity_even(group):
    """Assumes Hermitian and all terms even in number"""
    result = 0
    for a in collect_terms(group):
        result += (a.scalar**2)*(1-len(a.product)%4)
    return expand(result)

iofvars = []
split_orders = []
normdict = {}
psi = load_group(argv[1], iofvars = iofvars, split_orders = split_orders, normdict=normdict)
psi = substitute_group(psi, normdict)
psi = [N(a.scalar,a.product[1:]) for a in psi if 1 in a.product]
prop = square_to_find_identity_even(psi)
print(str(prop)+'\n\n')
print(mstr(prop)+'\n\n')
print(latex(prop).replace('\\\\', '\\'))
