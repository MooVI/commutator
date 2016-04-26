#!/home/kempj/py343ve/bin/python

from sys import argv
from commutator import load_group, save_to_c_poly, substitute_group
from sympy import symbols

iofvars = []
split_orders = []
normdict = {}
f, V = symbols('f V')
psi = load_group(argv[1], iofvars = iofvars, split_orders = split_orders, normdict=normdict)
psi = substitute_group(psi, normdict)
save_to_c_poly(psi, [f, V], argv[2])
