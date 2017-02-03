#!/home/kempj/py343ve/bin/python

from sys import argv
from commutator import load_group, save_to_c_general, substitute_group
from sympy import symbols

iofvars = []
split_orders = []
normdict = {}
f, J2 = symbols('f J2')
psi = load_group(argv[1], iofvars = iofvars, split_orders = split_orders, normdict=normdict)
psi = substitute_group(psi, normdict)
save_to_c_general(psi, [f, J2], argv[2])
