#!/home/kempj/py343ve/bin/python

from sys import argv
from commutator import load_group, print_group, substitute_group
import importlib


iofvars = []
split_orders = []
normdict = {}
psi = load_group(argv[1], iofvars = iofvars, split_orders = split_orders, normdict=normdict)
psi = substitute_group(psi, normdict)
head = importlib.import_module(argv[2])
print_group(head.c(head.H, psi))
