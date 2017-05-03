#!/home/kempj/py343ve/bin/python

import sys, getopt
from commutator import load_group, print_group, substitute_group, collect_terms, convert_group, order_group, texify_group, remove_zeros
from sympy import symbols, Symbol

def main(argv):
    try:
        opts,args = getopt.getopt(argv, 'scvtbofzu')
    except getopt.GetoptError:
        print ("usage: print_options.py [-scvtbofzu]")
        sys.exit(2)
    iofvars = []
    split_orders = []
    normdict = {}
    sub = False
    collect = False
    convert = False
    texify = False
    extract_free = False
    bare = False
    order = False
    zero_free = False
    unitary = False
    for opt, arg in opts:
        if opt =='-s':
            sub = True
        if opt == '-c':
            collect = True
        if opt == '-v':
            convert = True
        if opt == '-t':
            texify = True
            bare = False
        if opt == '-b':
            bare = True
            texify = False
        if opt == '-o':
            order = True
        if opt == '-f':
            extract_free = True
        if opt == '-z':
            zero_free = True
        if opt == '-u':
            unitary = True
            
    if sub or zero_free:
        psi = load_group(argv[1], iofvars = iofvars, normdict=normdict)
    elif unitary:
        gs = load_group(argv[1])
        psi = []
        for g in gs:
            psi+=g
    else:
        psi = load_group(argv[1], iofvars = iofvars)
    if extract_free:
        psi = [el for el in psi if any(i in iofvars for i in el.scalar.atoms(Symbol))]
    if sub:
        psi = substitute_group(psi, normdict)
    if zero_free:
        psi = substitute_group(psi, dict(zip(iofvars,[0]*len(iofvars))))
        remove_zeros(psi)
    if collect:
        psi = collect_terms(psi)
        remove_zeros(psi)
    if convert:
        psi = convert_group(psi)
    if order:
        V, f, V1, V2, X, Y, Vy = symbols('V f V1 V2 X Y Vy')
        orders = {V:1,f:1,V1:1,V2:1, X:1, Y:1, Vy:1}
        psi = order_group(psi, orders)

    if texify:
        print('\\documentclass{article}\n'
              '\\usepackage{amsmath, amssymb, graphics, setspace}\n'
              '\\allowdisplaybreaks\n'
              '\\begin{document}')
        print(texify_group(psi, newlines = True))
        print('\\end{document}')
    elif bare:
        for el in psi:
            print(repr(el))
    else:
        print_group(psi)

if __name__ == "__main__":
    main(sys.argv[1:])
