from sympy import I, symbols, simplify
from mathematica_printer import mstr
from commutator import load_group, substitute_group, calculate_commutator, square_to_find_identity
from commutator import Ncproduct, trace_inner_product
import sys

def main(argv):
    N = Ncproduct
    NAME = argv[0]
    ORDER = 0
    END_ORDER = int(argv[1])
    iofvars = []
    split_orders = []
    normdict = {}
    psi = load_group(NAME, iofvars = iofvars, split_orders = split_orders,
                 normdict = normdict)
    sz1 = [Ncproduct(1,[1])]
    psi = substitute_group(psi, normdict)
    for order in range(ORDER, END_ORDER+1):
        norm = trace_inner_product(sz1, psi[split_orders[0]:split_orders[order+1]])
        print(mstr(simplify(norm)))

if __name__ == "__main__":
    main(sys.argv[1:])
