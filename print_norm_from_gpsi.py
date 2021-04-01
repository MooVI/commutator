from sympy import I, symbols, simplify
from mathematica_printer import mstr
from commutatorbitarray import load_group, substitute_group, calculate_commutator, square_to_find_identity
from commutatorbitarray import trace_inner_product
import sys
import pickle

def main(argv):
    """Print norm of psi."""
    if len(argv) < 1 or len(argv) > 2:
        print("Usage: python3 print_error.py psi_name [max_order]")
    NAME = argv[0]
    iofvars = []
    split_orders = []
    normdict = {}
    psi = load_group(NAME,
                     iofvars=iofvars,
                     split_orders=split_orders,
                     normdict=normdict)
    max_order = int(str(argv[1])) if len(argv) > 2 else len(split_orders)-2
    zeroth_order = psi[split_orders[0]:split_orders[1]]
    psi = substitute_group(psi, normdict)
    for order in range(max_order+1):
        norm = trace_inner_product(zeroth_order, psi[split_orders[0]:split_orders[order+1]])
        filename = f'norm_{order}_{NAME}.p'
        with open(filename, "wb" ) as fz:
            pickle.dump(norm, fz, protocol=pickle.HIGHEST_PROTOCOL)
        #print(mstr(norm), flush=True)

if __name__ == "__main__":
    main(sys.argv[1:])
