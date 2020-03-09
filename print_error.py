import sys
import importlib.util
from pathlib import Path

from sympy import I


def main(argv):
    """ Print trace norm square of commutator of psi with H, the error squared."""

    if len(argv) < 2 or len(argv) > 3:
        print("Usage: python3 print_error.py XXX_hamil.py psi_name [max_order]")

    # Import the Hamiltonian module from the expanded absolute path into H.
    spec = importlib.util.spec_from_file_location("hamiltonian", str(Path(argv[0]).resolve()))
    H = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(H)
    comm = H.comm

    iofvars = []
    split_orders = []
    normdict = {}
    psi = comm.load_group(Path(argv[1]).resolve(),
                          iofvars=iofvars,
                          split_orders=split_orders,
                          normdict=normdict)
    max_order = int(str(argv[2])) if len(argv) > 2 else len(split_orders-1)
    psi = comm.substitute_group(psi, normdict)
    for order in range(0, max_order+1):
        error = comm.calculate_commutator(H.small, psi[split_orders[order]:split_orders[order+1]])
        #print(psi[split_orders[order]:split_orders[order+1]])
        #print(error)
        error_norm = comm.square_to_find_identity_scalar(error)
        print(comm.mstr(error_norm))
        #print(order)

if __name__ == "__main__":
    main(sys.argv[1:])
