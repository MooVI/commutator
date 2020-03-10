import sys
import importlib.util
from pathlib import Path

def main(argv):
    """ Checks the commutator with H is the claimed order in perturbation theory.
    This is maximally suspicious check, unlike print_error.py which assumes this is true
    already for optimisation. This means it might be slow..."""

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
    psi = comm.load_group(argv[1],
                          iofvars=iofvars,
                          split_orders=split_orders,
                          normdict=normdict)
    max_order = int(str(argv[2])) if len(argv) > 2 else len(split_orders)-2
    psi = comm.substitute_group(psi, normdict)
    for order in range(0, max_order+1):
        print("At order " + str(order) + ": " +
              str(comm.check_group_at_least_order(
                  comm.calculate_commutator(H.H, psi[:split_orders[order+1]]),
                  order+1,
                  H.orders))
        )

if __name__ == "__main__":
    main(sys.argv[1:])
