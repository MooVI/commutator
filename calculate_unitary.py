import sys
import importlib.util
from pathlib import Path

from sympy import I


START_ORDER = 1
END_ORDER = 8

START_GS = []

def main(argv):
    """ Calculate unitary need to transform zeroth order such that it
    is conserved by H order by order in small. Hamiltonian can be either Paulis
    or Majoranas.
    """

    if len(argv) != 2:
        print("Usage: python3 calculate_unitary.py XXX_hamil.py output_name")

    # Import the Hamiltonian module from the expanded absolute path into H.
    spec = importlib.util.spec_from_file_location("hamiltonian", str(Path(argv[0]).resolve()))
    H = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(H)

    #Additional setup
    FILEHEAD = argv[1]
    comm = H.comm
    c = comm.calculate_commutator
    Gs = START_GS
    H.large = comm.NCSum(H.large)
    H.small = comm.NCSum(H.small)

    if START_GS:
        psi_order = comm.unitary_transform_to_order(H.zeroth_order, Gs, START_ORDER-1)
    else:
        psi_order = comm.NCSum(H.zeroth_order[:])

    #Solve order by order for the unitary
    for test_order in range(START_ORDER, END_ORDER+1):

        print("Solving for term to cancel commutator order: {}.".format(test_order))

        Hcomm = c(H.small, psi_order)
        psi_order = comm.unitary_transform_to_order(H.zeroth_order,
                                                    Gs,
                                                    test_order,
                                                    not_single_comm=True)
        Hcomm = c(H.large, psi_order) + Hcomm

        fvarname = 'F' + str(test_order) + '_'
        iofvars = []

        g_comm = comm.solve_commutator_equation(H.large,
                                                Hcomm,
                                                H.vardict,
                                                iofvars=iofvars,
                                                fvarname=fvarname,
                                                verbose=True)

        print("Found term to add at order {}. Finding unitary...".format(test_order))

        inv = comm.solve_single_term_commutator_equation(H.zeroth_order,
                                                         g_comm,
                                                         iofvars,
                                                         H.vardict,
                                                         verbose=True)
        Gs.append(I*inv)

        psi_order.add_no_simplify(I*comm.commute_group(H.zeroth_order, Gs[-1]))

        comm.save_group(Gs, FILEHEAD + '_r' + str(test_order), zeroth_order=comm.NCSum(H.zeroth_order))


print('Done!')

if __name__ == "__main__":
    main(sys.argv[1:])
