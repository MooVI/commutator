from commutator import load_group, unitary_transform_to_order, save_group
import sys

def main(argv):
    """Convert unitary of form U = e^i(Gs[0] + Gs[1]+...) to operator psi."""
    gs_name = argv[0]
    psi_name = argv[1]
    max_order = int(argv[2])
    zeroth_order = []
    Gs = load_group(gs_name, zeroth_order)
    split_orders = [0, 1]
    psi = zeroth_order[:]
    for order in range(1, max_order+1):
        psi += unitary_transform_to_order(zeroth_order, Gs, order)
        split_orders.append(len(psi))
        save_group(psi, psi_name+'_r'+ str(order), split_orders=split_orders)

if __name__ == "__main__":
    main(sys.argv[1:])
