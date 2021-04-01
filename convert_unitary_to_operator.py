import sys
from commutatorbitarray import load_group, unitary_transform_to_order, save_group, NCSum


def main(argv):
    """Convert unitary of form U = e^i(Gs[0] + Gs[1]+...) to operator psi."""
    if len(argv) < 2 or len(argv) > 3:
        print("Usage: python3 print_error.py unitary_name output_name [max_order]")
    gs_name = argv[0]
    psi_name = argv[1]
    zeroth_order = NCSum([])
    Gs = load_group(gs_name, zeroth_order=zeroth_order)
    max_order = int(str(argv[2])) if len(argv) > 2 else len(Gs)
    split_orders = [0, 1]
    psi = NCSum(zeroth_order[:])
    for order in range(1, max_order+1):
        psi.add_no_simplify(unitary_transform_to_order(zeroth_order, Gs, order))
        split_orders.append(len(psi))
        save_group(psi, psi_name+'_r'+ str(order), split_orders=split_orders)

if __name__ == "__main__":
    main(sys.argv[1:])
