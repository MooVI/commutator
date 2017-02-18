from commutator import *
import sys

def main(argv):
    gs_name = argv[0]
    psi_name = argv[1]
    max_order = int(argv[2])
    Gs = load_group(gs_name)
    split_orders = [0,1]
    sz1 = [Ncproduct(1, [1])]
    psi = sz1[:]
    for order in range(1, max_order+1):
        psi += unitary_transform_to_order(sz1, Gs, order)
        split_orders.append(len(psi))
        save_group(psi, psi_name+'_r'+ str(order), split_orders = split_orders)

if __name__ == "__main__":
    main(sys.argv[1:])
