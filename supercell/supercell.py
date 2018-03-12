import argparse

from pymatgen import Structure, Lattice
from pymatgen.io.vasp.inputs import Poscar

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='POSCAR', type=str,
                        help='path to input file')
    parser.add_argument('-o', '--output', default='poscar',
                        help='output file format')
    parser.add_argument('-g', '--grid', nargs='+', type=int)

    args = parser.parse_args()

    struct = Structure.from_file(args.file)

    supercell = struct * (args.grid)

    supercell.to(filename='supercell'.format(args.file),fmt=args.output)



if __name__ == "__main__":
    main()
