import argparse
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='vasprun.xml', type=str,
                        help='path to input file')
    parser.add_argument('-o', '--output', default='totdos.dat',
                        help='output file format')
    parser.add_argument('-s', '--spin', help='yes, if spin polarised, otherwise leave blank')

    args = parser.parse_args()
    dosrun = Vasprun(args.file)
    
    
    totdos = dosrun.tdos.energies - dosrun.efermi
    
    sc_input = np.column_stack((totdos, dosrun.tdos.densities[Spin.up]))

    np.savetxt("totdos.dat", sc_input)

    if args.spin is not None:
        sc_input = np.column_stack((totdos, dosrun.tdos.densities[Spin.up], dosrun.tdos.densities[Spin.up]))
        np.savetxt("totdos.dat", sc_input)

if __name__ == "__main__":
    main()
