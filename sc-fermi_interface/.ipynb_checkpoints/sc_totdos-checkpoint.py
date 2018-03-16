import argparse
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
import numpy as np

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='vasprun.xml', type=str, help='path to input file')
    parser.add_argument('-o', '--output', default='totdos.dat', help='output file format')
    parser.add_argument('-s', '--spin', help='yes, if spin polarised, otherwise leave blank')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    
    args = parse_command_line_arguments()
    
    dosrun = Vasprun(args.file)
    totdens = dosrun.tdos.densities[Spin.up] + dosrun.tdos.densities[Spin.down]
    
    totdos = dosrun.tdos.energies -  dosrun.eigenvalue_band_properties[2]    #Set VBM to 0 eV
    
    
    if args.spin is not None:                                                                                   # Use this if dos is spin polarised
    sc_input = np.column_stack((totdos, dosrun.tdos.densities[Spin.up], dosrun.tdos.densities[Spin.down]))  # Create array of energy against spin polarised density
    np.savetxt("totdos.dat", sc_input)                                                                      # Save as totdos.dat
    
    
    else:
    sc_input = np.column_stack((totdos, totdens))                            #Create Array of Energy against density
    np.savetxt("totdos.dat", sc_input)                                       #Save as totdos.dat file

    

