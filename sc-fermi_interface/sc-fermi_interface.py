import argparse
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
import numpy as np

def parse_command_line_arguments():
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='vasprun.xml', type=str, help='path to vasprun')
    parser.add_argument('-s', '--spin', help='yes, if spin polarised, otherwise leave blank')
    parser.add_argument('-t', '--temperature', default='1000', type=str, help='Temperature in Kelvin?')
    parser.add_argument('-n', '--nelect', type=str, help='Number of Electrons')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    
    args = parse_command_line_arguments()
    
    vasprun = Vasprun(args.file)
    bandgap = vasprun.eigenvalue_band_properties[0]    # bandgap from vasprun 

    if args.spin is not None:            
        line1 = '2'
        line2 = args.nelect
        line3 = str(bandgap)
        line4 = args.temperature
                
    else:
        line1 = '1'
        line2 = args.nelect
        line3 = str(bandgap)
        line4 = args.temperature 
        

    lines = [line1, line2, line3, line4]
    with open('sc-fermi-params.dat', 'w') as f:
    	f.writelines("%s\n" % l for l in lines) 


        
