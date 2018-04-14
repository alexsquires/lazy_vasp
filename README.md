# lazy_vasp (WIP)
Scripts for things people better at this than me haven't made/already have a better way of doing that I didn't think of. 

## Python

- SC-Fermi interface, to generate input files for https://github.com/jbuckeridge/sc-fermi currently very, very hands on and hacky and needs serious generalisation

- CPLAPer! (or, Python CPLAP wrapper) automated generation of inputs for https://sourceforge.net/projects/cplap/files/. Should work for binary, ternary and quaternary materials systems.

Usage:

``` to_cplap( composition_of_interest, competing_phases, elemental_references, dependant_variable ) ```

Will give an `input.dat` file, that should be ready to be fed into CPLAP. Composition of interest, competing phases, and elemental references should be fed in as .yaml files containing stoichiometry and energy information. To do this, use https://github.com/bjmorgan/vasppy/blob/master/scripts/vasp_summary.py. Dependant variable should be fed in as a string of chemical symbol, or 'none'. 


 ## Shell Scripts

cpcalc.sh
when.sh  


### DEPRECATED

supercell.py
 -takes a poscar and multiplies it in three dimensions to return desired supercell (written with charged defect calculations in mind) 
