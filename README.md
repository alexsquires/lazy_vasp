# lazy_vasp (WIP)
Ph.D. project related scripts. 

## Python

### SC-Fermi interface
To generate input files for https://github.com/jbuckeridge/sc-fermi currently very, very hands on and hacky and needs serious generalisation

### CPLAPer! 

#### Currently some issues in how the input file is complied.

(or, Python CPLAP wrapper) automated generation of inputs for https://sourceforge.net/projects/cplap/files/. Should work for binary, ternary and quaternary materials systems.

Usage:

``` to_cplap( composition_of_interest, competing_phases, elemental_references, dependant_variable ) ```

Will give an `input.dat` file, that should be ready to be fed into CPLAP. Composition of interest, competing phases, and elemental references should be fed in as .yaml files containing stoichiometry and energy information. To do this, use https://github.com/bjmorgan/vasppy/blob/master/scripts/vasp_summary.py. Dependant variable should be fed in as a string of chemical symbol, or 'none'. 

TODO:
 - tests
 - document
 - option to automatically reduce stoichiometry to 1 formula unit

## Shell Scripts

cpcalc.sh
when.sh  


### DEPRECATED

supercell.py
 -takes a poscar and multiplies it in three dimensions to return desired supercell (written with charged defect calculations in mind) 
