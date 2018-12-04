# lazy_vasp (WIP)
Ph.D. project related scripts. 

## Python

### CPLAPer! 

(or, a Python CPLAP wrapper) automated generation of inputs for https://sourceforge.net/projects/cplap/files/. Should work for binary, ternary and quaternary materials systems. 

Usage:

``` to_cplap( composition_of_interest, competing_phases, elemental_references, dependant_variable ) ```

Will give an `input.dat` file, that should be ready to be fed into CPLAP. Composition of interest, competing phases, and elemental references should be fed in as .yaml files containing stoichiometry and energy information. To do this, use https://github.com/bjmorgan/vasppy/blob/master/scripts/vasp_summary.py. Dependant variable should be fed in as a string of chemical symbol, or 'none'. 

TODO:
 - tests
 - document
 - option to automatically reduce stoichiometry to 1 formula unit
 
 ### SC-Pyper
 
Python wrapper for the Fortran code https://github.com/jbuckeridge/sc-fermi
 
### structure_makers.py

Takes input crystal structure and can give:

 - all symetrically distinct vacancies
 - all symetrically distinct anti-sites

### Para

Testing a script to apply standard vasp band and kpoint parallelism based off number of bands, kpoints and system archetecture. 

### Neb script

Basic NEB post processing

## Shell Scripts

cpcalc.sh: copy a calculation to a new folder, copying CONTCAR to POSCAR

when.sh: grep for calculation timestamp

clean_calc: delete all VASP output files for a 'clean' re-run.


### DEPRECATED

supercell.py
 -takes a poscar and multiplies it in three dimensions to return desired supercell (written with charged defect calculations in mind) 
