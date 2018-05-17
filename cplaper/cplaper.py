import vasppy
from vasppy.calculation import *
import fileinput
import csv
from itertools import islice
import os

def cplap_mkinput():
    filenames = ['cplap_params.dat', 'cplap_input.dat' ]
    with open('input.dat', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())

    os.remove('cplap_input.dat')
    os.remove('cplap_params.dat')
    os.remove('interim.dat')
    os.remove('params.dat')

def to_cplap( compound_of_interest, competing_phases, elemental_references, dependant_variable ):

    interest = import_calculations_from_file( compound_of_interest )
    ele_ref  = import_calculations_from_file( elemental_references )
    comps    = import_calculations_from_file( competing_phases )
    
    elements = cplaper_elements(ele_ref)
    
    
    constutient_elements   = len( ele_ref )
    total_competing_phases = len( comps )
    
    compound = []
    interest_phase_fomula = []
    interest_info = []

    
    for material in interest:                                 # iterate over the different compounds in the reference phase calculations
        interest_phase_key = str(material)                    # convert the name of the compound to a string
        interest_phase_fomula.append(interest_phase_key)      # make a list of strings of the names of competing phases
    
    for i in interest_phase_fomula:
        number_of_elements = (len (interest['{}'.format(i)].stoichiometry))
        stoich = interest['{}'.format(i)].stoichiometry
        formatted_stoich = [(v,k) for k,v in stoich.items()]
                                                                                                                                   # while iterating over the competing phases, 

        #interest_info.append(number_of_elements)
        interest_info.append([formatted_stoich, cplaper_energy(interest['{}'.format(i)], elements)])
    
    parameters = []
    
    parameters.append(constutient_elements)
    parameters.append(interest_info)
    parameters.append(dependant_variable)
    parameters.append(total_competing_phases)
    

    with open('params.dat', 'w') as file:
        for item in parameters:
            file.write("%s\n" % item)
    
    f = open('params.dat','r')
    filedata = f.read()
    f.close()
    
    
    remove_formatting = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("{","").replace("}","").replace(":","").replace("(","").replace(")","")


    f = open('cplap_params.dat','w')
    f.write(remove_formatting)
    f.close()
    
    cplap_writer(comps, ele_ref)
    
    cplap_mkinput()
    

    
def cplap_writer( competing_phases, elemental_references ):
    
    elements = cplaper_elements( elemental_references )
    
    
    
    competing_phase_energies = []  # define empty list
    interim_list = []              # define empty list  
    competing_phase_fomula = []    # define empty list 
    
    
    for compound in competing_phases:                                 # iterate over the different compounds in the reference phase calculations
        competing_phase_key = str(compound)                   # convert the name of the compound to a string
        competing_phase_fomula.append(competing_phase_key)    # make a list of strings of the names of competing phases
    
    cplap = []
    rad = []
    competing = []
    final = []
    res = []
    blank = []
    
    for i in competing_phase_fomula:
        number_of_elements = (len (competing_phases['{}'.format(i)].stoichiometry))
        stoich = competing_phases['{}'.format(i)].stoichiometry
        formatted_stoich = [(v,k) for k,v in stoich.items()]
                                                                                                                                   # while iterating over the competing phases, 

        rad.append(number_of_elements)
        rad.append([formatted_stoich, cplaper_energy(competing_phases['{}'.format(i)], elements)])
                                            


    with open('interim.dat', 'w') as file:
        for item in rad:
            file.write("%s\n" % item)
        

    f = open('interim.dat','r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("(","").replace(")","").replace(":","")


    f = open('cplap_input.dat','w')
    f.write(newdata)
    f.close()


def cplaper_elements( elemental_references ):
    
         
    elements = []                                                      # defining an empty list, 'elements'
    for key in elemental_references:                                   # iterating over the keys (ions) in dictionary 'elemental_references' 
        ion = str(key)                                                 # for key in ele_ref, make a string of elemental symbol 
        elements.append(ion)                                           # appened each elemental symbol in turn to list 'elements'   
                                                                       # elements should now be a list of strings of the ions in elemental references
    
    elemental_reference_energies = []                                                                                     # define a new empty list, elemental                                   
    
    for ion in elements:                                                                                                  # iterating over the ions in elements lists
        single_element_energy = float (elemental_references['{}'.format(ion)].energy 
                                       / elemental_references['{}'.format(ion)].stoichiometry[ion] )                      # divide calculation energy by number of ions to get energy per ion  
        individual_ion_dict = { ion : single_element_energy }                                                             # for each ion, create a dictionary of the symbol and the energy per ion                  # y       
        elemental_reference_energies.append(individual_ion_dict)               
    
    normalised_elemental_references = { k: v for d in elemental_reference_energies for k, v in d.items() }                # convert the individual dictionaries to one dictionary              
    return normalised_elemental_references    


def cplaper_energy(compound, elemental_references):
    
    blank = []
    

    for element in elemental_references:                                                                                                   # also iterate within the elements in that competing phase
            single_energy = compound.stoichiometry['{}'.format(element)] * elemental_references['{}'.format(element)] 
            blank.append(single_energy)
            sums = sum(blank)
            energy =( ( compound.energy) - sums )
    return energy