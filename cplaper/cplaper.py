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



def to_cplap( compound_of_interest, competing_phases, elemental_references, dependant_variable ):

    interest = import_calculations_from_file( compound_of_interest )
    
    constutient_elements = len ( cplaper_elements( elemental_references ) )
    
    compound = []
    interest_phase_fomula = []
    res = []
    
    for material in interest:                                 # iterate over the different compounds in the reference phase calculations
        interest_phase_key = str(material)                    # convert the name of the compound to a string
        interest_phase_fomula.append(interest_phase_key)      # make a list of strings of the names of competing phases
    
    for i in interest_phase_fomula:
        x = interest['{}'.format(i)].stoichiometry
        j = {v:k for k,v in x.items()}
        res.append(j)
    
    interest_energy = cplaper_competing(compound_of_interest, elemental_references)
    
    total_competing_phases = len (cplaper_competing(competing_phases, elemental_references))
    
    parameters = []
    
    parameters.append(constutient_elements)
    parameters.append([res, interest_energy])
    parameters.append(dependant_variable)
    parameters.append(total_competing_phases)
    

    with open('params.dat', 'w') as file:
        for item in parameters:
            file.write("%s\n" % item)
    
    f = open('params.dat','r')
    filedata = f.read()
    f.close()
    
    
    remove_formatting = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("{","").replace("}","").replace(":","")


    f = open('cplap_params.dat','w')
    f.write(remove_formatting)
    f.close()
    
    cplap_writer(competing_phases, elemental_references)
    
    cplap_mkinput()
    
    os.remove('cplap_input.dat')
    os.remove('cplap_params.dat')
    os.remove('interim.dat')
    os.remove('params.dat')
    
def cplap_writer( competing_phases, elemental_references ):
    
    elements = cplaper_elements( elemental_references )
    #energies = cplaper_competing( competing_phases, elemental_references )
    
    competing_phase_energies = []  # define empty list
    interim_list = []              # define empty list  
    competing_phase_fomula = []    # define empty list 
    ref_phas = import_calculations_from_file( competing_phases )
    
    for compound in ref_phas:                                 # iterate over the different compounds in the reference phase calculations
        competing_phase_key = str(compound)                   # convert the name of the compound to a string
        competing_phase_fomula.append(competing_phase_key)    # make a list of strings of the names of competing phases
    
    cplap = []
    rad = []
    competing = []
    final = []
    res = []
    
    for i in competing_phase_fomula:
        j = (len (ref_phas['{}'.format(i)].stoichiometry))
        x = ref_phas['{}'.format(i)].stoichiometry
        q = [(v,k) for k,v in x.items()]
                                                                                                # while iterating over the competing phases, 
        for element in elements:                                                                                                   # also iterate within the elements in that competing phase
            single_energy =  ref_phas['{}'.format(i)].stoichiometry['{}'.format(element)] * elements['{}'.format(element)]  # calculate energy per element 
            result = ( ( ref_phas['{}'.format(i)].energy) - single_energy )
        
        rad.append(j)
        rad.append([q, result])
                                            


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
    
    ele_ref  = import_calculations_from_file(elemental_references)     # given .yaml file, elements will read it in as vasppy calculations
    elements = []                                                      # defining an empty list, 'elements'
    for key in ele_ref:                                                # iterating over the keys (ions) in dictionary 'ele_ref' 
        ion = str(key)                                                 # for key in ele_ref, make a string of elemental symbol 
        elements.append(ion)                                           # appened each elemental symbol in turn to list 'elements'   
                                                                       # elements should now be a list of strings of the ions in elemental references
    
    elemental_reference_energies = []                                                                                     # define a new empty list, elemental                                   
    for ion in elements:                                                                                                  # iterating over the ions in elements lists
        single_element_energy = float (ele_ref['{}'.format(ion)].energy / ele_ref['{}'.format(ion)].stoichiometry[ion] )  # divide calculation energy by number of ions to get energy per ion  
        individual_ion_dict = { ion : single_element_energy }                                                             # for each ion, create a dictionary of the symbol and the energy per ion                  # y       
        elemental_reference_energies.append(individual_ion_dict)               
    normalised_elemental_references = { k: v for d in elemental_reference_energies for k, v in d.items() }                # convert the individual dictionaries to one dictionary              
    
    return normalised_elemental_references    
    
def cplaper_competing( competing_phases, elemental_references ):
 
    competing_phase_energies = []  # define empty list
    interim_list = []              # define empty list  
    competing_phase_fomula = []    # define empty list 

    elements = cplaper_elements(elemental_references)          # take .yaml elemental reference file and get energy per ion
    ref_phas = import_calculations_from_file(competing_phases) # read in .yaml competing phases as calculations
    
    for compound in ref_phas:                                 # iterate over the different compounds in the reference phase calculations
        competing_phase_key = str(compound)                   # convert the name of the compound to a string
        competing_phase_fomula.append(competing_phase_key)    # make a list of strings of the names of competing phases
    
    for compound in competing_phase_fomula:                                                                                        # while iterating over the competing phases, 
        for element in elements:                                                                                                   # also iterate within the elements in that competing phase
            single_energy =  ref_phas['{}'.format(compound)].stoichiometry['{}'.format(element)] * elements['{}'.format(element)]  # calculate energy per element 
            interim_list.append(single_energy)
            iterateable_list = iter(interim_list)
            k = [a + b + c + d for a,b,c,d in zip(iterateable_list,iterateable_list,iterateable_list,iterateable_list)]            # sum energy per element to get total energy
            
    for compound,i in zip (competing_phase_fomula,k):
         
        result = ( ( ref_phas['{}'.format(compound)].energy) - i )
            
        competing_phase_energies.append(result)    
    
          
    
    return competing_phase_energies