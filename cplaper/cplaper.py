import vasppy
from vasppy.calculation import *
import os

def to_cplap( compound_of_interest, competing_phases, elemental_references, dependant_variable ):

    interest = import_calculations_from_file( compound_of_interest )
    ele_ref  = import_calculations_from_file( elemental_references )
    comps    = import_calculations_from_file( competing_phases )
    
    elements = cplap_elements(ele_ref) 
    
    cplap_mkinput ( cplap_interest( interest, elements, comps, dependant_variable), cplap_competing(comps, elements) )
    
def cplap_interest( interest, elements, competing_phases, dependant_variable):  

    """
    Compiles interest phase information in a CPLAP-friendly format
    
    Args:
        
        interest (dict(vasppy.Calculation)): material to be considered
        competing_phases (dict(vasppy.Calculation)): competing phases to be considered
        elements (dict{Str:Float}): dictionary of {elemental symbol : energy per atom}
        dependant_variable (str): element CPLAP should consider as the dependant variable 
        
    Returns:
        Stoichiometry and energy formatted suitable for CPLAP
    
    """
    interest_phase_fomula = []
    interest_info = []

    for material in interest:                                
        interest_phase_key = str(material)                    
        interest_phase_fomula.append(interest_phase_key)      
    
    for i in interest_phase_fomula:
        number_of_elements = (len (interest['{}'.format(i)].stoichiometry))
        stoich = interest['{}'.format(i)].stoichiometry
        formatted_stoich = [(v,k) for k,v in stoich.items()]
        interest_info.append([formatted_stoich, cplap_energy(interest['{}'.format(i)], elements)])
    
    interest_parameters = []
    
    interest_parameters.append(len( elements ))
    interest_parameters.append(interest_info)
    interest_parameters.append(dependant_variable)
    interest_parameters.append(len( competing_phases ))
    
    return interest_parameters


def cplap_competing( competing_phases, elements ):
     
    
    """
    Compiles competing phase information in a CPLAP-friendly format
    
    Args:
        competing_phases (dict(vasppy.Calculation)): competing phases to be considered
        elements (dict{Str:Float}): dictionary of {elemental symbol : energy per atom}
        
    Returns:
        Stoichiometry and energy formatted suitable for CPLAP
    
    """
    
    competing_phase_fomula = []    
    
    for compound in competing_phases:                                 
        competing_phase_key = str(compound)                 
        competing_phase_fomula.append(competing_phase_key) 
    
    competing_phase_info = []

    for i in competing_phase_fomula:
        number_of_elements = (len (competing_phases['{}'.format(i)].stoichiometry))
        stoich = competing_phases['{}'.format(i)].stoichiometry
        formatted_stoich = [(v,k) for k,v in stoich.items()]
                                                                                                
        competing_phase_info.append(number_of_elements)
        competing_phase_info.append([formatted_stoich, cplap_energy(competing_phases['{}'.format(i)], elements)])
        
    return competing_phase_info
 
def cplap_mkinput( interest_info, competing_info ):
    
    """
    Compiles interest_info and competing_info into the input file
    
    Args:
        interest_info (str): CPLAP-formatted information about material of interest
        competing_info (str): CPLAP-formatted information about competing phases
    
    """
    
    with open('interim.dat', 'w') as file:
        for item in interest_info:
            file.write("%s\n" % item)
        for item in competing_info:
            file.write("%s\n" % item)
        

    f = open('interim.dat','r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("(","").replace(")","").replace(":","")


    f = open('input.dat','w')
    f.write(newdata)
    f.close()
    
    os.remove('interim.dat')



def cplap_elements( elemental_references ):
    
    """
    Takes a set of element vasppy.Calculation objects and returns a dictionary of scaled energies (per atom)
    
    Args:
        elemental_references (dict(vasppy.Calculations)): vasppy calculations for elements to be considered
        
    Returns:
        normalised_elemental_references (dict{Str:Float}): dictionary of {elemental symbol : energy per atom}
    
    """
    
    elements = []     
    elemental_reference_energies = []
    
    for key in elemental_references:                                   
        ion = str(key)                                                  
        elements.append(ion)                                              
                                                                    
    for ion in elements:                                                                                                  
        single_element_energy = float (elemental_references['{}'.format(ion)].energy 
                                       / elemental_references['{}'.format(ion)].stoichiometry[ion] )                        
        individual_ion_dict = { ion : single_element_energy }                                                                              
        elemental_reference_energies.append(individual_ion_dict)               
    
    normalised_elemental_references = { k: v for d in elemental_reference_energies for k, v in d.items() }                             
    return normalised_elemental_references    


def cplap_energy(compound, elemental_references):
    
    """
    Calulates the formation energy of a compound
    
    Args: 
        compound (vasppy.Calculation): compound for which formation energy should be calculated
        elemental_references (dict{str:float}): dictionary of elemental symbol:energy
    
    Returns:
        (float) formation energy: formation energy
    
    """
    individual_energies = []
    
    for element in elemental_references:
            energy_per_element = compound.stoichiometry['{}'.format(element)] * elemental_references['{}'.format(element)] 
            individual_energies.append(energy_per_element)
            sum_of_elemental_energies = sum(individual_energies)
            formation_energy =( ( compound.energy) - sum_of_elemental_energies )
   
    return formation_energy
