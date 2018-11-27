import pymatgen as pmg
from pymatgen.io.vasp import Poscar
from bsym.interface.pymatgen import unique_structure_substitutions
from collections import Counter

def mk_cancies(struc):
    arg = pmg.Structure.from_file(struc)
    for i,j in arg.composition.items():
        vac = unique_structure_substitutions(arg, '{}'.format(i), {'{}'.format(i):(int(j)-1),'V':1})
        for p,q in zip(vac,range(len(vac))):
            p.remove_species('V')
            p.to(filename='POSCAR_V_'+'{}'.format(q)+'{}'.format(i))

def mk_antis(struc):
    
    arg = pmg.Structure.from_file(struc)
    arg.add_oxidation_state_by_guess()
    coop = arg.as_dict()
    
    anions = ([i['species'][0]['element'] for i in coop['sites'] if i['species'][0]['oxidation_state'] < 0])
    cations = ([i['species'][0]['element'] for i in coop['sites'] if i['species'][0]['oxidation_state'] > 0])
    
    anion = Counter(anions)
    cation = Counter(cations)
    
    for i,j in anion.items():
        for x in anion.keys():
            if x is not i:
                sub = unique_structure_substitutions(arg, '{}'.format(i), {'{}'.format(i):(int(j)-1),'{}'.format(x):1})
                for p,q in zip(sub,range(len(sub))):
                    p.to(fmt='POSCAR', filename='POSCAR_'+'{}'.format(x)+'_'+'{}'.format(i)+'_'+'{}'.format(q))
    for i,j in cation.items():
        for x in cation.keys():
            if x is not i:
                sub = unique_structure_substitutions(arg, '{}'.format(i), {'{}'.format(i):(int(j)-1),'{}'.format(x):1})
                for p,q in zip(sub,range(len(sub))):
                    p.to(fmt='POSCAR', filename='POSCAR_'+'{}'.format(x)+'_'+'{}'.format(i)+'_'+'{}'.format(q))
    