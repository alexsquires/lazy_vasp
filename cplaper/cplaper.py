from vasppy.calculation import import_calculations_from_file

class Phase:

    def __init__( self, formula, energy):
        self.formula = formula
        self.energy = energy
        
    @property
    def n_elements( self ):
        return len( self.formula )
    
class Cplap_in:

    def __init__( self, competing, dependent, interest):
        self.competing = competing
        self.dependent = dependent
        self.interest = interest

    @property
    def n_competing( self ):
        return len( self.competing )

    def output( self ):
            print(self.interest)

            with open('input.dat', 'w') as f:
                for i in self.interest:
                    f.write( str(int(i.n_elements/2)) + '\n' )
                    f.write( str(i.formula) + ' ' + str(i.energy) + '\n')
                f.write( str(self.dependent) + '\n')
                f.write( str(self.n_competing) + '\n' )
                for p in self.competing:
                    f.write( str(int(p.n_elements/2)) + '\n')
                    f.write( str(p.formula) + ' ' + str(p.energy) + '\n')
                
            f.close()

def formation_energy(phase,elements,calc_type):
        single_el_energies = {i : elements[i].energy/elements[i].stoichiometry[i] for i in elements}
        competing_and_reference = [single_el_energies,dict(calc_type[phase].stoichiometry)]
        collected_info = {
            k: [d.get(k) for d in competing_and_reference]
            for k in set().union(*competing_and_reference)
        }
        formation_energy = (calc_type[phase].energy - sum([i*j for i,j in collected_info.values() if j is not None]))
        return formation_energy
            
            
def mk_entry(phase,calc_type,elements):
    format_formula = [(v,k) for k,v in phase.stoichiometry.items()]
    flat_formula = [g for j in format_formula for g in j ]
    energy = formation_energy(phase.title,elements,calc_type)
    n_els = len(phase.stoichiometry.items())
    return Phase(flat_formula,energy)
            
def mk_input(phase_of_interest,competing_phases,elemental_references,dependent):
    inter = import_calculations_from_file(phase_of_interest)
    comp = import_calculations_from_file(competing_phases)
    els = import_calculations_from_file(elemental_references)
    interest = [mk_entry(inter[i],inter,els) for i in inter]
    competing = [mk_entry(comp[i],comp,els) for i in comp]
    Cplap_in(competing,dependent,interest).output()
    tidy()
           

def tidy():
    # Read in the file
    with open('input.dat', 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace("'","").replace("[","").replace("[","").replace("]","").replace(",","").replace("(","").replace(")","").replace(":","")

    # Write the file out again
    with open('input.dat', 'w') as file:
        file.write(filedata)       

