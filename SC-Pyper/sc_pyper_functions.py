import re
from vasppy.calculation import *
import subprocess
import pandas as pd
from scipy.stats import linregress
from scipy.constants import physical_constants
import numpy as np

kb_e = physical_constants['Boltzmann constant in eV/K'][0]
j_to_ev = 1/physical_constants['electron volt-joule relationship'][0]
S_0 = 205 * j_to_ev / physical_constants['Avogadro constant'][0]
Cp = (7/2)*kb_e

def dependance(P,T):
    chem_pot = 0.5 * ( (Cp * (T - 298))
                      - T * ( (S_0 + (Cp * np.log(T/298)) + (kb_e * np.log((1/P)) ) ) ))  #### This function gives dependance of mu_O(T,P)
    return chem_pot                                                                       #### BEWARE, CURRENTLY BACKWARDS

E_vbm=0.527

def iden_defect(string):
       return re.findall(r"[^]_[^_]+",string)

def sc_fermi_vacancy_wrap( calc, lattice_site, stoich, limit, cor, charge ):
        formation_energy = ( calc.energy - stoich.energy ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( E_vbm )) + cor
        #print(calc.title, round(formation_energy,1), charge)
        #print(limit,cor)
        return ChargeState(  charge, formation_energy, 3)

def sc_fermi_interstitial_wrap( calc, lattice_site, stoich, limit, cor, charge ):
        formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( E_vbm )) + cor
        #print(calc.title, round(formation_energy,1), charge)
        return ChargeState(  charge, formation_energy, 1)

def sc_fermi_sub_wrap( calc, non_native, native, stoich, non_native_limit, native_limit, cor, charge ):
        formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (charge * ( E_vbm )) + cor
        #print(calc.title, round(formation_energy,1), charge)
        return ChargeState(  charge, formation_energy, 1)

class ChargeState:

    def __init__( self, charge, formation_energy, degeneracy ):
        self.charge = charge
        self.formation_energy = formation_energy
        self.degeneracy = degeneracy

class ChargeState:

    def __init__( self, charge, formation_energy, degeneracy ):
        self.charge = charge
        self.formation_energy = formation_energy
        self.degeneracy = degeneracy


class Defect:

    def __init__( self, label, charge_states, n_sites ):
        self.label = label
        self.charge_states = charge_states
        self.n_sites = n_sites

    @property
    def n_charge_states( self ):
        return len( self.charge_states )

class SCFermi:

    def __init__( self, defects, nelect, e_gap, temperature, spin_polarised=False ):
        self.defects = defects
        self.nelect = nelect
        self.e_gap = e_gap
        self.temperature = temperature
        self.spin_polarised = spin_polarised

    @property
    def n_defects( self ):
        return len( self.defects )

    def output( self ):

            with open('input-fermi.dat', 'w') as f:

                if self.spin_polarised:
                    f.write( '2' + '\n')
                else:
                    f.write( '1' + '\n' )
                f.write( str(self.nelect) + '\n' )
                f.write( str(self.e_gap) + '\n')
                f.write( str(self.temperature) + '\n')
                f.write( str(self.n_defects) + '\n' )
                for d in self.defects:
                    f.write( '{} {} {}'.format( d.label, d.n_charge_states, d.n_sites ) + '\n')
                    for c in d.charge_states:
                        f.write( '{} {} {}'.format( c.charge, c.formation_energy, c.degeneracy ) + '\n')


                f.close()

def make_defect(defects,elements,stoich,corr=0,delta_mu={'mu_Li':0,'mu_Cl':0,'mu_O':0}, sites=1):
    chg_states = []
    for i in defects:
        defect = iden_defect(i.title)
        defect.append('0')
        if defect[0] is 'v':
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich,  delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich,  delta_mu['mu_'+'{}'.format(defect[1])], 0, 0)
        elif defect[1] is 'i':
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], 0, 0)
        else:
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], 0, 0)
        chg_states.append(formation_energy)
        if len(defect) > 4:
            label = str (defect[0] + '_' + defect[1] + '_' + defect[2])
        else:
            label = str (defect[0] + '_' + defect[1])

    return Defect( label, chg_states, sites)

def run_some_fermi(defects,T,nelect, e_gap, spin_polarised):
    out = []
    scf = SCFermi( defects, nelect, e_gap, T, spin_polarised)
    scf.output()
    with open("out.txt", 'w') as f:
        sp = subprocess.run(["./sc-fermi"],stdout=f)
        text_file = open(("out.txt") , "r")
        lines =  text_file.readlines()
        for i in defects:
                print(i.label)
                for line in lines:
                    if re.search(str(i.label)+' '+r'.*?Charge',line) is not None:
                        for h in range(i.n_charge_states):
                            joop = lines[lines.index(line)+(h+1)]
                            coop = joop.split()
                            x = coop[1]
                            y = float(coop[-2])
                            a_dict = {(i.label)+'_'+str(x) : y}
                            out.append(a_dict)
                            print(a_dict)
                            flat = {k: v for d in out for k, v in d.items()}
        return flat


import itertools
import numpy as np

from math import exp, erfc


def get_image_charge_correction(lattice, dielectric_matrix, conv=0.3,
                                factor=30, motif=[0.0, 0.0, 0.0],
                                verbose=False):
    """Calculates the anisotropic image charge correction by Sam Murphy in eV.
    This a rewrite of the code 'madelung.pl' written by Sam Murphy (see [1]).
    The default convergence parameter of conv = 0.3 seems to work perfectly
    well. However, it may be worth testing convergence of defect energies with
    respect to the factor (i.e. cut-off radius).
    References:
        [1] S. T. Murphy and N. D. H. Hine, Phys. Rev. B 87, 094111 (2013).
    Args:
        lattice (list): The defect cell lattice as a 3x3 matrix.
        dielectric_matrix (list): The dielectric tensor as 3x3 matrix.
        conv (float): A value between 0.1 and 0.9 which adjusts how much real
                      space vs reciprocal space contribution there is.
        factor: The cut-off radius, defined as a mutliple of the longest cell
            parameter.
        motif: The defect motif (doesn't matter for single point defects, but
            included in case we include the extended code for defect clusters).
        verbose (bool): If True details of the correction will be printed.
    Returns:
        The image charge correction as {charge: correction}
    """
    inv_diel = np.linalg.inv(dielectric_matrix)
    det_diel = np.linalg.det(dielectric_matrix)
    latt = np.sqrt(np.sum(lattice**2, axis=1))

    # calc real space cutoff
    longest = max(latt)
    r_c = factor * longest

    # Estimate the number of boxes required in each direction to ensure
    # r_c is contained (the tens are added to ensure the number of cells
    # contains r_c). This defines the size of the supercell in which
    # the real space section is performed, however only atoms within rc
    # will be conunted.
    axis = np.array([int(r_c/a + 10) for a in latt])

    # Calculate supercell parallelpiped and dimensions
    sup_latt = np.dot(np.diag(axis), lattice)

    # Determine which of the lattice parameters is the largest and determine
    # reciprocal space supercell
    recip_axis = np.array([int(x) for x in factor * max(latt)/latt])
    recip_volume = abs(np.dot(np.cross(lattice[0], lattice[1]), lattice[2]))

    # Calculatate the reciprocal lattice vectors (need factor of 2 pi)
    recip_latt = np.linalg.inv(lattice).T * 2 * np.pi

    real_space = _get_real_space(conv, inv_diel, det_diel, latt, longest,
                                 r_c, axis, sup_latt)
    reciprocal = _get_recip(conv, inv_diel, det_diel, latt, recip_axis,
                            recip_volume, recip_latt, dielectric_matrix)

    # calculate the other terms and the final Madelung potential
    third_term = -2*conv/np.sqrt(np.pi*det_diel)
    fourth_term = -3.141592654/(recip_volume*conv**2)
    madelung = -(real_space + reciprocal + third_term + fourth_term)

    # convert to atomic units
    conversion = 14.39942
    real_ev = real_space * conversion / 2
    recip_ev = reciprocal * conversion / 2
    third_ev = third_term * conversion / 2
    fourth_ev = fourth_term * conversion / 2
    madelung_ev = madelung * conversion / 2

    correction = {}
    for q in range(1, 8):
        makov = 0.5 * madelung * q**2 * conversion
        lany = 0.65 * makov
        correction[q] = lany

# TODO: Use tabulate
    if verbose:
        logging.info("""
    Results                      v_M^scr    dE(q=1) /eV
    -----------------------------------------------------
    Real space contribution    =  {:.6f}     {:.6f}
    Reciprocal space component =  {:.6f}     {:.6f}
    Third term                 = {:.6f}    {:.6f}
    Neutralising background    = {:.6f}    {:.6f}
    -----------------------------------------------------
    Final Madelung potential   = {:.6f}     {:.6f}
    -----------------------------------------------------""".format(
            real_space, real_ev, reciprocal, recip_ev, third_term, third_ev,
            fourth_term, fourth_ev, madelung, madelung_ev))

        logging.info("""
    Here are your final corrections:
    +--------+------------------+-----------------+
    | Charge | Point charge /eV | Lany-Zunger /eV |
    +--------+------------------+-----------------+""")
        for q in range(1, 8):
            makov = 0.5 * madelung * q**2 * conversion
            lany = 0.65 * makov
            correction[q] = makov
            logging.info("|   {}    |     {:10f}   |    {:10f}   |".
                         format(q, makov, lany))
        logging.info("+--------+------------------+-----------------+")
    return correction


def _get_real_space(conv, inv_diel, det_diel, latt, longest, r_c, axis,
                    sup_latt):
    # Calculate real space component
    real_space = 0.0
    axis_ranges = [range(-a, a) for a in axis]

    # Pre-compute square of cutoff distance for cheaper comparison than
    # separation < r_c
    r_c_sq = r_c**2

    def _real_loop_function(mno):
        # Calculate the defect's fractional position in extended supercell
        d_super = np.array(mno, dtype=float) / axis
        d_super_cart = np.dot(d_super, sup_latt)

        # Test if the new atom coordinates fall within r_c, then solve
        separation_sq = np.sum(np.square(d_super_cart))
        # Take all cases within r_c except m,n,o != 0,0,0
        if separation_sq < r_c_sq and any(mno):
            mod = np.dot(d_super_cart, inv_diel)
            dot_prod = np.dot(mod, d_super_cart)
            N = np.sqrt(dot_prod)
            contribution = 1/np.sqrt(det_diel) * erfc(conv * N)/N
            return contribution
        else:
            return 0.
    real_space = sum(_real_loop_function(mno) for mno in
                     itertools.product(*axis_ranges))
    return real_space


def _get_recip(conv, inv_diel, det_diel, latt, recip_axis, recip_volume,
               recip_latt, dielectric_matrix):
    # convert factional motif to reciprocal space and
    # calculate reciprocal space supercell parallelpiped
    recip_sup_latt = np.dot(np.diag(recip_axis), recip_latt)

    # Calculate reciprocal space component
    axis_ranges = [range(-a, a) for a in recip_axis]

    def _recip_loop_function(mno):
        # Calculate the defect's fractional position in extended supercell
        d_super = np.array(mno, dtype=float) / recip_axis
        d_super_cart = np.dot(d_super, recip_sup_latt)

        if any(mno):
            mod = np.dot(d_super_cart, dielectric_matrix)
            dot_prod = np.dot(mod, d_super_cart)
            contribution = (exp(-dot_prod / (4 * conv**2)) / dot_prod)
            return contribution
        else:
            return 0.
    reciprocal = sum(_recip_loop_function(mno) for mno in
                     itertools.product(*axis_ranges))
    scale_factor = 4 * np.pi / recip_volume
    return reciprocal * scale_factor
