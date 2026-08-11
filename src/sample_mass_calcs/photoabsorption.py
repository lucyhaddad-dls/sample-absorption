"""
Sample properties as related to f2 to go here
"""
import numpy as np
from xraylib import (CompoundParser, EdgeEnergy, AtomicNumberToSymbol,
                     CS_Photo, SymbolToAtomicNumber, AtomicWeight)
from .formulas import _edges

def CS_Photo_Formula(formula:str,
                     absorber:str=None,
                     edge:str=None,
                     energy:float|int|np.ndarray=None)\
                    ->tuple[dict, np.ndarray, np.ndarray]:
    """
    Get mass absorption coefficient [cm^2.g^-1] for a given formula. \
    To convert to an x-ray attenuation length [cm^-1], this needs to be multiplied
    by density. \\
    Then from mu total (mu_t), can find sample thickness via: \
        mu_t = mu_m*density*thickness. = mu_l*thickness \\
    (mu_m = mass absorption coefficient, mu_l = linear absorption coeffient)

    Arguments:
        formula (str): Formula of material.
        absorber (str | None, Optional): Absorbing atom.
        edge (str | None, Optional): Absorption edge (e.g. "L1").
        energy (float | int | np.ndarray | None, Optional): Energy axis\
              (if absorber and edge not provided).

    Returns:
        tuple (tuple[dict, np.ndarray, np.ndarray]): tuple containing:

        photo_dict (dict): Dictionary of key-value pairs for:
            element (str): Nested dictionary for element i of key-value pairs for:
                Z (int): Atomic number.
                A (float): Atomic weight.
                N: Number of elements of type Z.
                massFraction: N/total mass.
                mu_m: Mass absorption coefficient for type Z.
        energy (np.ndarray): Energy [GeV].
        mu_m (np.ndarray): Mass absorption coefficient [cm^2.g^-1].
        

    Note:
        Currently if (absorber, edge) or single energy value is provided,\
        energy is taken to be -100,+200 eV about this value. The range itself \
        is arbitrary and may change but this is to later facilitate edge-dectection \
        steps in a similar way to xafsMass.
    """

    unit = "GeV"
    # update this to work with Measurement objects?

    compound = CompoundParser(formula)
    if edge is not None:
        if absorber is None:
            raise ValueError("Edge and absorber must be provided.")
        energy = EdgeEnergy(SymbolToAtomicNumber(absorber), _edges.index(edge.upper()))
    else:
        if energy is None: raise ValueError("Either edge + absorber or energy must be provided.")

        if isinstance(energy, (float, int)):
            # assume if energy > 1000 it's in eV (otherwise assume GeV)
            if energy*0.001 > 1:
                unit = "eV"
        else:
            if energy[0]*0.001 > 1: unit = "eV"
        if unit == "eV":
            energy *=1e-3; unit = "GeV"

    # make array for energy - 100 and + 200 eV about E:
    if isinstance(energy, (float, int)):
        energy = np.linspace(energy - .1, energy + .2, 100)

    photo_dict = {}
    mu_m = np.zeros_like(energy)

    for i in range(len(compound["Elements"])):
        z = compound["Elements"][i]
        cs_photo = np.empty_like(energy)
        for e in range(len(energy)):
            #  CS_photo gives linear absorption / density , units = [cm^2.gr^-1]
            cs_photo[e] = CS_Photo(z, energy[e])

        a = AtomicWeight(z); n = compound["nAtoms"][i]
        f = compound["massFractions"][i]
        mu_m_tmp = (a*n*cs_photo)/compound["molarMass"]
        
        photo_dict[AtomicNumberToSymbol(z)] = {"Z": z,
                                                "A": a,
                                                "N": n,
                                                "massFraction": f,
                                                "mu_m": mu_m_tmp
                                                }
        mu_m += mu_m_tmp
    # change mu_m and energy to be Measurements?

    return photo_dict, energy, mu_m

