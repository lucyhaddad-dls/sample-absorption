"""
Different ways to get sample information to formula strings (borrowing from xas_toolbox).
"""
import numpy as np
import re
from xraylib import (CompoundParser, SymbolToAtomicNumber, AtomicNumberToSymbol, AtomicWeight)
from itertools import chain


_edges = ('K','L1','L2','L3','M1','M2','M3','M4','M5')

N_A = 6.022*1e23 # Avogadro const in mol**-1
hbar = 6.582*1e-16 # reduced Planck const in eV.s
r_electron = 2.818*1e-16 # classical electron radius in cm.
hp = 4.136*1e-15 # planck constant 


def reorder_formula(formulae:list[str],
                    formula_string:str)->str:
    """
    Re-order a formula based on the initial input order. \
    This is purely cosmetic!
    """
    pattern = re.compile("[A-Z][a-z]?")
    pattern_num = re.compile("\\d+.\\d+")

    new_elements = re.findall(pattern, formula_string)
    new_conc = re.findall(pattern_num, formula_string)
    old_elements = re.findall(pattern, str(formulae))

    outf = ""

    for element in old_elements:
        new_idx = new_elements.index(element)
        part = "{0:s}{1:0.6}".format(element, new_conc[new_idx])
        outf = outf.join(["", part])
        
    return outf


def formula_from_ratios(formulae:list[str],
                        ratios:list[float|int],
                        keep_order:bool=True)\
                        -> str:
    """
    Get a formula from list of formulae/elements and their
    relative ratios. \\
    Returned formula will have mass fractions as close as possible \
    to the given ratios.

    Arguments:
        formulae (list[str]): List of samples in mixture.
        ratios (list[float | int]): List of ratios of each sample.
        keep_order (bool, Optional): Whether to maintain original order \
        of input in resulting formula.

    Returns:
        out (str): Stoichiometric(?) formula of mixture.
    """
    for r in ratios:
        if r <= 0:
            raise ValueError("Only positive, non-zero ratios valid.")
        
    full_compound = (CompoundParser(f) for f in formulae)

    mass_fractions = []
    nmol = np.inf

    c = 0
    for compound in full_compound:
        for i in range(len(compound["Elements"])):
            Z = compound["Elements"][i]
            frac = np.divide(compound["massFractions"][i]*ratios[c], np.sum(ratios),
                             casting="same_kind")
            moles = frac/AtomicWeight(Z)
            mass_fractions.append((Z, moles))

            if moles < nmol: nmol = moles
        c += 1

    out = ""
    for i in range(len(mass_fractions)):
        out += "{0:s}{1:0.6f}".format(AtomicNumberToSymbol(mass_fractions[i][0]), mass_fractions[i][-1]/nmol)

    if keep_order is True:
        out = reorder_formula(formulae, out)

    return out

def formula_metals_support(support:str, 
                           metal_sites:list[str],
                           metal_loadings:list[float|int])\
                            ->str:
    """
    For N metals on a support make a formula.

    Arguments:
        support (str): Formula of support.
        metal_sites (list[str]): List of metal sites on support.
        metal_loadings (list[float | int]): List of %wt\
              loadings of each site.

    Returns:
        out (str): Chemical formula of mixture.
    """
    atoms = []
    fractions = []
    support = CompoundParser(support)

    for i in range(len(metal_sites)):
        site = metal_sites[i]
        loading = metal_loadings[i]
        comp = CompoundParser(site)
        for n in range(len(comp["Elements"])):
            z = AtomicNumberToSymbol(comp["Elements"][n])
            m = comp["massFractions"][n]
            fractions.append(m*loading/100)
            atoms.append(z)

    for i in range(len(support["Elements"])):
        z = AtomicNumberToSymbol(support["Elements"][i])
        n = support["massFractions"][i]
        fractions.append(n*(100-np.sum(metal_loadings))/100)
        atoms.append(z)

    fractions /=np.min(fractions)

    return formula_from_ratios(atoms, fractions)
