"""
sample class (?) to go here
"""
from .photoabsorption import CS_Photo_Formula
from .formulas import _edges
from xraylib import SymbolToAtomicNumber, AtomicWeight, EdgeEnergy
from .edge_step import get_edge_vals, get_bounds, fit_bkg
import numpy as np

class Element:
    """
    Class for x-ray absorption data related to a 
    single element in a sample.
    """
    def __init__(self, name:str, **xray_vals):
        """
        Create `Element`, makes attributes:
            Z (int): Atomic number.
            A (float | int): Atomic weight.
            N (int): Number of `Element` in parent sample.
            mass_fraction (float | int): Mass fraction of `Element` in parent sample.
            energy (np.ndarray | float | int): Energy axis for element's mass absorption\
                  value(s).
            mass_absorption (np.ndarray | float | int): Mass absorption coefficient(s)\
                  for `Element`.
        """
        self.name = name

        if "Z" in xray_vals.keys():
            self.Z = xray_vals["Z"]
        else:
            self.Z = SymbolToAtomicNumber(self.name)

        if "A" in xray_vals.keys():
            self.A = xray_vals["A"]
        else:
            self.A = AtomicWeight(self.Z)
        
        if "N" in xray_vals.keys():
            self.N = xray_vals["N"]
        else: self.N = 1

        if "massFraction" in xray_vals.keys():
            self.mass_fraction = xray_vals["massFraction"]
        else: self.mass_fraction = 1

        if "energy" in xray_vals.keys():
            self.energy = xray_vals["energy"]
        else:
            self.energy = None

        if "mu_m" in xray_vals.keys():
            self.mass_absorption = xray_vals["mu_m"]
        else:
            self.mass_absorption = None
            print("Input edge/energy to calculate mass absorption coefficient.")


    def mass_absorption_calc(self, 
                             edge:str|None=None, 
                             energy:float|int|None=None):
        """
        Calculate mass absorption for single atom.
        """
        out = CS_Photo_Formula(self.name, 
                               self.name,
                               edge,
                               energy)
        
        self.mass_absorption = out["mu_m"]
        self.energy = out["energy"]
        self.edge = edge

    def get_jump_height(self,
                        bounds:tuple[int, int],
                        pre,
                        edge,
                        include:bool=True):
  
        self.step, self.abs_max, self.abs_min = fit_bkg(self.energy, 
                                                        self.mass_absorption, 
                                                        bounds,
                                                        pre, edge,
                                                        include)


class XRaySample:
    def __init__(self, 
                 formula:str,
                 absorber:str|None = None,
                 edge:str|None = None,
                #  energy:np.ndarray|float|int|None = None,
                 density:float|int|None = None,
                 surface_density:float|int|None=None,
                 area:float|int|None = None,
                 volume:float|int|None = None,
                 thickness:float|int|None = None,
                 mu_total:float|int|None = 2.6):
        """
        NOTE: single/multiple energy inputs not yet supported.
        """

        out = CS_Photo_Formula(formula,
                               absorber,
                               edge)
        self.unit = "GeV"
        self.energy = out["energy"]
        self.mass_absorption = out["mu_m"]
        self.mu_total = mu_total
        self.edge = edge
        self.absorber = absorber

        if self.edge is not None:
            self.edge_energy = EdgeEnergy(SymbolToAtomicNumber(absorber),
                                          _edges.index(self.edge.upper()))

        self.density = density
        self.surface_density = surface_density
        self.area = area
        self.volume = volume
        self.thickness = thickness

        # need to put in the case where energy provided and edges could be
        # multiple..
        
        # list of Elements for individual atom-type info.
        self.elements = []
        for k in out.keys():
            if k not in ["energy", "mu_m"]:
                out[k]["energy"] = self.energy
                out[k]["edge"] = self.edge
                self.elements.append(Element(k, **out[k]))


    def calulate_step_sum(self):
        """
        Calculate step in absorption coefficient at an edge
        energy.
        """
        absorber = [e for e in self.elements if e.name == self.absorber][0]
        
        edge, bounds, pre = get_bounds(self.energy,
                                        absorber.mass_absorption,
                                        self.edge_energy)
        
        self.step, self.abs_max, self.abs_min = fit_bkg(self.energy, 
                                                        self.mass_absorption, 
                                                        bounds,
                                                        pre,
                                                        edge,
                                                        include=True)

    def calculate_total_step(self):
        """
        Calculate absorption step based on mass, volume, thickness.
        """
        self.abs_step = self.step*(self.density)*self.thickness
        # this needs to be multiplied by density to get mu (mu_T = mu*thickness)
           
    def calculate_thickness(self):
        """
        Calculate sample thickness given: \\
        `density`, `mu_total`, `formula`
        """
        if self.mu_total is not None:
            if self.density is None:
                self.calculate_density()
            self.thickness = (self.mu_total/self.abs_max*self.density)


    def calculate_mass(self):
        if self.volume is not None:
            if self.density is not None:
                self.mass = self.volume * self.density
            if self.density is None:
                self.calculate_density()

    def calculate_density(self):
        if self.mass is not None:
            if self.volume is not None:
                self.density = self.mass/self.volume
            else:
                if self.area is not None:
                    self.surface_density = self.mass/self.area
  
        if self.surface_density is not None and self.thickness is not None:
                self.density = self.surface_density/self.thickness


    def calculate_mu_total():
        pass