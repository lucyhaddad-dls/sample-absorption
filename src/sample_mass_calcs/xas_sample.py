# set all input params as properties:
from .photoabsorption import CS_Photo_Formula
from xraylib import EdgeEnergy, SymbolToAtomicNumber
import numpy as np
from .formulas import _edges
from .edge_step import get_bounds, fit_bkg
from .measurements import Unit, units, Measurement

# do i want to make two objects with different logic for (possibly) multiple-edges and
# single edges?

class PhotoElement:
    def __init__(self,
                 name:str,
                 Z:int,
                 A:int|float,
                 N:int,
                 massFraction:float,
                 mu_m:np.ndarray,
                 mass_unit:str|Unit,
                 length_unit:str|Unit
                 ):
        """
        read input from CS_Photo_Formula and make 
        it easier to query per-atom information..?
        """
    
        # properties
        if isinstance(mass_unit, Unit):
            self._mass_unit = mass_unit
        elif isinstance(mass_unit, str):
            self._mass_unit = getattr(units, mass_unit)

        if isinstance(length_unit, Unit):
            self._length_unit = length_unit
        elif isinstance(length_unit, str):
            self._length_unit = getattr(units, length_unit)

        # measurements
        self.mass_absorption = Measurement(value = mu_m,
                                                _unit = units.cm**2/units.g)
        if self.mass_absorption.unit != self.length_unit**2 / self.mass_unit:
                    self.mass_absorption.to(self.length_unit**2 / self.mass_unit)

        # attributes
        self.name = name
        self.Z = Z
        self.A = A
        self.N = N
        self.massFraction = massFraction

    @property
    def length_unit(self):
        return self._length_unit

    @length_unit.setter
    def length_unit(self, val:str|Unit):
        if isinstance(val, str):
            self._length_unit = getattr(units, val)
        else:
            self._length_unit = val
        self.mass_absorption.to(self.length_unit**2/self.mass_unit)
        return self._length_unit

    @property
    def mass_unit(self):
        return self._mass_unit

    @mass_unit.setter
    def mass_unit(self, val:str|Unit):
        if isinstance(val, str):
            self._mass_unit = getattr(units, val)
        else:
            self._mass_unit = val
        self.mass_absorption.to(self.length_unit**2/self.mass_unit)
        return self._mass_unit


    def make_linear_absorption(density:Measurement):
        """
        Linear absorption coefficient = mass absorption coefficient * density.
        """
        pass

    def scale_linear_absorption(thickness:Measurement):
        """
        mu_total = linear absorption coefficient * thickness at value after
        the edge.
        """
        pass
      
class XRaySample:
    def __init__(self,
                formula:str,
                absorber:str|None = None,
                edge:str|None = None,
                density:float|int|None = None,
                surface_density:float|int|None = None,
                mass:float|int|None = None,
                area:float|int|None = None,
                thickness:float|int|None = None,
                mu_total:float|int|None = 2.6,
                mass_unit:Unit|str = "g",
                length_unit:Unit|str = "cm",
                energy_unit:Unit|str = "gev"):
        """
        """
        self.formula = formula
        self.absorber = absorber
        self.edge = edge
        self.mu_total = mu_total

        if isinstance(mass_unit, Unit):
            self._mass_unit = mass_unit
        elif isinstance(mass_unit, str):
            self._mass_unit = getattr(units, mass_unit)

        if isinstance(length_unit, Unit):
            self._length_unit = length_unit
        elif isinstance(length_unit, str):
            self._length_unit = getattr(units, length_unit)

        if isinstance(energy_unit, Unit):
            self._energy_unit = energy_unit
        elif isinstance(energy_unit, str):
            self._energy_unit = getattr(units, energy_unit)

        # do i want these all to be properties so that 
        # they can all be calculated wrt. each other?
        self.thickness = Measurement(value = thickness, _unit = self.length_unit)
        self.area = Measurement(value = area, _unit = self.length_unit**2)

        self.mass = Measurement(value = mass, _unit=self.mass_unit)
        self.surface_density = Measurement(value = surface_density, _unit = self.mass_unit/self.length_unit**2)
        self.density = Measurement(value = density, _unit = self.mass_unit/self.length_unit**3)
        
        if self.absorber != None and self.edge != None:
            self.elements = []
            atoms_dict, energy, mass_absorption = CS_Photo_Formula(self.formula,
                                   self.absorber,
                                   self.edge)
            

            edge_energy = EdgeEnergy(SymbolToAtomicNumber(self.absorber),
                                     _edges.index(self.edge.upper()))
            self.edge_energy = Measurement(value=edge_energy, _unit= units.gev)

            if self.edge_energy.unit != self.energy_unit:
                self.edge_energy.to(self.energy_unit)

            for element, values in atoms_dict.items():
                self.elements.append(PhotoElement(
                                                  name = element,
                                                  **values,
                                                  mass_unit= self.mass_unit,
                                                  length_unit=self.length_unit))
        else:
            raise NotImplementedError("Energy values leading to (possibly)\
                                      multiple edges not yet implemented.")

        self.energy = Measurement(value=energy, 
                                            _unit = units.gev)
        if self.energy.unit != self.energy_unit:
            self.energy.to(self.energy_unit)

        self.mass_absorption = Measurement(value = mass_absorption, 
                                           _unit = units.cm**2/units.g)
        if self.mass_absorption.unit != self.length_unit**2 / self.mass_unit:
            self.mass_absorption.to(self.length_unit**2 / self.mass_unit)

        self.calculate_step()

    @property
    def mass_unit(self):
        return self._mass_unit

    @mass_unit.setter
    def mass_unit(self, val:Unit|str):
        if isinstance(val, str):
            self._mass_unit = getattr(units, val)
        else:
            self._mass_unit = val
        # convert all mass-dependent values:
        self.mass.to(self.mass_unit)
        self.density.to(self.mass_unit/self.length_unit**3)
        self.surface_density.to(self.mass_unit/self.length_unit**2)
        self.mass_absorption.to(self.length_unit**2 / self.mass_unit)

        return self._mass_unit

    @property
    def length_unit(self):
        return self._length_unit

    @length_unit.setter
    def length_unit(self, val:Unit|str):
        if isinstance(val, Unit):
            self._length_unit = val
        else:
            self._length_unit = getattr(units, val)

        # convert all length-dependent values:
        self.thickness.to(self.length_unit)
        self.area.to(self.length_unit**2)
        self.density.to(self.mass_unit/self.length_unit**3)
        self.surface_density.to(self.mass_unit/self.length_unit**2)
        self.mass_absorption.to(self.length_unit**2 / self.mass_unit)
    
        return self._length_unit

    @property
    def energy_unit(self):
        return self._energy_unit

    @energy_unit.setter
    def energy_unit(self, val:str|Unit):
        if isinstance(val, Unit):
            self._energy_unit = val
        else:
            self._energy_unit = getattr(units, val)

        self.edge_energy.to(self.energy_unit)
        self.energy.to(self.energy_unit)
        

    def calculate_step(self):
        if isinstance(self.edge_energy.value, (int, float)):
            edge, bounds, pre = get_bounds(self.energy.value,
                                           self.mass_absorption.value,
                                           self.edge_energy.value)

            mass_abs_step, mass_abs_max, mass_abs_min = fit_bkg(self.energy.value,
                                                                self.mass_absorption.value,
                                                                bounds,
                                                                pre,
                                                                edge,
                                                                include=True)
            
            self.mass_abs_step = Measurement(value=mass_abs_step, _unit=self.mass_absorption.unit)
            self.mass_abs_max = Measurement(value=mass_abs_max, _unit=self.mass_absorption.unit)
            self.mass_abs_min = Measurement(value=mass_abs_min, _unit=self.mass_absorption.unit)

        else:
            raise NotImplementedError("Energy ranges not yet implemented.")

    def calculate_density(self):
        if self.density.value is not None:
            return
    
        if self.thickness.value is not None and self.area.value is not None:
            self.density = self.mass / (self.thickness * self.area)
        else:
            print("Not able to calculate density.")

    def calculate_thickness(self):
        """
        Calculate thickness for a foil.
        """
        self.calculate_density()
        # need density, mu_t for this.
        if self.density.value is None:
            raise ValueError("Density must not be None.")
        # make into a Measurement...
        self.linear_absorption = self.mass_absorption * self.density
        self.linear_step = self.mass_abs_step * self.density
        self.linear_abs_max = self.mass_abs_max * self.density
        self.linear_abs_min = self.mass_abs_min * self.density

        self.thickness.value = self.mu_total / self.linear_abs_max.value

        self.absorbance_step = self.linear_step * self.thickness

    def calculate_mass(self):
        """
        Calculate sample mass for a pellet.
        """
        # be careful all units are the same!
        self.mass = self.area * self.thickness * self.density