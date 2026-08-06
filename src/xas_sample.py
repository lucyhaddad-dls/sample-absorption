# set all input params as properties:
from .photoabsorption import CS_Photo_Formula
from xraylib import EdgeEnergy, SymbolToAtomicNumber
from dataclasses import dataclass
import numpy as np
from .formulas import _edges
from .edge_step import get_edge_vals, get_bounds, fit_bkg

# make units for everything:
from pint import UnitRegistry, Unit

import pint
units = pint.UnitRegistry(system='mks')

units.default_system

units.define("eV = [energy] = EV = ev")
units.define("GeV = 1e3 eV = GEV = gev")

# make a Measurement subclass that is dependent on a parent's
# value?

@dataclass
class Measurement:
    value: int | float | np.ndarray | None
    _unit: Unit

    @property
    def unit(self)->None:
        """
        Set unit value.
        """
        return self._unit

    @unit.setter
    def unit(self, value:str|Unit):
        if isinstance(value, Unit):
            self._unit = value
        elif isinstance(value, str):
            value = getattr(units, value)
            self._unit = value

    def to(self, unit:str|Unit):
        """
        Convert value to new unit.
        """
        if isinstance(unit, Unit):
            unit = str(unit)

        new = (self.value * self.unit).to(unit)
        self.value = new._magnitude
        self.unit = getattr(units, unit)

    def __mul__(self, other):
        return Measurement(value = self.value * other.value,
                          unit = self.unit*other.unit)
    
    def __add__(self, other):
        if self.unit != other.unit:
            try:
                other.to(self.unit)
            except:
                raise ValueError(f"Units {self.unit} {other.unit} \
                                 not able to be summed.")

        return Measurement(value = self.value + other.value,
                           unit = self.unit)

    def __sub__(self, other):
        if self.unit != other.unit:
            try:
                other.to(self.unit)
            except:
                raise ValueError(f"Units {self.unit} {other.unit} \
                                    not able to be subtracted.")
        return Measurement(value = self.value - other.value,
                           unit = self.unit)

    def __truediv__(self, other):
        if self.unit != other.unit:
            try:
                # this will return a scalar - no longer a measurement.
                other.to(self.unit)
                return self.value/other.value

            except:
                return Measurement(value = self.value/other.value,
                                   unit = self.unit/other.unit)

    def __pow__(self, exp):
        return Measurement(value= self.value**exp, 
                           unit = self.unit**exp)
    
class XRaySample:
    def __init__(self,
                formula:str,
                absorber:str|None = None,
                edge:str|None = None,
                density:float|int|None = None,
                surface_density:float|int|None = None,
                mass:float|int|None = None,
                area:float|int|None = None,
                volume:float|int|None = None,
                thickness:float|int|None = None,
                mu_total:float|int|None = 2.6,
                mass_unit:Unit|str|None = None,
                length_unit:Unit|str|None = None):
        """
        Assumed units are: mass [g], length [cm], energy [GeV]. To change unit:
        `sample.thickness.unit = "mm"`
        """

        self.formula = formula
        self.absorber = absorber
        self.edge = edge
        self.mu_total = mu_total

        if isinstance(mass_unit, Unit):
            self._mass_unit = mass_unit
        elif isinstance(mass_unit, str):
            self._mass_unit = getattr(units, mass_unit)
        else:
            self._mass_unit = units.g

        if isinstance(length_unit, Unit):
            self._length_unit = length_unit
        elif isinstance(length_unit, str):
            self._length_unit = getattr(units, length_unit)
        else:
            self._length_unit = units.cm
        
        # probs. want to make these properties.... because if one unit changes !
        # all other units are based on these two.
        self.mass = Measurement(value = mass, _unit=self.mass_unit)
        self.thickness = Measurement(value = thickness, _unit = self.length_unit)
        self.volume = Measurement(value = volume, _unit = self.length_unit**3)
        self.density = Measurement(value = density, _unit = self.mass_unit/self.length_unit**3)
        self.area = Measurement(value = area, _unit = self.length_unit**2)
        self.surface_density = Measurement(value = surface_density, _unit = self.mass_unit/self.length_unit**2)
        

        if self.absorber != None and self.edge != None:
            out = CS_Photo_Formula(self.formula,
                                   self.absorber,
                                   self.edge)
            self.edge_energy = Measurement(value= EdgeEnergy(SymbolToAtomicNumber(self.absorber),
                                          _edges.index(self.edge.upper())),
                                          _unit = units.GEV)
        else:
            raise NotImplementedError("Energy values leading to (possibly)\
                                      multiple edges not yet implemented.")

        self.energy = Measurement(value=out["energy"], _unit=units.GEV)


        self.mass_absorption = Measurement(value=out["mu_m"], _unit=units.cm**2/units.g)

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
        self.volume.to(self.length_unit**3)
        self.density.to(self.mass_unit/self.length_unit**3)
        self.surface_density.to(self.mass_unit/self.length_unit**2)
    
        return self._length_unit


    def calculate_step(self):
        if isinstance(self.edge_energy.value, (int, float)):
            edge, bounds, pre = get_bounds(self.energy,
                                           self.mass_absorption,
                                           self.edge_energy)

            mass_abs_step, mass_abs_max, mass_abs_min = fit_bkg(self.energy,
                                                                self.mass_absorption,
                                                                bounds,
                                                                pre,
                                                                edge,
                                                                include=True)
            
            self.mass_abs_step = Measurement(value=mass_abs_step, _unit=self.mass_absorption.unit)
            self.mass_abs_max = Measurement(value=mass_abs_max, _unit=self.mass_absorption.unit)
            self.mass_abs_min = Measurement(value=mass_abs_min, __name__unit=self.mass_absorption.unit)

        else:
            raise NotImplementedError("Energy ranges not yet implemented.")