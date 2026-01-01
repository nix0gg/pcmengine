# Physics/__init__.py

from . import conversion_factors
from . import dimensional_formulas
from . import gravitation
from . import laws_of_motion
from . import mechanical_properties_solids
from . import motion_in_one_dimension
from . import motion_in_plane
from . import physics_constants
from . import system_of_particles
from . import work_power_energy
from . import thermal_properties_of_matter
from . import mechanical_properties_fluids

__all__ = [
    # Core physics areas
    *conversion_factors.__all__,
    *dimensional_formulas.__all__,
    *gravitation.__all__,
    *laws_of_motion.__all__,
    *mechanical_properties_solids.__all__,
    *motion_in_one_dimension.__all__,
    *motion_in_plane.__all__,
    *physics_constants.__all__,
    *system_of_particles.__all__,
    *work_power_energy.__all__,
    *thermal_properties_of_matter.__all__,
    *mechanical_properties_fluids.__all__,
]
