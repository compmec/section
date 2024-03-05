"""
This module analyse sections of beams by using the boundary element method
It uses mainly curves as boundary to compute the elements.
"""

from .curve import Curve
from .dataio import JsonIO
from .geometry import Geometry
from .material import Isotropic, Material
from .section import Section

__version__ = "0.4.0"

if __name__ == "__main__":
    pass
