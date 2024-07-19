"""
This module analyse sections of beams by using the boundary element method
It uses mainly curves as boundary to compute the elements.
"""

from .curve import Curve
from .dataio import JsonIO
from .field import plot_field, plot_section
from .geometry import ConnectedGeometry
from .material import Isotropic, Material
from .section import HomogeneousSection

__version__ = "0.4.0"

if __name__ == "__main__":
    pass
