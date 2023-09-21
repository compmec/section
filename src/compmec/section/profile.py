from typing import Optional, Tuple

from compmec.shape.shape import DefinedShape

from compmec.section.material import Material


class Profile:
    def __init__(self, shapes: Tuple[DefinedShape], materials: Tuple[Material]):
        for shape in shapes:
            assert isinstance(shape, DefinedShape)
        for material in materials:
            assert isinstance(material, Material)
        self.__shapes = tuple(shapes)
        self.__materials = tuple(materials)

    def __iter__(self) -> Tuple[DefinedShape, Material]:
        for shape, material in zip(self.__shapes, self.__materials):
            yield (shape, material)
