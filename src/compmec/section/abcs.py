"""
File that contains base classes
"""

from abc import ABC, abstractmethod
from typing import Any, Optional, Tuple, Union


class Tracker(ABC):
    """
    Parent class to track all the instances
    """

    instances = {}

    @classmethod
    def clear(cls, keys: Optional[Tuple[int]] = None):
        """
        Remove node labels
        """
        if keys is None:
            cls.instances.clear()
            return
        for key in keys:
            if key in cls.instances:
                cls.instances.pop(key)

    @classmethod
    @abstractmethod
    def _next_available_key(cls) -> Any:
        """
        Protected method that gives next key
        """
        raise NotImplementedError

    def _get_key(self) -> Any:
        """
        Protected method that gives the key of the instance by iterating
        on the instances dictionary

        Parameters
        ----------

        :return: The instances key from the dictionary
        :rtype: Any

        """
        for key, instance in self.instances.items():
            if id(instance) == id(self):
                return key
        return None

    def _set_key(self, new_key: Any):
        """
        Protected method that sets a new key for the instance.
        Iterates over the dictionary and changes the name.

        If `None` is given, it will attribute the next available key

        Parameters
        ----------

        :param new_key: The new instance key"
        :type new_key: Any

        """
        cur_key = self._get_key()
        if new_key is None:
            new_key = self._next_available_key()
        elif cur_key == new_key:
            return
        elif new_key in self.instances:
            msg = f"Key {new_key} is already used"
            raise ValueError(msg)
        elif cur_key is not None:
            self.instances.pop(cur_key)
        self.instances[new_key] = self


class LabeledTracker(Tracker):
    """
    Labeled Base Tracker, to track instances of Curve and Node
    """

    @classmethod
    def _next_available_key(cls) -> str:
        """
        Protected method that gives next key
        """
        index = 1
        while index in cls.instances:
            index += 1
        return index

    @property
    def label(self) -> int:
        """
        Gives the instance label

        :getter: Returns the instance's label
        :setter: Attribuates a new label for instance
        :type: str

        """
        return int(self._get_key())

    @label.setter
    def label(self, new_label: Union[int, None]):
        if new_label is not None:
            new_label = int(new_label)
        self._set_key(new_label)


class NamedTracker(Tracker):
    """
    Named Base Tracker, to track instances of Material, Geometry and Section
    """

    @classmethod
    def _next_available_key(cls) -> str:
        """
        Protected method that gives next key
        """
        index = 1
        while f"instance-{index}" in cls.instances:
            index += 1
        return f"instance-{index}"

    @property
    def name(self) -> str:
        """
        Gives the instance name

        :getter: Returns the instance's name
        :setter: Attribuates a new name for instance
        :type: str

        """
        return str(self._get_key())

    @name.setter
    def name(self, new_name: Union[str, None]):
        if new_name is not None:
            new_name = str(new_name)
        self._set_key(new_name)
