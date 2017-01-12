import os
import pathlib
from .paths import DESIGNS


class Design:
    def __init__(self, name):
        self.path = DESIGNS / name

    @property
    def name(self):
        return str(self.path.relative_to(DESIGNS))

    @classmethod
    def find_by_name(cls, name: str) -> 'Design':
        """
        Returns a Design with the given name, or None if it does not exist
        :param name: A design name
        """
        found = [design for design in cls.list_all() if design.name == name]
        return found[0] if found else None

    @staticmethod
    def list_all():
        """
            Prints the name of all designs in the designs directory
        """

        # Find all directories that contain a samples.txt and add them to a list
        designs = []

        for file in DESIGNS.iterdir():
            if file.is_dir():
                designs.append(Design(file))

        # Sort by name
        return sorted(designs, key=lambda x: x.name)
