import os
import pathlib
import typing

from .paths import BATCHES

class Batch:
    def __init__(self, directory):
        self.path = pathlib.Path(directory).resolve()

    @property
    def name(self):
        return str(self.path.relative_to(BATCHES))

    @classmethod
    def find_by_name(cls, name: str) -> 'Batch':
        """
        Returns a Batch with the given name, or None if it does not exist
        :param name: A batch name
        """
        found = [batch for batch in cls.list_all() if batch.name == name]
        return found[0] if found else None

    @staticmethod
    def list_all() -> typing.Iterable['Batch']:
        """
            Prints the name of all batches that contain a samples.txt file
        """

        # Find all directories that contain a samples.txt and add them to a list
        batches = []

        for root, dirs, files in os.walk(str(BATCHES)):
            if 'samples.txt' in files:
                batches.append(Batch(root))

        # Sort by batch name
        return sorted(batches, key=lambda x: x.name)
