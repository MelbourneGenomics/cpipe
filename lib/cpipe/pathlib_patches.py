import pathlib
import shutil
import os

def path_copy(self, target):
    shutil.copy(str(self), str(target))

def path_rmtree(self):
    shutil.rmtree(str(self))

def path_symlink_from(self, link):
    os.symlink(str(self), str(link))

pathlib.Path.copy = path_copy
pathlib.Path.rmtree = path_rmtree
pathlib.Path.symlink_from = path_symlink_from
