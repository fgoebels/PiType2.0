import sys
import os
from cx_Freeze import setup, Executable, hooks

import scipy

scipy_path = os.path.dirname(scipy.__file__)

build_exe_options = {"packages": ['numpy.core._methods','numpy.lib.format'],"excludes": ['scipy'],"include_files":[scipy_path]}

# Dependencies are automatically detected, but it might need fine tuning.


setup(  name = "guifoo",
        version = "0.1",
        description = "My GUI application!",
        options = {"build_exe": build_exe_options},
        executables = [Executable("Pitype2.py")])