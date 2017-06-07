import sys
import os
from cx_Freeze import setup, Executable, hooks

import scipy

scipy_path = os.path.dirname(scipy.__file__)

build_exe_options = {"packages": ['numpy.core._methods','numpy.lib.format','numpy.matlib'],"excludes": ['scipy']}

# Dependencies are automatically detected, but it might need fine tuning.


setup(  name = "Pitype Test",
        version = "0.1",
        description = "My GUI application!",
        options = {"build_exe": build_exe_options},
        executables = [Executable("Pitype2.py")])