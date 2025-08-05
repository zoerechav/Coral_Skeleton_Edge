#packages needed for this work
import numpy as np
import h5py
import matplotlib
import pandas as pd
import scipy
import importlib.metadata as importlib_metadata  # Python 3.8+


def get_version(package_name):
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return "Not installed"

# Print versions of third-party packages
print("Library Versions:")
print(f"NumPy: {np.__version__}")
print(f"h5py: {h5py.__version__}")
print(f"Matplotlib: {matplotlib.__version__}")
print(f"Pandas: {pd.__version__}")
print(f"SciPy: {scipy.__version__}")

# Note: glob and csv are standard library modules
print("\nStandard Library Modules:")
print("glob: part of Python standard library")
print("csv: part of Python standard library")




##OUTPUT
# Library Versions:
# NumPy: 1.24.3
# h5py: 3.8.0
# Matplotlib: 3.7.1
# Pandas: 2.0.1
# SciPy: 1.10.1

# Standard Library Modules:
# glob: part of Python standard library
# csv: part of Python standard library

