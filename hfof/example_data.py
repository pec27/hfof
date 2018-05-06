import numpy as np
from os import path

_data = None

def get_pos():
    """
    Load cosmo simulation data (unit cube, (32768,3) array)
    """
    global _data
    if _data is not None:
        return _data

    name = path.join(path.dirname(path.abspath(__file__)), '../tests/cosmo32768.dat')
    if not path.exists(name):
        raise Exception('File '+str(name)+' does not exist. Maybe I dont understand your OS?')

    pos = np.reshape(np.fromfile(name, np.uint16, 32768*3), (32768,3)).astype(np.float64)
    pos *= 1.0/65536

    _data = pos
    return _data
