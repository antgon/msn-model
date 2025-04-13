from pathlib import Path

_root = Path(__path__[0])

files = {
    'mechanisms': _root.joinpath('mechanisms'),
    'morphology': {
        'dmsn': _root.joinpath(
            'morphology', 'WT-dMSN_P270-20_1.02_SGA1-m24.swc'),
        'imsn': _root.joinpath(
            'morphology', 'WT-iMSN_P270-09_1.01_SGA2-m1.swc')},
    'conductances': _root.joinpath('parameters', 'conductances.csv'),
    'params': {
        'dmsn': _root.joinpath(
            'parameters', 'D1_71bestFit_updRheob.pkl'),
        'imsn': _root.joinpath(
            'parameters', 'D2_34bestFit_updRheob.pkl')}}
            
# These imports require the `files` dictionary defined above so should
# not be moved to the top.
from .params import ModelParameters as ModelParameters
from . import cell as cell
