#! /usr/bin/env python3

import numpy as np
import pandas as pd

from pymatgen.analysis.transition_state import NEBAnalysis
# http://pymatgen.org/_modules/pymatgen/analysis/transition_state.html

neb = NEBAnalysis.from_dir( '.' )

output = np.vstack( ( neb.energies, neb.r, neb.forces ) ).T
df = pd.DataFrame( output, columns=[ 'energy', 'distance', 'force' ] )
df['relative energy'] = df.energy - df.energy.min()
print( df )
