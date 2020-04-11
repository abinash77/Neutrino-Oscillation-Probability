# -*- coding: utf-8 -*-
r"""Run the LIV 3nu example shown in README.md.

Runs the three-neutrino example of oscillations with LIV shown in
README.md

References
----------

.. [1] Mauricio Bustamante, "Exact neutrino oscillation probabilities:
   a fast general-purpose computation method for two and three neutrino
   flavors", arXiv:1904.XXXXX.

Created: 2019/04/30 00:50
Last modified: 2019/04/30 00:50
"""


from __future__ import print_function

__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import sys
sys.path.append('../src')

import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = \
    hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_NO_BF,
                                                                S23_NO_BF,
                                                                S13_NO_BF,
                                                                DCP_NO_BF,
                                                                D21_NO_BF,
                                                                D31_NO_BF)

# The values of the LIV parameters (SXI12, SXI23, SXI13, DXICP, B1, B2,
# B3, LAMBDA) are read from globaldefs
h_liv = hamiltonians3nu.hamiltonian_3nu_liv(h_vacuum_energy_indep, energy,
                                            SXI12, SXI23, SXI13, DXICP,
                                            B1, B2, B3, LAMBDA)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
    oscprob3nu.probabilities_3nu(h_liv, baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
