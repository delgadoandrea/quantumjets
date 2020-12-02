# Taken from main06.cc Pythia example. A study of event properties of LEP1 events

from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID

import numpy as np
import sys

# Pythia generator settings. Process selection e+e- collision.
# WeakSingleBoson:ffbar2gmZ : Produce Z bosons and interference with photon
# Consider only Z (PDGID=23) decays to light quarks

def generator(cut = 0.0):
    params = {"Beams:idA" : "11", "Beams:idB" : "-11","Beams:eCM": 91.1876, "WeakSingleBoson:ffbar2gmZ" : "on", "23:onMode": "off", "23:onIfAny": "1 2 3 4 5", "PDF:lepton": "off"}
    pythia = Pythia(params =params, random_state=1)

    # Consider only final-state particles
    selection = ((STATUS == 1) & ~HAS_END_VERTEX)

    # Array of pseudo-jets or particles to cluster
    particles = []
    array = []

    # Begin event loop. Generate events
    for event in pythia(events=1):
        array = event.all(selection)

    #For testing: filter out particles by energy
    particles = array[array["E"] > cut]

    return particles

