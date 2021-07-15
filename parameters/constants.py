# -*- coding: utf-8 -*-
# @Author: mmilde
# @Date:   2018-02-10 09:09:20
# @Last Modified by:   mmilde
# @Last Modified time: 2018-02-10 10:08:10
# Taken from https://code.ini.uzh.ch/ncs/teili

"""
This file contains shared constants for synapses and neuron models.
These constants were empirically estimated. For further information
have a look at Liu, Kramer, Indiveri, Delbruck and Douglas Analog VLSI: Circuits and Principles
"""

from brian2 import mV, pA

kappa_n = 0.75  # Subthreshold slope factor (n-type transistor)
kappa_p = 0.66  # Subthreshold slope factor (p-type transistor)
Ut = 25. * mV  # Thermal voltage
Io = 0.5 * pA  # Dark current
