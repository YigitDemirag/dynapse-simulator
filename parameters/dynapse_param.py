#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Modified from https://code.ini.uzh.ch/ncs/teili

"""
This dictionary contains the default parameters for the
Differential Pair Integrator (DPI) neuron including synapses as implemented on
the DYNAP-SE chip, developed at the Institute of Neuroinformatics.
Circuit and equations were published in Chicca et al. 2014
These values were estimated empircally based on the example in the tutorial.
These parameters may serve you as a starting point for a given experiment.

Variable Naming Convention: 
"""

from brian2 import pF, ms, pA, nA, mV
from parameters import constants


dynapse_param = {
    #  SUBSTRATE  #########################################################################################
    "kn": constants.kappa_n,                # Subthreshold slope factor, do not change
    "kp": constants.kappa_p,                # Subthreshold slope factor, do not change
    "Ut": constants.Ut,                     # Thermal voltage, do not change
    # Dark current, can be from 2x to 10x default Io if necessary
    "I0": constants.I0,
    ##################

    "alpha_soma": 4,                        # Scaling factor equal to Ig/Itau
    "alpha_ahp": 4,                         # Scaling factor equal to Ig/Itau
    "alpha_nmda": 4,                        # Scaling factor equal to Ig/Itau
    "alpha_ampa": 4,                        # Scaling factor equal to Ig/Itau
    "alpha_gaba_a": 4,                      # Scaling factor equal to Ig/Itau
    "alpha_gaba_b": 4,                      # Scaling factor equal to Ig/Itau

    #  Neuron parameters  ###############
    #  SOMA  ##############################################################################################
    # Membrane capacitance, fixed at layout time (see chip for details)
    "Csoma_mem": 2 * pF,
    # Membrane time constant current, the time constant is inversely proportional to Itau
    "Isoma_dpi_tau": 5 * constants.I0,
    # Spiking threshold current, depends on layout (see chip for details)
    "Isoma_th": 2000 * constants.I0,
    "Isoma_reset": 1.2 * constants.I0,      # Reset current after spike generation
    # Initialize constant current injection to Io
    "Isoma_const": constants.I0,
    # New Default according to Refractory Period experiments: Refractory period, limits maximum firing rate at 200Hz
    "soma_refP": 5. * ms,

    #  ADAPTATION  ########################################################################################
    "Csoma_ahp": 4 * pF,                    # Spike-frequency adaptation capacitance
    # Spike-frequency adaptation time constant current
    "Isoma_ahp_tau": 1 * constants.I0,
    # Spike-frequency adaptation weight current
    "Isoma_ahp_w": 1 * constants.I0,

    #  POSITIVE FEEDBACK ##################################################################################
    # Feedback gain current, fitting parameter
    "Isoma_pfb_gain": 100 * constants.I0,
    # Feedback threshold current, typically a fraction of Isoma_th
    "Isoma_pfb_th": 1000 * constants.I0,
    # Feedback normalization current, fitting parameter
    "Isoma_pfb_norm": 20 * constants.I0,
    ##################

    # Synapse parameters ################
    # SLOW_EXC, NMDA ########################################################################################
    # Synaptic capacitance, fixed at layout time (see chip for details)
    'Cnmda': 2 * pF,
    # New Default according to guidance from Giacomo: Synaptic time constant current, the time constant is inversely proportional to I_tau
    'Inmda_tau': 2 * constants.I0,
    # Base synaptic weight current which can be scaled by the .weight parameter
    'Inmda_w0': 100 * constants.I0,
    'Inmda': constants.I0,                  # Output current initial value
    'Vnmda': 10 * mV,                       # Voltage NMDA DPI slow

    #FAST_EXC, AMPA ########################################################################################
    # Synaptic capacitance, fixed at layout time (see chip for details)
    'Campa': 2 * pF,
    # Synaptic time constant current, the time constant is inversely proportional to I_tau
    'Iampa_tau': 20 * constants.I0,
    # Base synaptic weight current which can be scaled by the .weight parameter
    'Iampa_w0': 100 * constants.I0,
    'Iampa': constants.I0,                  # Output current initial value

    #INH, SLOW_INH, GABA_B, subtractive ##################################################################
    # Synaptic capacitance, fixed at layout time (see chip for details)
    'Cgaba_b': 2 * pF,
    # Synaptic time constant current, the time constant is inversely proportional to I_tau
    'Igaba_b_tau': 5 * constants.I0,
    # Base synaptic weight current which can be scaled by the .weight parameter
    'Igaba_b_w0': 100 * constants.I0,
    'Igaba_b': constants.I0,                # Output current initial value

    #FAST_INH, GABA_A, shunting, a mixture of subtractive and divisive ############################################
    # Synaptic capacitance, fixed at layout time (see chip for details)
    'Cgaba_a': 2 * pF,
    # Synaptic time constant current, the time constant is inversely proportional to I_tau
    'Igaba_a_tau': 5 * constants.I0,
    # Base synaptic weight current which can be scaled by the .weight parameter
    'Igaba_a_w0': 100 * constants.I0,
    'Igaba_a': constants.I0                 # Output current initial value
    ##################
}
