#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alpha Renner, Yigit Demirag, Ioan Fodorut, Giacomo Indiveri

Modified from https://code.ini.uzh.ch/ncs/teili
"""
# neuron


def dynapse_eq():
    return{'model': '''
                    # Neuronal dynamics #########################################

                    # Differential equations
                    dIsoma_mem/dt = (((Isoma_dpi_g / Isoma_dpi_tau) *\
                    (Iin_clip + Isoma_pfb_clip - Igaba_a_clip)) -\
                    Isoma_dpi_g - ((1 + ((Igaba_a_clip - Isoma_pfb_clip) /\
                    Isoma_dpi_tau)) * Isoma_mem_clip)) / (tau_soma *\
                    (1+(Isoma_dpi_g/Isoma_mem_clip))) : amp (unless refractory)

                    dIsoma_ahp/dt = (- Isoma_ahp_g - Isoma_ahp_clip) / (tau_soma_ahp *\
                    (1+(Isoma_ahp_g/Isoma_ahp_clip))) : amp # adaptation current

                    Isoma_pfb = Isoma_pfb_gain / (1 + exp(-(Isoma_mem_clip - Isoma_pfb_th) / Isoma_pfb_norm)) : amp       # Positive feedback current

                    Iin_clip=clip(Inmda_dp + Iampa_clip - Igaba_b_clip - Isoma_ahp_clip +\
                    Isoma_const, I0,1*amp) : amp
                    Isoma_mem_clip=clip(Isoma_mem,I0,1*amp) : amp
                    Isoma_pfb_clip=clip(Isoma_pfb,I0,1*amp) : amp
                    Isoma_ahp_clip=clip(Isoma_ahp,I0,1*amp) : amp

                    tau_soma_ahp = (Csoma_ahp * Ut) / (kappa * Isoma_ahp_tau) : second              # Time constant of adaptation
                    tau_soma = (Csoma_mem * Ut) / (kappa * Isoma_dpi_tau) : second             # Membrane time constant                    
                    kappa = (kn + kp) / 2 : 1
                    Vsoma_mem = Ut / kappa * log(Isoma_mem / I0) : volt           # Membrane voltage

                    # Substrate constants
                    kn              : 1     (shared, constant)                  # Subthreshold slope factor for nFETs
                    kp              : 1     (shared, constant)                  # Subthreshold slope factor for pFETs
                    Ut              : volt  (shared, constant)                  # Thermal voltage
                    I0              : amp   (shared, constant)                  # Dark current
                    alpha_soma      : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau
                    alpha_ahp       : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 
                    alpha_nmda      : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 
                    alpha_ampa      : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 
                    alpha_gaba_a    : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 
                    alpha_gaba_b    : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 

                    # Soma constants
                    Csoma_mem       : farad (shared, constant)                  # Membrane capacitance
                    Isoma_dpi_tau   : amp   (shared, constant)                  # Leakage current
                    Isoma_dpi_g = alpha_soma * Isoma_dpi_tau : amp              # Soma DPI gain current 
                    Isoma_th        : amp   (shared, constant)                  # Spiking threshold
                    Isoma_reset     : amp   (shared, constant)                  # Reset current
                    Isoma_const     : amp   (constant)                          # Additional input current similar to constant current injection
                    soma_refP       : second (shared, constant)                 # Refractory period

                    # Adaptation constants
                    Csoma_ahp       : farad (shared, constant)                  # Spike-frequency adaptation capacitance
                    Isoma_ahp_tau   : amp   (shared, constant)                  # Leakage current for spike-frequency adaptation
                    Isoma_ahp_g = alpha_ahp * Isoma_ahp_tau  : amp              # AHP gain current
                    Isoma_ahp_w     : amp   (constant)                          # AHP jump height, on post

                    # Positive feedback constants
                    Isoma_pfb_gain  : amp   (shared, constant)                  # Positive feedback gain
                    Isoma_pfb_th    : amp   (shared, constant)                  # Positive feedback threshold (since it is a DPI circuit)
                    Isoma_pfb_norm  : amp   (shared, constant)                  # Positive feedback normalization current

                    # Synaptic dynamics #########################################

                    # NMDA #######################################################
                    dInmda/dt = (-Inmda_clip - Inmda_g)/(tau_nmda*\
                    (1+(Inmda_g/Inmda_clip))) : amp

                    Inmda_clip = clip(Inmda,I0,1*amp) : amp
                    Inmda_dp = Inmda_clip / (1 + Inmda_thr / Isoma_mem_clip) : amp # Output current of the DPI Slow synapse 

                    Vnmda : volt  (constant)                                    # NMDA threshold voltage
                    Inmda_tau : amp (constant)                                  # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Inmda_g = alpha_nmda * Inmda_tau : amp                      # Current flowing through ?? sets the DPI's threshold
                    Inmda_thr = I0 * exp(kappa * Vnmda / Ut) : amp
                    Inmda_w0 : amp (constant)                                   # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_nmda = Cnmda * Ut /(kappa * Inmda_tau) : second    # Synaptic time-constant
                    Cnmda : farad (constant)                                    # Synapse's capacitance

                    # AMPA #######################################################
                    dIampa/dt = (-Iampa_clip - Iampa_g)/(tau_ampa*\
                    (1+(Iampa_g/Iampa_clip))) : amp

                    Iampa_clip = clip(Iampa,I0,1*amp) : amp
                    Iampa_tau : amp (constant)                                  # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Iampa_g = alpha_ampa * Iampa_tau : amp                      # Current flowing through ?? sets the DPI's threshold
                    Iampa_w0 : amp (constant)                                   # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_ampa = Campa * Ut / (kappa * Iampa_tau) : second   # Synaptic time-constant
                    Campa          : farad (constant)                           # Synapse's capacitance

                    # GABA B - inh #######################################################
                    # the ihn synapse does not actually decrease Imem, it just
                    # decreases the input current from other synapses
                    dIgaba_b/dt = (-Igaba_b_clip - Igaba_b_g)/(tau_gaba_b*\
                    (1+(Igaba_b_g/Igaba_b_clip))) : amp

                    Igaba_b_clip = clip(Igaba_b,I0,1*amp) : amp
                    Igaba_b_tau      : amp (constant)                           # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_b_g = alpha_gaba_b * Igaba_b_tau : amp                # Current flowing through ?? sets the DPI's threshold
                    Igaba_b_w0 : amp (constant)                                 # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_gaba_b = Cgaba_b * Ut / (kappa * Igaba_b_tau) : second    # Synaptic time-constant
                    Cgaba_b          : farad (constant)                         # Synapse's capacitance

                    # GABA A - shunt #####################################################
                    dIgaba_a/dt = (-Igaba_a_clip - Igaba_a_g)/(tau_gaba_a*\
                    (1+(Igaba_a_g/Igaba_a_clip))) : amp

                    Igaba_a_clip = clip(Igaba_a,I0,1*amp) : amp
                    Igaba_a_tau     : amp (constant)                            # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_a_g = alpha_gaba_a * Igaba_a_tau : amp                # Current flowing through ?? sets the DPI's threshold
                    Igaba_a_w0 : amp (constant)                                 # Synaptic weight, to convert unitless weight to current
                    tau_gaba_a = Cgaba_a * Ut / (kappa * Igaba_a_tau) : second     # Synaptic time-constant
                    Cgaba_a         : farad (constant)                          # Synapse's capacitance
                    ''',
           'threshold': '''Isoma_mem > Isoma_th''',
           'reset': '''
                    Isoma_ahp += Isoma_ahp_w
                    Isoma_mem = Isoma_reset
                    ''',
           'refractory': 'soma_refP',
           'method': 'euler'}

# synapses
# Inmda_post += Inmda_w0*weight*Inmda_g_post/(Inmda_tau_post*((Inmda_g_post/Inmda_post)+1))


def dynapse_nmda_syn_eq():  # SLOW_EXC
    """This function returns the slow excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Inmda_post += Inmda_w0*weight
                    """,
           'on_post': """ """,
           'method': 'euler'}


def dynapse_ampa_syn_eq():  # FAST_EXC
    """This function returns the fast excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Iampa_post += Iampa_w0*weight
                    """,
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_b_syn_eq():  # SLOW_INH
    """This function returns the subtractive inhibitory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                     Igaba_b_post += Igaba_b_w0_post*weight
                     """,
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_a_syn_eq():  # FAST_INH
    """This function returns the shunting synapse (a mixture of subtractive and divisive) equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Igaba_a_post += Igaba_a_w0_post*weight
                    """,            # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}
