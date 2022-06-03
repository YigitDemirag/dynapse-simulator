#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alpha Renner, Yigit Demirag, Ioan Fodorut, Giacomo Indiveri

Modified from https://code.ini.uzh.ch/ncs/teili
"""

def dynapse_eq():
    return{'model': '''
                    # Neuronal dynamics #########################################

                    # Differential equations
                    dIsoma_mem/dt = (((Isoma_dpi_g_shunt / Isoma_dpi_tau_shunt) *\
                    (Iin_clip + Isoma_pfb_shunt - Ishunt_clip - Isoma_ahp_shunt)) -\
                    Isoma_dpi_g_shunt - ((1 + ((Ishunt_clip + Isoma_ahp_shunt - Isoma_pfb_shunt) /\
                    Isoma_dpi_tau_shunt)) * Isoma_mem)) / (tau_soma *\
                    ((Isoma_dpi_g_shunt / Isoma_mem_clip) + 1)) : amp (unless refractory)

                    dIsoma_ahp/dt = (- Isoma_ahp_g_shunt - Isoma_ahp + 2*I0*(Isoma_ahp<=I0)) /\
                    (tau_soma_ahp * (1 + (Isoma_ahp_g_shunt / Isoma_ahp))) : amp # Adaptation current

                    Isoma_pfb = Isoma_pfb_gain/(1+exp(-(Isoma_mem-Isoma_pfb_th)/Isoma_pfb_norm)) : amp  # Positive feedback current
                    
                    # The *_clip and *_shunt currents are needed to prevent current from going
                    # below I0, since negative currents are not possible on chips
                    Isoma_mem_clip = clip(Isoma_mem, I0, 1*amp) : amp
                    Iin_clip = clip(Inmda_dp + Iampa - Igaba_b + Isoma_const, I0, 1*amp) : amp
                    Ishunt_clip = clip(Igaba_a, I0, Isoma_mem) : amp

                    Isoma_dpi_g_shunt = Isoma_dpi_g*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp           # Shunt g current if Imem goes to I0
                    Isoma_dpi_tau_shunt = Isoma_dpi_tau*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp       # Shunt tau current if Imem goes to I0
                    Isoma_pfb_shunt = Isoma_pfb*(Isoma_mem>I0) + I0*(Isoma_mem<=I0)    : amp
                    Isoma_ahp_shunt = Isoma_ahp*(Isoma_mem>I0) + I0*(Isoma_mem<=I0)  : amp
                    Isoma_ahp_g_shunt = Isoma_ahp_g*(Isoma_ahp>I0) + I0*(Isoma_ahp<=I0) : amp           # Shunt g current if Iahp goes to I0
                    Isoma_ahp_tau_shunt = Isoma_ahp_tau*(Isoma_ahp>I0) + I0*(Isoma_ahp<=I0) : amp       # Shunt g current if Iahp goes to I0
                  
                    tau_soma_ahp = (Csoma_ahp * Ut) / (kappa * Isoma_ahp_tau_shunt) : second            # Time constant of adaptation
                    tau_soma = (Csoma_mem * Ut) / (kappa * Isoma_dpi_tau_shunt) : second                # Membrane time constant                    
                    kappa = (kn + kp) / 2 : 1

                    # Substrate constants
                    kn              : 1     (shared, constant)                  # Subthreshold slope factor for nFETs
                    kp              : 1     (shared, constant)                  # Subthreshold slope factor for pFETs
                    Ut              : volt  (shared, constant)                  # Thermal voltage
                    I0              : amp   (shared, constant)                  # Dark current
                    alpha_soma      : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau
                    alpha_ahp       : 1     (shared, constant)                  # Scaling factor equal to Ig/Itau 
                    alpha_nmda      : 1     (constant)                          # Scaling factor equal to Ig/Itau 
                    alpha_ampa      : 1     (constant)                          # Scaling factor equal to Ig/Itau 
                    alpha_gaba_a    : 1     (constant)                          # Scaling factor equal to Ig/Itau 
                    alpha_gaba_b    : 1     (constant)                          # Scaling factor equal to Ig/Itau 

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
                    dInmda/dt = (- Inmda - Inmda_g_shunt + 2*I0*(Inmda <= I0))/(tau_nmda * ((Inmda_g_shunt / Inmda) + 1)) : amp

                    Inmda_dp = Inmda/(1 + Inmda_thr / Isoma_mem_clip) : amp             # Voltage gating differential pair block 

                    Inmda_tau_shunt = Inmda_tau*(Inmda>I0) + I0*(Inmda<=I0) : amp       # Shunt tau current if Inmda goes to I0
                    Inmda_g_shunt = Inmda_g*(Inmda>I0) + I0*(Inmda<=I0) : amp           # Shunt g current if Inmda goes to I0         
                    
                    Cnmda : farad (constant)                                            # Synapse's capacitance
                    Inmda_tau : amp (constant)                                          # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Inmda_thr : amp (constant)                                          # NMDA voltage-gating threshold
                    Inmda_w0 : amp (constant)                                           # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    Inmda_g = alpha_nmda * Inmda_tau : amp                              # NMDA synapse gain term expressed in terms of its tau current
                    tau_nmda = Cnmda * Ut /(kappa * Inmda_tau_shunt) : second           # Synaptic time-constant

                    # AMPA #######################################################
                    dIampa/dt = (- Iampa - Iampa_g_shunt + 2*I0*(Iampa <= I0))/(tau_ampa * ((Iampa_g_shunt / Iampa) + 1)) : amp
                    
                    Iampa_tau_shunt = Iampa_tau*(Iampa>I0) + I0*(Iampa<=I0) : amp       # Shunt tau current if Iampa goes to I0
                    Iampa_g_shunt = Iampa_g*(Iampa>I0) + I0*(Iampa<=I0) : amp           # Shunt g current if Iampa goes to I0

                    Campa : farad (constant)                                            # Synapse's capacitance
                    Iampa_tau : amp (constant)                                          # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Iampa_w0 : amp (constant)                                           # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    Iampa_g = alpha_ampa * Iampa_tau : amp                              # AMPA synapse gain expressed in terms of its tau current
                    tau_ampa = Campa * Ut / (kappa * Iampa_tau_shunt) : second          # Synaptic time-constant

                    # GABA B - inh #######################################################
                    # the ihn synapse does not actually decrease Imem, it just
                    # decreases the input current from other synapses
                    dIgaba_b/dt = (- Igaba_b - Igaba_b_g_shunt + 2*I0*(Igaba_b <= I0))/(tau_gaba_b * ((Igaba_b_g_shunt / Igaba_b) + 1)) : amp

                    Igaba_b_tau_shunt = Igaba_b_tau*(Igaba_b > I0) + I0*(Igaba_b <= I0) : amp   # Shunt tau current if Igaba_b goes to I0
                    Igaba_b_g_shunt = Igaba_b_g*(Igaba_b>I0) + I0*(Igaba_b<=I0) : amp           # Shunt g current if Igaba_b goes to I0

                    Cgaba_b : farad (constant)                                          # Synapse's capacitance
                    Igaba_b_tau : amp (constant)                                        # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_b_w0 : amp (constant)                                         # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    Igaba_b_g = alpha_gaba_b * Igaba_b_tau : amp                        # GABA A synapse gain expressed in terms of its tau current
                    tau_gaba_b = Cgaba_b * Ut / (kappa * Igaba_b_tau_shunt) : second    # Synaptic time-constant

                    # GABA A - shunt #####################################################
                    dIgaba_a/dt =(- Igaba_a - Igaba_a_g_shunt + 2*I0*(Igaba_a <= I0))/(tau_gaba_a * ((Igaba_a_g_shunt / Igaba_a) + 1)) : amp
                    
                    Igaba_a_tau_shunt = Igaba_a_tau*(Igaba_a > I0) + I0*(Igaba_a <= I0) : amp   # Shunt tau current if Iampa goes to I0
                    Igaba_a_g_shunt = Igaba_a_g*(Igaba_a>I0) + I0*(Igaba_a<=I0) : amp           # Shunt g current if Igaba_a goes to I0

                    Cgaba_a : farad (constant)                                          # Synapse's capacitance
                    Igaba_a_tau : amp (constant)                                        # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_a_w0 : amp (constant)                                         # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    Igaba_a_g = alpha_gaba_a * Igaba_a_tau : amp                        # GABA A synapse gain expressed in terms of its tau current
                    tau_gaba_a = Cgaba_a * Ut / (kappa * Igaba_a_tau_shunt) : second    # Synaptic time-constant
                    ''',
           'threshold': '''Isoma_mem > Isoma_th''',
           'reset': '''
                    Isoma_ahp += Isoma_ahp_w
                    Isoma_mem = Isoma_reset
                    ''',
           'refractory': 'soma_refP',
           'method': 'euler'}

# Warning: do NOT write comments in the on_pre / on_post equations,
# to prevent problems with GPU code generation.

def dynapse_nmda_syn_eq():  # SLOW_EXC
    """This function returns the slow excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Inmda_post += Inmda_w0_post * weight * alpha_nmda_post
                    """, # On pre-synaptic spike adds current to state variable of DPI synapse.
           'on_post': """ """,
           'method': 'euler'}


def dynapse_ampa_syn_eq():  # FAST_EXC
    """This function returns the fast excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Iampa_post += Iampa_w0_post * weight * alpha_ampa_post
                    """, # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_b_syn_eq():  # SLOW_INH
    """This function returns the subtractive inhibitory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Igaba_b_post += Igaba_b_w0_post * weight * alpha_gaba_b_post
                    """, # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_a_syn_eq():  # FAST_INH
    """This function returns the shunting synapse (a mixture of subtractive and divisive) equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Igaba_a_post += Igaba_a_w0_post*weight * alpha_gaba_a_post
                    """, # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}
