#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alpha Renner, Yigit Demirag

Modified from https://code.ini.uzh.ch/ncs/teili
"""
#neuron
def dynapse_eq():
    return{'model': '''
                    # Neuronal dynamics #########################################

                    # Differential equations
                    dIsoma_mem/dt = (((Isoma_dpi_g_clip / Isoma_dpi_tau_clip) * (Isoma_in_clip + Isoma_pfb_clip -\
                    Isoma_shunt_clip - Isoma_ahp_clip)) - Isoma_dpi_g_clip - ((1 + ((Isoma_shunt_clip +\
                    Isoma_ahp_clip - Isoma_pfb_clip) / Isoma_dpi_tau_clip)) * Isoma_mem)) / (tau_soma *\
                    ((Isoma_dpi_g_clip/(Isoma_mem + I0)) + 1)) : amp (unless refractory)

                    dIsoma_ahp/dt = (- Isoma_ahp_th_clip - Isoma_ahp + 2*I0*(Isoma_ahp<=I0)) / (tau_soma_ahp *\
                    (Isoma_ahp_th_clip / Isoma_ahp + 1)) : amp # adaptation current

                    # The *_clip currents are needed to prevent current from going
                    # below I0, since negative currents are not possible on chips
                    Isoma_dpi_tau_clip = Isoma_dpi_tau*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp
                    Isoma_dpi_g_clip = Isoma_dpi_g*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp
                    Isoma_in_clip = clip(Inmda_dp + Iampa - Igaba_b + Isoma_const, I0, 1*amp) : amp
                    Isoma_ahp_clip = Isoma_ahp*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp
                    Isoma_pfb_clip = Isoma_pfb*(Isoma_mem>I0) + 2*I0*(Isoma_mem<=I0) : amp
                    Isoma_ahp_th_clip = Isoma_ahp_th*(Isoma_ahp>I0) + I0*(Isoma_ahp<=I0) : amp
                    Isoma_shunt_clip = clip(Igaba_a, I0, Isoma_mem) : amp
                    Isoma_mem_clip = Isoma_mem*(Isoma_mem>I0) + I0*(Isoma_mem<=I0) : amp

                    Isoma_ahp_max = (Isoma_ahp_w / Isoma_ahp_tau) * Isoma_ahp_th_clip : amp         # Ratio of currents through diffpair and adaptation block
                    Isoma_pfb = Isoma_pfb_gain / (1 + exp(-(Isoma_mem - Isoma_pfb_th) / Isoma_pfb_norm)) : amp       # Positive feedback current

                    tau_soma_ahp = (Csoma_ahp * Ut) / (kappa * Isoma_ahp_tau) : second              # Time constant of adaptation
                    tau_soma = (Csoma_mem * Ut) / (kappa * Isoma_dpi_tau_clip) : second             # Membrane time constant
                    kappa = (kn + kp) / 2 : 1
                    Vsoma_mem = Ut / kappa * log(Isoma_mem / I0*(Isoma_mem>I0) + 1*(Isoma_mem<=I0)) : volt           # Membrane voltage

                    # Substrate constants
                    kn      : 1     (shared, constant)                          # Subthreshold slope factor for nFETs
                    kp      : 1     (shared, constant)                          # Subthreshold slope factor for pFETs
                    Ut      : volt  (shared, constant)                          # Thermal voltage
                    I0      : amp   (shared, constant)                          # Dark current
                    alpha   : 1     (shared, constant)                          # Scaling factor equal to Ig/Itau 

                    # Soma constants
                    Csoma_mem       : farad (shared, constant)                  # Membrane capacitance
                    Isoma_dpi_tau   : amp   (shared, constant)                  # Leakage current
                    Isoma_dpi_g = alpha * Isoma_dpi_tau : amp                   # DPI threshold (low pass filter) (Igmem)
                    Isoma_th        : amp   (shared, constant)                  # Spiking threshold
                    Isoma_reset     : amp   (shared, constant)                  # Reset current
                    Isoma_const     : amp   (constant)                          # Additional input current similar to constant current injection
                    soma_refP       : second (shared, constant)                 # Refractory period

                    # Adaptation constants
                    Csoma_ahp       : farad (shared, constant)                  # Spike-frequency adaptation capacitance
                    Isoma_ahp_tau   : amp   (shared, constant)                  # Leakage current for spike-frequency adaptation
                    Isoma_ahp_th    : amp   (shared, constant)                  # Threshold for spike-frequency adaptation
                    Isoma_ahp_w     : amp   (constant)                          # Calcium current

                    # Positive feedback constants
                    Isoma_pfb_gain  : amp   (shared, constant)                  # Positive feedback gain
                    Isoma_pfb_th    : amp   (shared, constant)                  # Positive feedback threshold (since it is a DPI circuit)
                    Isoma_pfb_norm  : amp   (shared, constant)                  # Positive feedback normalization current

                    # Synaptic dynamics #########################################

                    # NMDA #######################################################
                    dInmda/dt = (-Inmda - Inmda_g_clip +\
                    2*I0*(Inmda<=I0))/(tau_nmda*((Inmda_g_clip/Inmda)+1)) : amp

                    Inmda_g_clip = I0*(Inmda<=I0) + Inmda_g*(Inmda>I0) : amp
                    Inmda_tau_clip = I0*(Inmda<=I0) + Inmda_tau*(Inmda>I0) : amp
                    Inmda_dp = Inmda / (1 + Inmda_thr / Isoma_mem_clip) : amp # Output current of the DPI Slow synapse 

                    Vnmda : volt  (constant)                                    # NMDA threshold voltage
                    Inmda_tau : amp (constant)                                  # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Inmda_g = alpha * Inmda_tau : amp                           # Current flowing through ?? sets the DPI's threshold
                    Inmda_thr = I0 * exp(kappa * Vnmda / Ut) : amp
                    Inmda_w0 : amp (constant)                                   # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_nmda = Cnmda * Ut /(kappa * Inmda_tau_clip) : second    # Synaptic time-constant
                    Cnmda : farad (constant)                                    # Synapse's capacitance

                    # AMPA #######################################################
                    dIampa/dt = (-Iampa - Iampa_g_clip + 2*I0*(Iampa<=I0))/(tau_ampa*((Iampa_g_clip/Iampa)+1)) : amp

                    Iampa_g_clip = I0*(Iampa<=I0) + Iampa_g*(Iampa>I0) : amp
                    Iampa_tau_clip = I0*(Iampa<=I0) + Iampa_tau*(Iampa>I0) : amp

                    Iampa_tau : amp (constant)                                  # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Iampa_g = alpha * Iampa_tau : amp                           # Current flowing through ?? sets the DPI's threshold
                    Iampa_w0 : amp (constant)                                   # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_ampa = Campa * Ut / (kappa * Iampa_tau_clip) : second   # Synaptic time-constant
                    Campa          : farad (constant)                           # Synapse's capacitance

                    # GABA B - inh #######################################################
                    # the ihn synapse does not actually decrease Imem, it just
                    # decreases the input current from other synapses
                    dIgaba_b/dt = (-Igaba_b - Igaba_b_g_clip +\
                        2*I0*(Igaba_b<=I0))/(tau_gaba_b * ((Igaba_b_g_clip/Igaba_b)+1)) : amp

                    Igaba_b_g_clip = I0*(Igaba_b<=I0) + Igaba_b_g*(Igaba_b>I0) : amp
                    Igaba_b_tau_clip  = I0*(Igaba_b<=I0) + Igaba_b_tau*(Igaba_b>I0) : amp

                    Igaba_b_tau      : amp (constant)                           # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_b_g = alpha * Igaba_b_tau : amp                       # Current flowing through ?? sets the DPI's threshold
                    Igaba_b_w0 : amp (constant)                                 # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    tau_gaba_b = Cgaba_b * Ut / (kappa * Igaba_b_tau_clip) : second    # Synaptic time-constant
                    Cgaba_b          : farad (constant)                         # Synapse's capacitance

                    # GABA A - shunt #####################################################
                    dIgaba_a/dt =(-Igaba_a - Igaba_a_g_clip +\
                        2*I0*(Igaba_a<=I0))/(tau_gaba_a * ((Igaba_a_g_clip/Igaba_a)+1)) : amp

                    Igaba_a_g_clip = I0*(Igaba_a<=I0) + Igaba_a_g*(Igaba_a>I0) : amp  # DPI's gain factor
                    Igaba_a_tau_clip = I0*(Igaba_a<=I0) + Igaba_a_tau*(Igaba_a>I0) : amp

                    Igaba_a_tau     : amp (constant)                            # Leakage current, i.e. how much current is constantly leaked away (time-constant)
                    Igaba_a_g = alpha * Igaba_a_tau : amp                       # Current flowing through ?? sets the DPI's threshold
                    Igaba_a_w0 : amp (constant)                                 # Synaptic weight, to convert unitless weight to current
                    tau_gaba_a = Cgaba_a * Ut / (kappa * Igaba_a_tau_clip) : second     # Synaptic time-constant
                    Cgaba_a         : farad (constant)                          # Synapse's capacitance
                    ''',
           'threshold': '''Isoma_mem > Isoma_th''',
           'reset': '''
                    Isoma_mem = Isoma_reset
                    Isoma_ahp += Isoma_ahp_max
                    ''',
           'refractory': 'soma_refP',
           'method': 'euler'}

#synapses

def dynapse_nmda_syn_eq(): # SLOW_EXC
    """This function returns the slow excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Inmda_post += Inmda_w0*weight*Inmda_g_post/(Inmda_tau_post*((Inmda_g_post/Inmda_post)+1))
                    """,
           'on_post': """ """,
           'method': 'euler'}

def dynapse_ampa_syn_eq(): # FAST_EXC
    """This function returns the fast excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Iampa_post += Iampa_w0_post*weight*Iampa_g_post/(Iampa_tau_post*((Iampa_g_post/Iampa_post)+1))
                    """,
           'on_post': """ """,
           'method': 'euler'}

def dynapse_gaba_b_syn_eq(): # SLOW_INH
    """This function returns the subtractive inhibitory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                     Igaba_b_post += Igaba_b_w0_post*weight*Igaba_b_g_post/(Igaba_b_tau_post*((Igaba_b_g_post/Igaba_b_post)+1))
                     """,
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_a_syn_eq(): # FAST_INH
    """This function returns the shunting synapse (a mixture of subtractive and divisive) equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    Igaba_a_post += Igaba_a_w0_post*weight*Igaba_a_g_post/(Igaba_a_tau_post*((Igaba_a_g_post/Igaba_a_post)+1))
                    """,            # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}


