#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alpha Renner, Yigit Demirag

Modified from https://code.ini.uzh.ch/ncs/teili
"""
#neuron
def dynapse_eq():
    return{'model': '''

                    dImem/dt = (((Ith_clip / Itau_clip) * (Iin_clip  + Ia_clip -\
                    Ishunt_clip - Iahp_clip)) - Ith_clip - ((1 + ((Ishunt_clip +\
                    Iahp_clip - Ia_clip) / Itau_clip)) * Imem)   ) / (tau *\
                    ((Ith_clip/(Imem + Io)) + 1)) : amp (unless refractory)

                    dIahp/dt = (- Ithahp_clip - Iahp + 2*Io*(Iahp<=Io)) / (tauahp *\
                    (Ithahp_clip / Iahp + 1)) : amp # adaptation current

                    # The *_clip currents are needed to prevent current from going
                    # below Io, since negative currents are not possible on chips
                    Itau_clip = Itau*(Imem>Io) + Io*(Imem<=Io)  : amp
                    Ith_clip = Ith*(Imem>Io) + Io*(Imem<=Io)    : amp
                    Iin_clip = clip(I_syn_nmda+I_syn_ampa-I_syn_gaba_a+Iconst,Io, 1*amp) : amp
                    Iahp_clip = Iahp*(Imem>Io) + Io*(Imem<=Io)  : amp
                    Ia_clip = Ia*(Imem>Io) + 2*Io*(Imem<=Io)    : amp
                    Ithahp_clip = Ithahp*(Iahp>Io) + Io*(Iahp<=Io) : amp
                    Ishunt_clip = clip(I_syn_gaba_b, Io, Imem) : amp

                    Iahpmax = (Ica / Itauahp) * Ithahp_clip : amp                # Ratio of currents through diffpair and adaptation block
                    Ia = Iagain / (1 + exp(-(Imem - Iath) / Ianorm)) : amp       # Positive feedback current

                    tauahp = (Cahp * Ut) / (kappa * Itauahp) : second       # Time constant of adaptation
                    tau = (Cmem * Ut) / (kappa * Itau_clip) : second        # Membrane time constant
                    kappa = (kn + kp) / 2 : 1

                    # constants
                    kn      : 1     (shared, constant)                 # Subthreshold slope factor for nFETs
                    kp      : 1     (shared, constant)                 # Subthreshold slope factor for pFETs
                    Ut      : volt  (shared, constant)                 # Thermal voltage
                    Io      : amp   (shared, constant)                 # Dark current
                    Cmem    : farad (shared, constant)                 # Membrane capacitance
                    Ispkthr : amp   (shared, constant)                 # Spiking threshold
                    Ireset  : amp   (shared, constant)                 # Reset current
                    refP    : second (shared, constant)                # Refractory period
                    Ith     : amp   (shared, constant)                 # DPI threshold (low pass filter)
                    Itau    : amp   (shared, constant)                 # Leakage current
                    Iconst  : amp   (constant)                         # Additional input current similar to constant current injection
                    Ithahp  : amp   (shared, constant)                 # Threshold for spike-frequency adaptation
                    Itauahp : amp   (shared, constant)                 # Leakage current for spike-frequency adaptation
                    Cahp    : farad (shared, constant)                 # Spike-frequency adaptation capacitance
                    tauca   : second (shared, constant)                # Calcium time-constant
                    Iagain  : amp   (shared, constant)                 # Positive feedback gain
                    Iath    : amp   (shared, constant)                 # Positive feedback threshold (since it is a DPI circuit)
                    Ianorm  : amp   (shared, constant)                 # Positive feedback normailzation current
                    Ishunt  : amp   (constant)                         # Shunting inhibitory synapic current
                    Ica     : amp   (constant)                         # Calcium current

                    #Synaptic dynamics #########################################

                    #exc #######################################################
                    dI_syn_nmda/dt = (-I_syn_nmda - I_g_clip_syn_exc +\
                    2*Io*(I_syn_nmda<=Io))/(tau_syn_exc*((I_g_clip_syn_exc/I_syn_nmda)+1)) : amp

                    I_g_clip_syn_exc = Io*(I_syn_nmda<=Io) + I_g_syn_nmda*(I_syn_nmda>Io) : amp
                    I_tau_clip_syn_exc = Io*(I_syn_nmda<=Io) + I_tau_syn_nmda*(I_syn_nmda>Io) : amp

                    I_wo_syn_nmda : amp (constant)                       # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    I_tau_syn_nmda : amp (constant)                      # Leakage current, i.e. how much current is constantly leaked away (time-cosntant)
                    I_g_syn_nmda       : amp (constant)                 # Current flowing through ?? sets the DPI's threshold
                    tau_syn_exc = C_syn_nmda * Ut /(kappa * I_tau_clip_syn_exc) : second    # Synaptic time-constant
                    C_syn_nmda          : farad (constant)               # Synapse's capacitance
                    
                    
                    #exc2 #######################################################
                    dI_syn_ampa/dt = (-I_syn_ampa - I_g_clip_syn_exc2 +\
                    2*Io*(I_syn_ampa<=Io))/(tau_syn_exc2*((I_g_clip_syn_exc2/I_syn_ampa)+1)) : amp

                    I_g_clip_syn_exc2 = Io*(I_syn_ampa<=Io) + I_g_syn_ampa*(I_syn_ampa>Io) : amp
                    I_tau_clip_syn_exc2 = Io*(I_syn_ampa<=Io) + I_tau_syn_ampa*(I_syn_ampa>Io) : amp

                    I_wo_syn_ampa : amp (constant)                       # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    I_tau_syn_ampa : amp (constant)                      # Leakage current, i.e. how much current is constantly leaked away (time-cosntant)
                    I_g_syn_ampa       : amp (constant)                 # Current flowing through ?? sets the DPI's threshold
                    tau_syn_exc2 = C_syn_ampa * Ut /(kappa * I_tau_clip_syn_exc2) : second    # Synaptic time-constant
                    C_syn_ampa          : farad (constant)               # Synapse's capacitance

                    #inh #######################################################
                    # the ihn synapse does not actually decrease Imem, it just
                    # decreases the input current from other synapses
                    dI_syn_gaba_a/dt = (-I_syn_gaba_a - I_g_clip_syn_inh +\
                    2*Io*(I_syn_gaba_a<=Io))/(tau_syn_inh *((I_g_clip_syn_inh/I_syn_gaba_a)+1)) : amp

                    I_g_clip_syn_inh = Io*(I_syn_gaba_a<=Io) + I_g_syn_gaba_a*(I_syn_gaba_a>Io) : amp
                    I_tau_clip_syn_inh  = Io*(I_syn_gaba_a<=Io) + I_tau_syn_gaba_a*(I_syn_gaba_a>Io) : amp

                    I_wo_syn_gaba_a : amp (constant)                       # Base synaptic weight, to convert unitless weight (set in synapse) to current
                    I_tau_syn_gaba_a      : amp (constant)                 # Leakage current, i.e. how much current is constantly leaked away (time-cosntant)
                    I_g_syn_gaba_a       : amp (constant)                 # Current flowing through ?? sets the DPI's threshold
                    tau_syn_inh  = C_syn_gaba_a * Ut /(kappa * I_tau_clip_syn_inh) : second    # Synaptic time-constant
                    C_syn_gaba_a          : farad (constant)               # Synapse's capacitance

                    #shunt #####################################################
                    dI_syn_gaba_b/dt =(-I_syn_gaba_b - I_g_clip_syn_shunt +\
                    2*Io*(I_syn_gaba_b<=Io))/(tau_syn_shunt*((I_g_clip_syn_shunt/I_syn_gaba_b)+1)) : amp

                    I_g_clip_syn_shunt = Io*(I_syn_gaba_b<=Io) +\
                    I_g_syn_gaba_b*(I_syn_gaba_b>Io) : amp  # DPI's gain factor

                    I_tau_clip_syn_shunt = Io*(I_syn_gaba_b<=Io) + I_tau_syn_gaba_b*(I_syn_gaba_b>Io) : amp
                    I_wo_syn_gaba_b : amp (constant)             # Synaptic weight, to convert unitless weight to current
                    tau_syn_shunt = C_syn_gaba_b * Ut /(kappa * I_tau_clip_syn_shunt) : second     # Synaptic time-constant

                    I_tau_syn_gaba_b       : amp (constant)    # Leakage current, i.e. how much current is constantly leaked away (time-cosntant)
                    I_g_syn_gaba_b        : amp (constant)    # Current flowing through ?? sets the DPI's threshold
                    C_syn_gaba_b         : farad (constant)    # Synapse's capacitance
                    ''',
           'threshold': '''Imem > Ispkthr''',
           'reset': '''
                    Imem = Ireset
                    Iahp += Iahpmax
                    ''',
           'refractory': 'refP',
           'method': 'euler'}

#synapses

def dynapse_nmda_syn_eq(): # SLOW_EXC
    """This function returns the slow excitatory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                    I_syn_nmda_post += I_wo_syn_nmda_post*weight*I_g_syn_nmda_post/(I_tau_syn_nmda_post*((I_g_syn_nmda_post/I_syn_nmda_post)+1))
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
                    I_syn_ampa_post += I_wo_syn_ampa_post*weight*I_g_syn_ampa_post/(I_tau_syn_ampa_post*((I_g_syn_ampa_post/I_syn_ampa_post)+1))
                    """,
           'on_post': """ """,
           'method': 'euler'}

def dynapse_gaba_a_syn_eq(): # SLOW_INH
    """This function returns the inhibitory synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                     I_syn_gaba_a_post += I_wo_syn_gaba_a_post*-weight*I_g_syn_gaba_a_post/(I_tau_syn_gaba_a_post*((I_g_syn_gaba_a_post/I_syn_gaba_a_post)+1))
                     """,
           'on_post': """ """,
           'method': 'euler'}


def dynapse_gaba_b_syn_eq(): # FAST_INH
    """This function returns the shunting synapse equation dictionary.
    """
    return{'model': """
                    weight : 1 # Can only be integer on the chip
                    """,
           'on_pre': """
                     I_syn_gaba_b_post += I_wo_syn_gaba_b_post*-weight*I_g_syn_gaba_b_post/(I_tau_syn_gaba_b_post*((I_g_syn_gaba_b_post/I_syn_gaba_b_post)+1))
                     """,            # On pre-synaptic spike adds current to state variable of DPI synapse
           'on_post': """ """,
           'method': 'euler'}


