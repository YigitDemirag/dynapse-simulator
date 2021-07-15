# -*- coding: utf-8 -*-
# Copyright 2020 Yigit Demirag, Elisa Donati, Giacomo Indiveri
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import brian2
from brian2 import *

from parameters.dynapse_param import dynapse_param
from parameters.set_params import set_params
from equations.dynapse_eq import *
import numpy as np 
import sys 

NEURONS_PER_CORE = 256
NUM_CORES = 4

class DynapSE:
    def __init__(self, Network, num_cores=NUM_CORES, neurons_per_core=NEURONS_PER_CORE):
        
        self.num_cores         = num_cores
        self.neurons_per_core  = neurons_per_core
        self.core_usage        = np.zeros(num_cores)
        
        self.Network           = Network
        self.cores             = [Core(neurons_per_core, Network, name='Core_'+str(i)) for i in range(num_cores)]
        self.num_connections   = 0 

    def get_neurons(self, num_neurons, core_name=''):
        '''
        Try to allocate neurons from single core, if not enough, across cores (TBC)
        '''
        assert num_neurons < self.get_num_free_neurons_chip(), \
            'Requested {} neurons but number of free neurons on DynapSE cores are {}, total {}.'\
            .format(num_neurons,self.neurons_per_core-self.core_usage, np.sum(self.neurons_per_core-self.core_usage))
            
        if core_name == '': 
            if np.sum(num_neurons<=(self.neurons_per_core-self.core_usage)): # First, check available single core
                available_core = np.where(num_neurons<=(self.neurons_per_core-self.core_usage))[0][0]
                neurons = self.cores[available_core].get_neurons(num_neurons)
                print('{} neurons are allocated from Core_{}.'.format(num_neurons, available_core))

            else: # Second, allocate neurons across cores
                sys.exit('Tried to allocate {} neurons but not enough free neurons per single core. Maximum left is {}.'.format(num_neurons, max(self.neurons_per_core-self.core_usage)))
        
        else:
            core_names = [self.cores[i].name for i in range(self.num_cores)]
            assert core_name in core_names, "Please provide a valid core name e.g. 'Core_1', 'Core_2', etc. "
            core_idx = core_names.index(core_name)
            assert (self.neurons_per_core-self.cores[core_idx].usage) >= num_neurons, "Not enough space in Core_{} for {} neurons.".format(core_idx, num_neurons)
            neurons = self.cores[core_idx].get_neurons(num_neurons)
            print('{} neurons are allocated from Core_{}.'.format(num_neurons, core_idx))

        self.update_core_usage()
        return neurons  
    
    def set_bias(self, biases=dict(), core_name=''):
        ''' Sets neuron or synapse biases for core
        '''
        core_names = [self.cores[i].name for i in range(self.num_cores)]
        assert core_name in core_names, "Please provide a valid core name e.g. 'Core_1', 'Core_2', etc. "
        core_idx = core_names.index(core_name)

        for key in biases.keys():
            assert key in dynapse_param, "Bias {} is not a valid name.".format(key)
            
        set_params(self.cores[core_idx].neurons, biases)
        print('New bias values are loaded to {}.'.format(core_names[core_idx]))
        
        
    def add_connection(self, NeuronGroup1, NeuronGroup2, synapse_type='NMDA'):
        ''' Add synaptic connections between different NeuronGroups
        '''
        assert synapse_type in ['NMDA', 'AMPA', 'GABA_A','GABA_B'],'Please provide a valid synapse type.'
        n11,n12 = NeuronGroup1.start, NeuronGroup1.stop
        n21,n22 = NeuronGroup2.start, NeuronGroup2.stop
        
        core_id = self.get_core_id_of_pop(NeuronGroup2)
        self.cores[core_id].blocks.update({synapse_type : [n21, n22]})
       
        if synapse_type == 'NMDA':
            synapse = Synapses(NeuronGroup1, NeuronGroup2,
                    **dynapse_nmda_syn_eq(), name='NMDA'+str(self.num_connections))
        elif synapse_type == 'AMPA':
            synapse = Synapses(NeuronGroup1, NeuronGroup2,
                    **dynapse_ampa_syn_eq(), name='AMPA'+str(self.num_connections))
        elif synapse_type == 'GABA_A':
            synapse = Synapses(NeuronGroup1, NeuronGroup2,
                    **dynapse_gaba_a_syn_eq(), name='GABA_A'+str(self.num_connections))
        else: 
            synapse = Synapses(NeuronGroup1, NeuronGroup2,
                    **dynapse_gaba_b_syn_eq(), name='GABA_B'+str(self.num_connections))
        
        self.num_connections += 1
            
        return synapse
    
    def connect(self, Synapse, j='i'):
        n11,n12 = Synapse.source.start,Synapse.source.stop
        n21,n22 = Synapse.target.start,Synapse.target.stop
        if isinstance(j, bool):
            if j:
                assert (n12 - n11) < 64 and (n22 - n11) < 64, 'Pre and Post neurons should have less then 64 connections.'
                Synapse.connect(True)
        if isinstance(j, str):
            if j=='i':
                Synapse.connect(j=j)
            
    
    def get_num_free_neurons_chip(self):
        ''' Returns the number of free neurons per single chip
        '''
        return self.neurons_per_core * self.num_cores - np.sum(self.core_usage)
    
    def update_core_usage(self):
        ''' Updates the number of neurons allocated in the core
        '''
        self.core_usage = np.array([self.cores[i].usage for i in range(self.num_cores)])
        
    def get_core_id_of_pop(self, NeuronGroup):
        ''' Given a Brian NeuronGroup, returns the core they are allocated from
        '''
        core_idx = [self.cores[i].name in NeuronGroup.name for i in range(self.num_cores)].index(True)
        return core_idx
    
class Core:
    def __init__(self, num_neurons, Network, name=''):
        
        self.neurons     = NeuronGroup(num_neurons, **dynapse_eq(), name=name)
        self.num_neurons = num_neurons
        self.name        = name
        self.usage       = 0
        self.blocks      = dict()
        self.Network     = Network
        
        # Init
        self.set_default_core_params()
        self.Network.add(self.neurons)
        
    def set_default_core_params(self):
        ''' Sets default bias parameters
        '''
        set_params(self.neurons, dynapse_param)
                
    def get_neurons(self, num_neurons):
        neuron_slice = slice(self.usage, self.usage+num_neurons)
        neurons = self.neurons[neuron_slice]
        self.usage += num_neurons
        return neurons
    