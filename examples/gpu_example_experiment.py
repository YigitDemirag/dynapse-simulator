"""
@author: Ioan Fodorut

"""

import os
import sys 
sys.path.append('..')


from brian2 import *
import brian2genn

from DynapSE import DynapSE
from equations.dynapse_eq import *
from parameters.dynapse_param import *
from parameters.set_params import set_params
from utils.utils import get_attrs

def run_experiment(params):
	
	now = datetime.datetime.now() # current date and time
	date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
	pid = os.getpid()
	print(f'Process {pid} launched, at ' + date_time)

	directory = f"gpu_run_{pid}_{date_time}"
		
	set_device('genn', directory=directory)

	BrianLogger.suppress_name('base')
	defaultclock.dt = 20 * us

	# Do NOT use the TimedArrays!
	# They are not supported by brian2genn
	spikes = np.zeros(params['_exp_duration_'])
	dt = int(1000/params['_input_rate_'])
	spikes[0::dt] = 1.0
	spike_timing = np.where(spikes==1)[0] * ms
	neuron_indices = np.zeros(len(spike_timing))
	input_spike_generator = SpikeGeneratorGroup(params['_nr_inputs_'], indices=neuron_indices, times=spike_timing, name='InpSpikeGenerator') 

	# Do NOT use the "chip" to "get_neurons"! 
	# The chip functionality is implemented using subgroups, which are not supported by brain2geen
	neuron = NeuronGroup(1, **dynapse_eq())
	set_params(neuron, dynapse_param)

	# Do NOT use the "chip" to "add connections" or "connect"! 
	# It won't work because you cannot use the neurons of the "chip" 
	synapses = Synapses(input_spike_generator, neuron, **dynapse_ampa_syn_eq())
	synapses.connect()
	synapses.weight = params['_ampa_weight_']

	monitor_neuron_state  = StateMonitor(neuron, ['Isoma_mem_clip', 'Iampa'], record=True, name="monitor_neuron_state")
	monitor_neuron_output = SpikeMonitor(neuron, name="monitor_neuron_output")

	run(1 * second)
	
	result = {
		'neuron_membrane_current': [monitor_neuron_state.t/ms, get_attrs(monitor_neuron_state, 'Isoma_mem_clip')],
		'neuron_spike_output': array([monitor_neuron_output.t/ms, monitor_neuron_output.i]),
		'synaptic_current': [monitor_neuron_state.t/ms, get_attrs(monitor_neuron_state, 'Iampa')],
	}

	device.reinit()
	
	return result
