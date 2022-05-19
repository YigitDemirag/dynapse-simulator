import os
import sys 
sys.path.append('..')

from brian2 import *

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

	directory = f"mp_run_{pid}_{date_time}"
		
	set_device('cpp_standalone', directory=directory)

	BrianLogger.suppress_name('base')
	defaultclock.dt = 20 * us
	
	net = Network()
	chip = DynapSE(net)

	neuron = chip.get_neurons(1, 'Core_1')
	
	net.add(neuron)

	net.run(params['_exp_duration_'])
	
	result = {
		#Monitors results
	}

	device.reinit()
	
	return result
