import brian2
import numpy as np

def get_attrs(all_objects, attribute, args=None): 
    """
	Returns the requested attribute of all the objects as a list
	Useful for getting a specific monitored parameter for all monitored objects
	Does not work with SpikeMonitors
	"""
    if(args==None):
        # If the requested attribute is a variable
        return np.array([getattr(obj, attribute) for obj in all_objects])
    else:
        # If the requested attribute is a method
        return np.array([getattr(obj, attribute)(*args) for obj in all_objects])
