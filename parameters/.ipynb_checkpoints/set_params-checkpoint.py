# -*- coding: utf-8 -*-
# @Author: alpren
# @Date:   2017-12-09 18:06:31
# @Last Modified by:   mmilde
# @Last Modified time: 2018-02-10 12:03:00

"""
This function allows to easily set parameters to neuron/synapse group
by providing a dictionary with keys (parameter names) and values
"""

import warnings


def set_params(group, params, raise_error = True):
    '''function that allows you to set parameters directly
    from the dictionary that is provided, instead of setting each parameter manually

    Args:
        group (brian2.groups.group): Neuron or synapse group
        params (dict): Dictionary with parameters to be set
    '''
    for par in params:
        if hasattr(group, par):
            setattr(group, par, params[par])
        else:
            warnings.warn("Group " + str(group.name) +
                          " has no state variable " + str(par))
            if raise_error:
                raise AttributeError("Group " + str(group.name) +
                          " has no state variable " + str(par) +
                          ', but you tried to set it with set_params '+
                          'if you want to ignore this error, pass raise_error = False')
