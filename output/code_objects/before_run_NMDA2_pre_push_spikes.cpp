#include "code_objects/before_run_NMDA2_pre_push_spikes.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>

void _before_run_NMDA2_pre_push_spikes()
{
    using namespace brian;
    ///// CONSTANTS ///////////
    double* const _array_NMDA2_delay = _dynamic_array_NMDA2_delay.empty()? 0 : &_dynamic_array_NMDA2_delay[0];
const size_t _numdelay = _dynamic_array_NMDA2_delay.size();
const size_t _num_source_dt = 1;
const size_t _num_spikespace = 201;
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_NMDA2_delay = _array_NMDA2_delay;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    int32_t* __restrict  _ptr_array_poissongroup_2__spikespace = _array_poissongroup_2__spikespace;

    std::vector<double> &real_delays = _dynamic_array_NMDA2_delay;
    double* real_delays_data = real_delays.empty() ? 0 : &(real_delays[0]);
    int32_t* sources = NMDA2_pre.sources.empty() ? 0 : &(NMDA2_pre.sources[0]);
    const size_t n_delays = real_delays.size();
    const size_t n_synapses = NMDA2_pre.sources.size();
    NMDA2_pre.prepare(200,
                           200,
                           real_delays_data, n_delays, sources,
                           n_synapses,
                           _ptr_array_defaultclock_dt[0]);

}
