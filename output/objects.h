
#ifndef _BRIAN_OBJECTS_H
#define _BRIAN_OBJECTS_H

#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>


namespace brian {

// In OpenMP we need one state per thread
extern std::vector< rk_state* > _mersenne_twister_states;

//////////////// clocks ///////////////////
extern Clock defaultclock;

//////////////// networks /////////////////
extern Network network_1;

//////////////// dynamic arrays ///////////
extern std::vector<int32_t> _dynamic_array_AMPA0__synaptic_post;
extern std::vector<int32_t> _dynamic_array_AMPA0__synaptic_pre;
extern std::vector<double> _dynamic_array_AMPA0_delay;
extern std::vector<double> _dynamic_array_AMPA0_delay_1;
extern std::vector<int32_t> _dynamic_array_AMPA0_N_incoming;
extern std::vector<int32_t> _dynamic_array_AMPA0_N_outgoing;
extern std::vector<double> _dynamic_array_AMPA0_weight;
extern std::vector<int32_t> _dynamic_array_InpSpikeGenerator__timebins;
extern std::vector<int32_t> _dynamic_array_InpSpikeGenerator_neuron_index;
extern std::vector<int32_t> _dynamic_array_InpSpikeGenerator_spike_number;
extern std::vector<double> _dynamic_array_InpSpikeGenerator_spike_time;
extern std::vector<int32_t> _dynamic_array_mon_neuron_input_i;
extern std::vector<double> _dynamic_array_mon_neuron_input_t;
extern std::vector<int32_t> _dynamic_array_mon_neuron_output_i;
extern std::vector<double> _dynamic_array_mon_neuron_output_t;
extern std::vector<double> _dynamic_array_statemonitor_1_t;
extern std::vector<double> _dynamic_array_statemonitor_2_t;

//////////////// arrays ///////////////////
extern int32_t *_array_AMPA0_N;
extern const int _num__array_AMPA0_N;
extern int32_t *_array_Core_0__spikespace;
extern const int _num__array_Core_0__spikespace;
extern double *_array_Core_0_alpha;
extern const int _num__array_Core_0_alpha;
extern double *_array_Core_0_Campa;
extern const int _num__array_Core_0_Campa;
extern double *_array_Core_0_Cgaba_a;
extern const int _num__array_Core_0_Cgaba_a;
extern double *_array_Core_0_Cgaba_b;
extern const int _num__array_Core_0_Cgaba_b;
extern double *_array_Core_0_Cnmda;
extern const int _num__array_Core_0_Cnmda;
extern double *_array_Core_0_Csoma_ahp;
extern const int _num__array_Core_0_Csoma_ahp;
extern double *_array_Core_0_Csoma_mem;
extern const int _num__array_Core_0_Csoma_mem;
extern int32_t *_array_Core_0_i;
extern const int _num__array_Core_0_i;
extern double *_array_Core_0_I0;
extern const int _num__array_Core_0_I0;
extern double *_array_Core_0_Iampa;
extern const int _num__array_Core_0_Iampa;
extern double *_array_Core_0_Iampa_tau;
extern const int _num__array_Core_0_Iampa_tau;
extern double *_array_Core_0_Iampa_w0;
extern const int _num__array_Core_0_Iampa_w0;
extern double *_array_Core_0_Igaba_a;
extern const int _num__array_Core_0_Igaba_a;
extern double *_array_Core_0_Igaba_a_tau;
extern const int _num__array_Core_0_Igaba_a_tau;
extern double *_array_Core_0_Igaba_a_w0;
extern const int _num__array_Core_0_Igaba_a_w0;
extern double *_array_Core_0_Igaba_b;
extern const int _num__array_Core_0_Igaba_b;
extern double *_array_Core_0_Igaba_b_tau;
extern const int _num__array_Core_0_Igaba_b_tau;
extern double *_array_Core_0_Igaba_b_w0;
extern const int _num__array_Core_0_Igaba_b_w0;
extern double *_array_Core_0_Inmda;
extern const int _num__array_Core_0_Inmda;
extern double *_array_Core_0_Inmda_tau;
extern const int _num__array_Core_0_Inmda_tau;
extern double *_array_Core_0_Inmda_w0;
extern const int _num__array_Core_0_Inmda_w0;
extern double *_array_Core_0_Isoma_ahp;
extern const int _num__array_Core_0_Isoma_ahp;
extern double *_array_Core_0_Isoma_ahp_tau;
extern const int _num__array_Core_0_Isoma_ahp_tau;
extern double *_array_Core_0_Isoma_ahp_th;
extern const int _num__array_Core_0_Isoma_ahp_th;
extern double *_array_Core_0_Isoma_ahp_w;
extern const int _num__array_Core_0_Isoma_ahp_w;
extern double *_array_Core_0_Isoma_const;
extern const int _num__array_Core_0_Isoma_const;
extern double *_array_Core_0_Isoma_dpi_tau;
extern const int _num__array_Core_0_Isoma_dpi_tau;
extern double *_array_Core_0_Isoma_mem;
extern const int _num__array_Core_0_Isoma_mem;
extern double *_array_Core_0_Isoma_pfb_gain;
extern const int _num__array_Core_0_Isoma_pfb_gain;
extern double *_array_Core_0_Isoma_pfb_norm;
extern const int _num__array_Core_0_Isoma_pfb_norm;
extern double *_array_Core_0_Isoma_pfb_th;
extern const int _num__array_Core_0_Isoma_pfb_th;
extern double *_array_Core_0_Isoma_reset;
extern const int _num__array_Core_0_Isoma_reset;
extern double *_array_Core_0_Isoma_th;
extern const int _num__array_Core_0_Isoma_th;
extern double *_array_Core_0_kn;
extern const int _num__array_Core_0_kn;
extern double *_array_Core_0_kp;
extern const int _num__array_Core_0_kp;
extern double *_array_Core_0_lastspike;
extern const int _num__array_Core_0_lastspike;
extern char *_array_Core_0_not_refractory;
extern const int _num__array_Core_0_not_refractory;
extern double *_array_Core_0_soma_refP;
extern const int _num__array_Core_0_soma_refP;
extern double *_array_Core_0_Ut;
extern const int _num__array_Core_0_Ut;
extern double *_array_Core_0_Vnmda;
extern const int _num__array_Core_0_Vnmda;
extern int32_t *_array_Core_1__spikespace;
extern const int _num__array_Core_1__spikespace;
extern double *_array_Core_1_alpha;
extern const int _num__array_Core_1_alpha;
extern double *_array_Core_1_Campa;
extern const int _num__array_Core_1_Campa;
extern double *_array_Core_1_Cgaba_a;
extern const int _num__array_Core_1_Cgaba_a;
extern double *_array_Core_1_Cgaba_b;
extern const int _num__array_Core_1_Cgaba_b;
extern double *_array_Core_1_Cnmda;
extern const int _num__array_Core_1_Cnmda;
extern double *_array_Core_1_Csoma_ahp;
extern const int _num__array_Core_1_Csoma_ahp;
extern double *_array_Core_1_Csoma_mem;
extern const int _num__array_Core_1_Csoma_mem;
extern int32_t *_array_Core_1_i;
extern const int _num__array_Core_1_i;
extern double *_array_Core_1_I0;
extern const int _num__array_Core_1_I0;
extern double *_array_Core_1_Iampa;
extern const int _num__array_Core_1_Iampa;
extern double *_array_Core_1_Iampa_tau;
extern const int _num__array_Core_1_Iampa_tau;
extern double *_array_Core_1_Iampa_w0;
extern const int _num__array_Core_1_Iampa_w0;
extern double *_array_Core_1_Igaba_a;
extern const int _num__array_Core_1_Igaba_a;
extern double *_array_Core_1_Igaba_a_tau;
extern const int _num__array_Core_1_Igaba_a_tau;
extern double *_array_Core_1_Igaba_a_w0;
extern const int _num__array_Core_1_Igaba_a_w0;
extern double *_array_Core_1_Igaba_b;
extern const int _num__array_Core_1_Igaba_b;
extern double *_array_Core_1_Igaba_b_tau;
extern const int _num__array_Core_1_Igaba_b_tau;
extern double *_array_Core_1_Igaba_b_w0;
extern const int _num__array_Core_1_Igaba_b_w0;
extern double *_array_Core_1_Inmda;
extern const int _num__array_Core_1_Inmda;
extern double *_array_Core_1_Inmda_tau;
extern const int _num__array_Core_1_Inmda_tau;
extern double *_array_Core_1_Inmda_w0;
extern const int _num__array_Core_1_Inmda_w0;
extern double *_array_Core_1_Isoma_ahp;
extern const int _num__array_Core_1_Isoma_ahp;
extern double *_array_Core_1_Isoma_ahp_tau;
extern const int _num__array_Core_1_Isoma_ahp_tau;
extern double *_array_Core_1_Isoma_ahp_th;
extern const int _num__array_Core_1_Isoma_ahp_th;
extern double *_array_Core_1_Isoma_ahp_w;
extern const int _num__array_Core_1_Isoma_ahp_w;
extern double *_array_Core_1_Isoma_const;
extern const int _num__array_Core_1_Isoma_const;
extern double *_array_Core_1_Isoma_dpi_tau;
extern const int _num__array_Core_1_Isoma_dpi_tau;
extern double *_array_Core_1_Isoma_mem;
extern const int _num__array_Core_1_Isoma_mem;
extern double *_array_Core_1_Isoma_pfb_gain;
extern const int _num__array_Core_1_Isoma_pfb_gain;
extern double *_array_Core_1_Isoma_pfb_norm;
extern const int _num__array_Core_1_Isoma_pfb_norm;
extern double *_array_Core_1_Isoma_pfb_th;
extern const int _num__array_Core_1_Isoma_pfb_th;
extern double *_array_Core_1_Isoma_reset;
extern const int _num__array_Core_1_Isoma_reset;
extern double *_array_Core_1_Isoma_th;
extern const int _num__array_Core_1_Isoma_th;
extern double *_array_Core_1_kn;
extern const int _num__array_Core_1_kn;
extern double *_array_Core_1_kp;
extern const int _num__array_Core_1_kp;
extern double *_array_Core_1_lastspike;
extern const int _num__array_Core_1_lastspike;
extern char *_array_Core_1_not_refractory;
extern const int _num__array_Core_1_not_refractory;
extern double *_array_Core_1_soma_refP;
extern const int _num__array_Core_1_soma_refP;
extern int32_t *_array_Core_1_subgroup_1__sub_idx;
extern const int _num__array_Core_1_subgroup_1__sub_idx;
extern double *_array_Core_1_Ut;
extern const int _num__array_Core_1_Ut;
extern double *_array_Core_1_Vnmda;
extern const int _num__array_Core_1_Vnmda;
extern int32_t *_array_Core_2__spikespace;
extern const int _num__array_Core_2__spikespace;
extern double *_array_Core_2_alpha;
extern const int _num__array_Core_2_alpha;
extern double *_array_Core_2_Campa;
extern const int _num__array_Core_2_Campa;
extern double *_array_Core_2_Cgaba_a;
extern const int _num__array_Core_2_Cgaba_a;
extern double *_array_Core_2_Cgaba_b;
extern const int _num__array_Core_2_Cgaba_b;
extern double *_array_Core_2_Cnmda;
extern const int _num__array_Core_2_Cnmda;
extern double *_array_Core_2_Csoma_ahp;
extern const int _num__array_Core_2_Csoma_ahp;
extern double *_array_Core_2_Csoma_mem;
extern const int _num__array_Core_2_Csoma_mem;
extern int32_t *_array_Core_2_i;
extern const int _num__array_Core_2_i;
extern double *_array_Core_2_I0;
extern const int _num__array_Core_2_I0;
extern double *_array_Core_2_Iampa;
extern const int _num__array_Core_2_Iampa;
extern double *_array_Core_2_Iampa_tau;
extern const int _num__array_Core_2_Iampa_tau;
extern double *_array_Core_2_Iampa_w0;
extern const int _num__array_Core_2_Iampa_w0;
extern double *_array_Core_2_Igaba_a;
extern const int _num__array_Core_2_Igaba_a;
extern double *_array_Core_2_Igaba_a_tau;
extern const int _num__array_Core_2_Igaba_a_tau;
extern double *_array_Core_2_Igaba_a_w0;
extern const int _num__array_Core_2_Igaba_a_w0;
extern double *_array_Core_2_Igaba_b;
extern const int _num__array_Core_2_Igaba_b;
extern double *_array_Core_2_Igaba_b_tau;
extern const int _num__array_Core_2_Igaba_b_tau;
extern double *_array_Core_2_Igaba_b_w0;
extern const int _num__array_Core_2_Igaba_b_w0;
extern double *_array_Core_2_Inmda;
extern const int _num__array_Core_2_Inmda;
extern double *_array_Core_2_Inmda_tau;
extern const int _num__array_Core_2_Inmda_tau;
extern double *_array_Core_2_Inmda_w0;
extern const int _num__array_Core_2_Inmda_w0;
extern double *_array_Core_2_Isoma_ahp;
extern const int _num__array_Core_2_Isoma_ahp;
extern double *_array_Core_2_Isoma_ahp_tau;
extern const int _num__array_Core_2_Isoma_ahp_tau;
extern double *_array_Core_2_Isoma_ahp_th;
extern const int _num__array_Core_2_Isoma_ahp_th;
extern double *_array_Core_2_Isoma_ahp_w;
extern const int _num__array_Core_2_Isoma_ahp_w;
extern double *_array_Core_2_Isoma_const;
extern const int _num__array_Core_2_Isoma_const;
extern double *_array_Core_2_Isoma_dpi_tau;
extern const int _num__array_Core_2_Isoma_dpi_tau;
extern double *_array_Core_2_Isoma_mem;
extern const int _num__array_Core_2_Isoma_mem;
extern double *_array_Core_2_Isoma_pfb_gain;
extern const int _num__array_Core_2_Isoma_pfb_gain;
extern double *_array_Core_2_Isoma_pfb_norm;
extern const int _num__array_Core_2_Isoma_pfb_norm;
extern double *_array_Core_2_Isoma_pfb_th;
extern const int _num__array_Core_2_Isoma_pfb_th;
extern double *_array_Core_2_Isoma_reset;
extern const int _num__array_Core_2_Isoma_reset;
extern double *_array_Core_2_Isoma_th;
extern const int _num__array_Core_2_Isoma_th;
extern double *_array_Core_2_kn;
extern const int _num__array_Core_2_kn;
extern double *_array_Core_2_kp;
extern const int _num__array_Core_2_kp;
extern double *_array_Core_2_lastspike;
extern const int _num__array_Core_2_lastspike;
extern char *_array_Core_2_not_refractory;
extern const int _num__array_Core_2_not_refractory;
extern double *_array_Core_2_soma_refP;
extern const int _num__array_Core_2_soma_refP;
extern double *_array_Core_2_Ut;
extern const int _num__array_Core_2_Ut;
extern double *_array_Core_2_Vnmda;
extern const int _num__array_Core_2_Vnmda;
extern int32_t *_array_Core_3__spikespace;
extern const int _num__array_Core_3__spikespace;
extern double *_array_Core_3_alpha;
extern const int _num__array_Core_3_alpha;
extern double *_array_Core_3_Campa;
extern const int _num__array_Core_3_Campa;
extern double *_array_Core_3_Cgaba_a;
extern const int _num__array_Core_3_Cgaba_a;
extern double *_array_Core_3_Cgaba_b;
extern const int _num__array_Core_3_Cgaba_b;
extern double *_array_Core_3_Cnmda;
extern const int _num__array_Core_3_Cnmda;
extern double *_array_Core_3_Csoma_ahp;
extern const int _num__array_Core_3_Csoma_ahp;
extern double *_array_Core_3_Csoma_mem;
extern const int _num__array_Core_3_Csoma_mem;
extern int32_t *_array_Core_3_i;
extern const int _num__array_Core_3_i;
extern double *_array_Core_3_I0;
extern const int _num__array_Core_3_I0;
extern double *_array_Core_3_Iampa;
extern const int _num__array_Core_3_Iampa;
extern double *_array_Core_3_Iampa_tau;
extern const int _num__array_Core_3_Iampa_tau;
extern double *_array_Core_3_Iampa_w0;
extern const int _num__array_Core_3_Iampa_w0;
extern double *_array_Core_3_Igaba_a;
extern const int _num__array_Core_3_Igaba_a;
extern double *_array_Core_3_Igaba_a_tau;
extern const int _num__array_Core_3_Igaba_a_tau;
extern double *_array_Core_3_Igaba_a_w0;
extern const int _num__array_Core_3_Igaba_a_w0;
extern double *_array_Core_3_Igaba_b;
extern const int _num__array_Core_3_Igaba_b;
extern double *_array_Core_3_Igaba_b_tau;
extern const int _num__array_Core_3_Igaba_b_tau;
extern double *_array_Core_3_Igaba_b_w0;
extern const int _num__array_Core_3_Igaba_b_w0;
extern double *_array_Core_3_Inmda;
extern const int _num__array_Core_3_Inmda;
extern double *_array_Core_3_Inmda_tau;
extern const int _num__array_Core_3_Inmda_tau;
extern double *_array_Core_3_Inmda_w0;
extern const int _num__array_Core_3_Inmda_w0;
extern double *_array_Core_3_Isoma_ahp;
extern const int _num__array_Core_3_Isoma_ahp;
extern double *_array_Core_3_Isoma_ahp_tau;
extern const int _num__array_Core_3_Isoma_ahp_tau;
extern double *_array_Core_3_Isoma_ahp_th;
extern const int _num__array_Core_3_Isoma_ahp_th;
extern double *_array_Core_3_Isoma_ahp_w;
extern const int _num__array_Core_3_Isoma_ahp_w;
extern double *_array_Core_3_Isoma_const;
extern const int _num__array_Core_3_Isoma_const;
extern double *_array_Core_3_Isoma_dpi_tau;
extern const int _num__array_Core_3_Isoma_dpi_tau;
extern double *_array_Core_3_Isoma_mem;
extern const int _num__array_Core_3_Isoma_mem;
extern double *_array_Core_3_Isoma_pfb_gain;
extern const int _num__array_Core_3_Isoma_pfb_gain;
extern double *_array_Core_3_Isoma_pfb_norm;
extern const int _num__array_Core_3_Isoma_pfb_norm;
extern double *_array_Core_3_Isoma_pfb_th;
extern const int _num__array_Core_3_Isoma_pfb_th;
extern double *_array_Core_3_Isoma_reset;
extern const int _num__array_Core_3_Isoma_reset;
extern double *_array_Core_3_Isoma_th;
extern const int _num__array_Core_3_Isoma_th;
extern double *_array_Core_3_kn;
extern const int _num__array_Core_3_kn;
extern double *_array_Core_3_kp;
extern const int _num__array_Core_3_kp;
extern double *_array_Core_3_lastspike;
extern const int _num__array_Core_3_lastspike;
extern char *_array_Core_3_not_refractory;
extern const int _num__array_Core_3_not_refractory;
extern double *_array_Core_3_soma_refP;
extern const int _num__array_Core_3_soma_refP;
extern double *_array_Core_3_Ut;
extern const int _num__array_Core_3_Ut;
extern double *_array_Core_3_Vnmda;
extern const int _num__array_Core_3_Vnmda;
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern int64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_InpSpikeGenerator__lastindex;
extern const int _num__array_InpSpikeGenerator__lastindex;
extern int32_t *_array_InpSpikeGenerator__period_bins;
extern const int _num__array_InpSpikeGenerator__period_bins;
extern int32_t *_array_InpSpikeGenerator__spikespace;
extern const int _num__array_InpSpikeGenerator__spikespace;
extern int32_t *_array_InpSpikeGenerator_i;
extern const int _num__array_InpSpikeGenerator_i;
extern double *_array_InpSpikeGenerator_period;
extern const int _num__array_InpSpikeGenerator_period;
extern int32_t *_array_mon_neuron_input__source_idx;
extern const int _num__array_mon_neuron_input__source_idx;
extern int32_t *_array_mon_neuron_input_count;
extern const int _num__array_mon_neuron_input_count;
extern int32_t *_array_mon_neuron_input_N;
extern const int _num__array_mon_neuron_input_N;
extern int32_t *_array_mon_neuron_output__source_idx;
extern const int _num__array_mon_neuron_output__source_idx;
extern int32_t *_array_mon_neuron_output_count;
extern const int _num__array_mon_neuron_output_count;
extern int32_t *_array_mon_neuron_output_N;
extern const int _num__array_mon_neuron_output_N;
extern int32_t *_array_statemonitor_1__indices;
extern const int _num__array_statemonitor_1__indices;
extern double *_array_statemonitor_1_Isoma_ahp;
extern const int _num__array_statemonitor_1_Isoma_ahp;
extern int32_t *_array_statemonitor_1_N;
extern const int _num__array_statemonitor_1_N;
extern int32_t *_array_statemonitor_2__indices;
extern const int _num__array_statemonitor_2__indices;
extern double *_array_statemonitor_2_Isoma_mem;
extern const int _num__array_statemonitor_2_Isoma_mem;
extern int32_t *_array_statemonitor_2_N;
extern const int _num__array_statemonitor_2_N;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_statemonitor_1_Isoma_ahp;
extern DynamicArray2D<double> _dynamic_array_statemonitor_2_Isoma_mem;

/////////////// static arrays /////////////
extern int32_t *_static_array__dynamic_array_InpSpikeGenerator__timebins;
extern const int _num__static_array__dynamic_array_InpSpikeGenerator__timebins;
extern double *_static_array__dynamic_array_InpSpikeGenerator_neuron_index;
extern const int _num__static_array__dynamic_array_InpSpikeGenerator_neuron_index;
extern int64_t *_static_array__dynamic_array_InpSpikeGenerator_spike_number;
extern const int _num__static_array__dynamic_array_InpSpikeGenerator_spike_number;
extern double *_static_array__dynamic_array_InpSpikeGenerator_spike_time;
extern const int _num__static_array__dynamic_array_InpSpikeGenerator_spike_time;

//////////////// synapses /////////////////
// AMPA0
extern SynapticPathway AMPA0_post;
extern SynapticPathway AMPA0_pre;

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif

