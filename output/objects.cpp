
#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>
#include<iostream>
#include<fstream>

namespace brian {

std::vector< rk_state* > _mersenne_twister_states;

//////////////// networks /////////////////
Network network;

//////////////// arrays ///////////////////
int32_t * _array_Core_0__spikespace;
const int _num__array_Core_0__spikespace = 257;
double * _array_Core_0_alpha;
const int _num__array_Core_0_alpha = 1;
double * _array_Core_0_Campa;
const int _num__array_Core_0_Campa = 256;
double * _array_Core_0_Cgaba_a;
const int _num__array_Core_0_Cgaba_a = 256;
double * _array_Core_0_Cgaba_b;
const int _num__array_Core_0_Cgaba_b = 256;
double * _array_Core_0_Cnmda;
const int _num__array_Core_0_Cnmda = 256;
double * _array_Core_0_Csoma_ahp;
const int _num__array_Core_0_Csoma_ahp = 1;
double * _array_Core_0_Csoma_mem;
const int _num__array_Core_0_Csoma_mem = 1;
int32_t * _array_Core_0_i;
const int _num__array_Core_0_i = 256;
double * _array_Core_0_I0;
const int _num__array_Core_0_I0 = 1;
double * _array_Core_0_Iampa;
const int _num__array_Core_0_Iampa = 256;
double * _array_Core_0_Iampa_tau;
const int _num__array_Core_0_Iampa_tau = 256;
double * _array_Core_0_Iampa_w0;
const int _num__array_Core_0_Iampa_w0 = 256;
double * _array_Core_0_Igaba_a;
const int _num__array_Core_0_Igaba_a = 256;
double * _array_Core_0_Igaba_a_tau;
const int _num__array_Core_0_Igaba_a_tau = 256;
double * _array_Core_0_Igaba_a_w0;
const int _num__array_Core_0_Igaba_a_w0 = 256;
double * _array_Core_0_Igaba_b;
const int _num__array_Core_0_Igaba_b = 256;
double * _array_Core_0_Igaba_b_tau;
const int _num__array_Core_0_Igaba_b_tau = 256;
double * _array_Core_0_Igaba_b_w0;
const int _num__array_Core_0_Igaba_b_w0 = 256;
double * _array_Core_0_Inmda;
const int _num__array_Core_0_Inmda = 256;
double * _array_Core_0_Inmda_tau;
const int _num__array_Core_0_Inmda_tau = 256;
double * _array_Core_0_Inmda_w0;
const int _num__array_Core_0_Inmda_w0 = 256;
double * _array_Core_0_Isoma_ahp;
const int _num__array_Core_0_Isoma_ahp = 256;
double * _array_Core_0_Isoma_ahp_tau;
const int _num__array_Core_0_Isoma_ahp_tau = 1;
double * _array_Core_0_Isoma_ahp_th;
const int _num__array_Core_0_Isoma_ahp_th = 1;
double * _array_Core_0_Isoma_ahp_w;
const int _num__array_Core_0_Isoma_ahp_w = 256;
double * _array_Core_0_Isoma_const;
const int _num__array_Core_0_Isoma_const = 256;
double * _array_Core_0_Isoma_dpi_tau;
const int _num__array_Core_0_Isoma_dpi_tau = 1;
double * _array_Core_0_Isoma_mem;
const int _num__array_Core_0_Isoma_mem = 256;
double * _array_Core_0_Isoma_pfb_gain;
const int _num__array_Core_0_Isoma_pfb_gain = 1;
double * _array_Core_0_Isoma_pfb_norm;
const int _num__array_Core_0_Isoma_pfb_norm = 1;
double * _array_Core_0_Isoma_pfb_th;
const int _num__array_Core_0_Isoma_pfb_th = 1;
double * _array_Core_0_Isoma_reset;
const int _num__array_Core_0_Isoma_reset = 1;
double * _array_Core_0_Isoma_th;
const int _num__array_Core_0_Isoma_th = 1;
double * _array_Core_0_kn;
const int _num__array_Core_0_kn = 1;
double * _array_Core_0_kp;
const int _num__array_Core_0_kp = 1;
double * _array_Core_0_lastspike;
const int _num__array_Core_0_lastspike = 256;
char * _array_Core_0_not_refractory;
const int _num__array_Core_0_not_refractory = 256;
double * _array_Core_0_soma_refP;
const int _num__array_Core_0_soma_refP = 1;
double * _array_Core_0_Ut;
const int _num__array_Core_0_Ut = 1;
double * _array_Core_0_Vnmda;
const int _num__array_Core_0_Vnmda = 256;
int32_t * _array_Core_1__spikespace;
const int _num__array_Core_1__spikespace = 257;
double * _array_Core_1_alpha;
const int _num__array_Core_1_alpha = 1;
double * _array_Core_1_Campa;
const int _num__array_Core_1_Campa = 256;
double * _array_Core_1_Cgaba_a;
const int _num__array_Core_1_Cgaba_a = 256;
double * _array_Core_1_Cgaba_b;
const int _num__array_Core_1_Cgaba_b = 256;
double * _array_Core_1_Cnmda;
const int _num__array_Core_1_Cnmda = 256;
double * _array_Core_1_Csoma_ahp;
const int _num__array_Core_1_Csoma_ahp = 1;
double * _array_Core_1_Csoma_mem;
const int _num__array_Core_1_Csoma_mem = 1;
int32_t * _array_Core_1_i;
const int _num__array_Core_1_i = 256;
double * _array_Core_1_I0;
const int _num__array_Core_1_I0 = 1;
double * _array_Core_1_Iampa;
const int _num__array_Core_1_Iampa = 256;
double * _array_Core_1_Iampa_tau;
const int _num__array_Core_1_Iampa_tau = 256;
double * _array_Core_1_Iampa_w0;
const int _num__array_Core_1_Iampa_w0 = 256;
double * _array_Core_1_Igaba_a;
const int _num__array_Core_1_Igaba_a = 256;
double * _array_Core_1_Igaba_a_tau;
const int _num__array_Core_1_Igaba_a_tau = 256;
double * _array_Core_1_Igaba_a_w0;
const int _num__array_Core_1_Igaba_a_w0 = 256;
double * _array_Core_1_Igaba_b;
const int _num__array_Core_1_Igaba_b = 256;
double * _array_Core_1_Igaba_b_tau;
const int _num__array_Core_1_Igaba_b_tau = 256;
double * _array_Core_1_Igaba_b_w0;
const int _num__array_Core_1_Igaba_b_w0 = 256;
double * _array_Core_1_Inmda;
const int _num__array_Core_1_Inmda = 256;
double * _array_Core_1_Inmda_tau;
const int _num__array_Core_1_Inmda_tau = 256;
double * _array_Core_1_Inmda_w0;
const int _num__array_Core_1_Inmda_w0 = 256;
double * _array_Core_1_Isoma_ahp;
const int _num__array_Core_1_Isoma_ahp = 256;
double * _array_Core_1_Isoma_ahp_tau;
const int _num__array_Core_1_Isoma_ahp_tau = 1;
double * _array_Core_1_Isoma_ahp_th;
const int _num__array_Core_1_Isoma_ahp_th = 1;
double * _array_Core_1_Isoma_ahp_w;
const int _num__array_Core_1_Isoma_ahp_w = 256;
double * _array_Core_1_Isoma_const;
const int _num__array_Core_1_Isoma_const = 256;
double * _array_Core_1_Isoma_dpi_tau;
const int _num__array_Core_1_Isoma_dpi_tau = 1;
double * _array_Core_1_Isoma_mem;
const int _num__array_Core_1_Isoma_mem = 256;
double * _array_Core_1_Isoma_pfb_gain;
const int _num__array_Core_1_Isoma_pfb_gain = 1;
double * _array_Core_1_Isoma_pfb_norm;
const int _num__array_Core_1_Isoma_pfb_norm = 1;
double * _array_Core_1_Isoma_pfb_th;
const int _num__array_Core_1_Isoma_pfb_th = 1;
double * _array_Core_1_Isoma_reset;
const int _num__array_Core_1_Isoma_reset = 1;
double * _array_Core_1_Isoma_th;
const int _num__array_Core_1_Isoma_th = 1;
double * _array_Core_1_kn;
const int _num__array_Core_1_kn = 1;
double * _array_Core_1_kp;
const int _num__array_Core_1_kp = 1;
double * _array_Core_1_lastspike;
const int _num__array_Core_1_lastspike = 256;
char * _array_Core_1_not_refractory;
const int _num__array_Core_1_not_refractory = 256;
double * _array_Core_1_soma_refP;
const int _num__array_Core_1_soma_refP = 1;
int32_t * _array_Core_1_subgroup__sub_idx;
const int _num__array_Core_1_subgroup__sub_idx = 200;
double * _array_Core_1_Ut;
const int _num__array_Core_1_Ut = 1;
double * _array_Core_1_Vnmda;
const int _num__array_Core_1_Vnmda = 256;
int32_t * _array_Core_2__spikespace;
const int _num__array_Core_2__spikespace = 257;
double * _array_Core_2_alpha;
const int _num__array_Core_2_alpha = 1;
double * _array_Core_2_Campa;
const int _num__array_Core_2_Campa = 256;
double * _array_Core_2_Cgaba_a;
const int _num__array_Core_2_Cgaba_a = 256;
double * _array_Core_2_Cgaba_b;
const int _num__array_Core_2_Cgaba_b = 256;
double * _array_Core_2_Cnmda;
const int _num__array_Core_2_Cnmda = 256;
double * _array_Core_2_Csoma_ahp;
const int _num__array_Core_2_Csoma_ahp = 1;
double * _array_Core_2_Csoma_mem;
const int _num__array_Core_2_Csoma_mem = 1;
int32_t * _array_Core_2_i;
const int _num__array_Core_2_i = 256;
double * _array_Core_2_I0;
const int _num__array_Core_2_I0 = 1;
double * _array_Core_2_Iampa;
const int _num__array_Core_2_Iampa = 256;
double * _array_Core_2_Iampa_tau;
const int _num__array_Core_2_Iampa_tau = 256;
double * _array_Core_2_Iampa_w0;
const int _num__array_Core_2_Iampa_w0 = 256;
double * _array_Core_2_Igaba_a;
const int _num__array_Core_2_Igaba_a = 256;
double * _array_Core_2_Igaba_a_tau;
const int _num__array_Core_2_Igaba_a_tau = 256;
double * _array_Core_2_Igaba_a_w0;
const int _num__array_Core_2_Igaba_a_w0 = 256;
double * _array_Core_2_Igaba_b;
const int _num__array_Core_2_Igaba_b = 256;
double * _array_Core_2_Igaba_b_tau;
const int _num__array_Core_2_Igaba_b_tau = 256;
double * _array_Core_2_Igaba_b_w0;
const int _num__array_Core_2_Igaba_b_w0 = 256;
double * _array_Core_2_Inmda;
const int _num__array_Core_2_Inmda = 256;
double * _array_Core_2_Inmda_tau;
const int _num__array_Core_2_Inmda_tau = 256;
double * _array_Core_2_Inmda_w0;
const int _num__array_Core_2_Inmda_w0 = 256;
double * _array_Core_2_Isoma_ahp;
const int _num__array_Core_2_Isoma_ahp = 256;
double * _array_Core_2_Isoma_ahp_tau;
const int _num__array_Core_2_Isoma_ahp_tau = 1;
double * _array_Core_2_Isoma_ahp_th;
const int _num__array_Core_2_Isoma_ahp_th = 1;
double * _array_Core_2_Isoma_ahp_w;
const int _num__array_Core_2_Isoma_ahp_w = 256;
double * _array_Core_2_Isoma_const;
const int _num__array_Core_2_Isoma_const = 256;
double * _array_Core_2_Isoma_dpi_tau;
const int _num__array_Core_2_Isoma_dpi_tau = 1;
double * _array_Core_2_Isoma_mem;
const int _num__array_Core_2_Isoma_mem = 256;
double * _array_Core_2_Isoma_pfb_gain;
const int _num__array_Core_2_Isoma_pfb_gain = 1;
double * _array_Core_2_Isoma_pfb_norm;
const int _num__array_Core_2_Isoma_pfb_norm = 1;
double * _array_Core_2_Isoma_pfb_th;
const int _num__array_Core_2_Isoma_pfb_th = 1;
double * _array_Core_2_Isoma_reset;
const int _num__array_Core_2_Isoma_reset = 1;
double * _array_Core_2_Isoma_th;
const int _num__array_Core_2_Isoma_th = 1;
double * _array_Core_2_kn;
const int _num__array_Core_2_kn = 1;
double * _array_Core_2_kp;
const int _num__array_Core_2_kp = 1;
double * _array_Core_2_lastspike;
const int _num__array_Core_2_lastspike = 256;
char * _array_Core_2_not_refractory;
const int _num__array_Core_2_not_refractory = 256;
double * _array_Core_2_soma_refP;
const int _num__array_Core_2_soma_refP = 1;
double * _array_Core_2_Ut;
const int _num__array_Core_2_Ut = 1;
double * _array_Core_2_Vnmda;
const int _num__array_Core_2_Vnmda = 256;
int32_t * _array_Core_3__spikespace;
const int _num__array_Core_3__spikespace = 257;
double * _array_Core_3_alpha;
const int _num__array_Core_3_alpha = 1;
double * _array_Core_3_Campa;
const int _num__array_Core_3_Campa = 256;
double * _array_Core_3_Cgaba_a;
const int _num__array_Core_3_Cgaba_a = 256;
double * _array_Core_3_Cgaba_b;
const int _num__array_Core_3_Cgaba_b = 256;
double * _array_Core_3_Cnmda;
const int _num__array_Core_3_Cnmda = 256;
double * _array_Core_3_Csoma_ahp;
const int _num__array_Core_3_Csoma_ahp = 1;
double * _array_Core_3_Csoma_mem;
const int _num__array_Core_3_Csoma_mem = 1;
int32_t * _array_Core_3_i;
const int _num__array_Core_3_i = 256;
double * _array_Core_3_I0;
const int _num__array_Core_3_I0 = 1;
double * _array_Core_3_Iampa;
const int _num__array_Core_3_Iampa = 256;
double * _array_Core_3_Iampa_tau;
const int _num__array_Core_3_Iampa_tau = 256;
double * _array_Core_3_Iampa_w0;
const int _num__array_Core_3_Iampa_w0 = 256;
double * _array_Core_3_Igaba_a;
const int _num__array_Core_3_Igaba_a = 256;
double * _array_Core_3_Igaba_a_tau;
const int _num__array_Core_3_Igaba_a_tau = 256;
double * _array_Core_3_Igaba_a_w0;
const int _num__array_Core_3_Igaba_a_w0 = 256;
double * _array_Core_3_Igaba_b;
const int _num__array_Core_3_Igaba_b = 256;
double * _array_Core_3_Igaba_b_tau;
const int _num__array_Core_3_Igaba_b_tau = 256;
double * _array_Core_3_Igaba_b_w0;
const int _num__array_Core_3_Igaba_b_w0 = 256;
double * _array_Core_3_Inmda;
const int _num__array_Core_3_Inmda = 256;
double * _array_Core_3_Inmda_tau;
const int _num__array_Core_3_Inmda_tau = 256;
double * _array_Core_3_Inmda_w0;
const int _num__array_Core_3_Inmda_w0 = 256;
double * _array_Core_3_Isoma_ahp;
const int _num__array_Core_3_Isoma_ahp = 256;
double * _array_Core_3_Isoma_ahp_tau;
const int _num__array_Core_3_Isoma_ahp_tau = 1;
double * _array_Core_3_Isoma_ahp_th;
const int _num__array_Core_3_Isoma_ahp_th = 1;
double * _array_Core_3_Isoma_ahp_w;
const int _num__array_Core_3_Isoma_ahp_w = 256;
double * _array_Core_3_Isoma_const;
const int _num__array_Core_3_Isoma_const = 256;
double * _array_Core_3_Isoma_dpi_tau;
const int _num__array_Core_3_Isoma_dpi_tau = 1;
double * _array_Core_3_Isoma_mem;
const int _num__array_Core_3_Isoma_mem = 256;
double * _array_Core_3_Isoma_pfb_gain;
const int _num__array_Core_3_Isoma_pfb_gain = 1;
double * _array_Core_3_Isoma_pfb_norm;
const int _num__array_Core_3_Isoma_pfb_norm = 1;
double * _array_Core_3_Isoma_pfb_th;
const int _num__array_Core_3_Isoma_pfb_th = 1;
double * _array_Core_3_Isoma_reset;
const int _num__array_Core_3_Isoma_reset = 1;
double * _array_Core_3_Isoma_th;
const int _num__array_Core_3_Isoma_th = 1;
double * _array_Core_3_kn;
const int _num__array_Core_3_kn = 1;
double * _array_Core_3_kp;
const int _num__array_Core_3_kp = 1;
double * _array_Core_3_lastspike;
const int _num__array_Core_3_lastspike = 256;
char * _array_Core_3_not_refractory;
const int _num__array_Core_3_not_refractory = 256;
double * _array_Core_3_soma_refP;
const int _num__array_Core_3_soma_refP = 1;
double * _array_Core_3_Ut;
const int _num__array_Core_3_Ut = 1;
double * _array_Core_3_Vnmda;
const int _num__array_Core_3_Vnmda = 256;
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
int64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_mon_neuron_input__source_idx;
const int _num__array_mon_neuron_input__source_idx = 200;
int32_t * _array_mon_neuron_input_count;
const int _num__array_mon_neuron_input_count = 200;
int32_t * _array_mon_neuron_input_N;
const int _num__array_mon_neuron_input_N = 1;
int32_t * _array_mon_neuron_output__source_idx;
const int _num__array_mon_neuron_output__source_idx = 200;
int32_t * _array_mon_neuron_output_count;
const int _num__array_mon_neuron_output_count = 200;
int32_t * _array_mon_neuron_output_N;
const int _num__array_mon_neuron_output_N = 1;
int32_t * _array_NMDA0_N;
const int _num__array_NMDA0_N = 1;
int32_t * _array_poissongroup__spikespace;
const int _num__array_poissongroup__spikespace = 201;
int32_t * _array_poissongroup_i;
const int _num__array_poissongroup_i = 200;

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_mon_neuron_input_i;
std::vector<double> _dynamic_array_mon_neuron_input_t;
std::vector<int32_t> _dynamic_array_mon_neuron_output_i;
std::vector<double> _dynamic_array_mon_neuron_output_t;
std::vector<int32_t> _dynamic_array_NMDA0__synaptic_post;
std::vector<int32_t> _dynamic_array_NMDA0__synaptic_pre;
std::vector<double> _dynamic_array_NMDA0_delay;
std::vector<double> _dynamic_array_NMDA0_delay_1;
std::vector<int32_t> _dynamic_array_NMDA0_N_incoming;
std::vector<int32_t> _dynamic_array_NMDA0_N_outgoing;
std::vector<double> _dynamic_array_NMDA0_weight;

//////////////// dynamic arrays 2d /////////

/////////////// static arrays /////////////
double * _timedarray_values;
const int _num__timedarray_values = 1000000;

//////////////// synapses /////////////////
// NMDA0
SynapticPathway NMDA0_post(
		_dynamic_array_NMDA0__synaptic_post,
		0, 200);
SynapticPathway NMDA0_pre(
		_dynamic_array_NMDA0__synaptic_pre,
		0, 200);

//////////////// clocks ///////////////////
Clock defaultclock;  // attributes will be set in run.cpp

// Profiling information for each code object
}

void _init_arrays()
{
	using namespace brian;

    // Arrays initialized to 0
	_array_Core_0__spikespace = new int32_t[257];
    
	for(int i=0; i<257; i++) _array_Core_0__spikespace[i] = 0;

	_array_Core_0_alpha = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_alpha[i] = 0;

	_array_Core_0_Campa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Campa[i] = 0;

	_array_Core_0_Cgaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Cgaba_a[i] = 0;

	_array_Core_0_Cgaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Cgaba_b[i] = 0;

	_array_Core_0_Cnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Cnmda[i] = 0;

	_array_Core_0_Csoma_ahp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Csoma_ahp[i] = 0;

	_array_Core_0_Csoma_mem = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Csoma_mem[i] = 0;

	_array_Core_0_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_0_i[i] = 0;

	_array_Core_0_I0 = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_I0[i] = 0;

	_array_Core_0_Iampa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Iampa[i] = 0;

	_array_Core_0_Iampa_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Iampa_tau[i] = 0;

	_array_Core_0_Iampa_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Iampa_w0[i] = 0;

	_array_Core_0_Igaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_a[i] = 0;

	_array_Core_0_Igaba_a_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_a_tau[i] = 0;

	_array_Core_0_Igaba_a_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_a_w0[i] = 0;

	_array_Core_0_Igaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_b[i] = 0;

	_array_Core_0_Igaba_b_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_b_tau[i] = 0;

	_array_Core_0_Igaba_b_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Igaba_b_w0[i] = 0;

	_array_Core_0_Inmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Inmda[i] = 0;

	_array_Core_0_Inmda_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Inmda_tau[i] = 0;

	_array_Core_0_Inmda_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Inmda_w0[i] = 0;

	_array_Core_0_Isoma_ahp = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Isoma_ahp[i] = 0;

	_array_Core_0_Isoma_ahp_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_ahp_tau[i] = 0;

	_array_Core_0_Isoma_ahp_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_ahp_th[i] = 0;

	_array_Core_0_Isoma_ahp_w = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Isoma_ahp_w[i] = 0;

	_array_Core_0_Isoma_const = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Isoma_const[i] = 0;

	_array_Core_0_Isoma_dpi_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_dpi_tau[i] = 0;

	_array_Core_0_Isoma_mem = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Isoma_mem[i] = 0;

	_array_Core_0_Isoma_pfb_gain = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_pfb_gain[i] = 0;

	_array_Core_0_Isoma_pfb_norm = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_pfb_norm[i] = 0;

	_array_Core_0_Isoma_pfb_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_pfb_th[i] = 0;

	_array_Core_0_Isoma_reset = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_reset[i] = 0;

	_array_Core_0_Isoma_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Isoma_th[i] = 0;

	_array_Core_0_kn = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_kn[i] = 0;

	_array_Core_0_kp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_kp[i] = 0;

	_array_Core_0_lastspike = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_lastspike[i] = 0;

	_array_Core_0_not_refractory = new char[256];
    
	for(int i=0; i<256; i++) _array_Core_0_not_refractory[i] = 0;

	_array_Core_0_soma_refP = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_soma_refP[i] = 0;

	_array_Core_0_Ut = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_0_Ut[i] = 0;

	_array_Core_0_Vnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_0_Vnmda[i] = 0;

	_array_Core_1__spikespace = new int32_t[257];
    
	for(int i=0; i<257; i++) _array_Core_1__spikespace[i] = 0;

	_array_Core_1_alpha = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_alpha[i] = 0;

	_array_Core_1_Campa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Campa[i] = 0;

	_array_Core_1_Cgaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Cgaba_a[i] = 0;

	_array_Core_1_Cgaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Cgaba_b[i] = 0;

	_array_Core_1_Cnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Cnmda[i] = 0;

	_array_Core_1_Csoma_ahp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Csoma_ahp[i] = 0;

	_array_Core_1_Csoma_mem = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Csoma_mem[i] = 0;

	_array_Core_1_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_1_i[i] = 0;

	_array_Core_1_I0 = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_I0[i] = 0;

	_array_Core_1_Iampa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Iampa[i] = 0;

	_array_Core_1_Iampa_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Iampa_tau[i] = 0;

	_array_Core_1_Iampa_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Iampa_w0[i] = 0;

	_array_Core_1_Igaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_a[i] = 0;

	_array_Core_1_Igaba_a_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_a_tau[i] = 0;

	_array_Core_1_Igaba_a_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_a_w0[i] = 0;

	_array_Core_1_Igaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_b[i] = 0;

	_array_Core_1_Igaba_b_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_b_tau[i] = 0;

	_array_Core_1_Igaba_b_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Igaba_b_w0[i] = 0;

	_array_Core_1_Inmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Inmda[i] = 0;

	_array_Core_1_Inmda_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Inmda_tau[i] = 0;

	_array_Core_1_Inmda_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Inmda_w0[i] = 0;

	_array_Core_1_Isoma_ahp = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Isoma_ahp[i] = 0;

	_array_Core_1_Isoma_ahp_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_ahp_tau[i] = 0;

	_array_Core_1_Isoma_ahp_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_ahp_th[i] = 0;

	_array_Core_1_Isoma_ahp_w = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Isoma_ahp_w[i] = 0;

	_array_Core_1_Isoma_const = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Isoma_const[i] = 0;

	_array_Core_1_Isoma_dpi_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_dpi_tau[i] = 0;

	_array_Core_1_Isoma_mem = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Isoma_mem[i] = 0;

	_array_Core_1_Isoma_pfb_gain = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_pfb_gain[i] = 0;

	_array_Core_1_Isoma_pfb_norm = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_pfb_norm[i] = 0;

	_array_Core_1_Isoma_pfb_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_pfb_th[i] = 0;

	_array_Core_1_Isoma_reset = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_reset[i] = 0;

	_array_Core_1_Isoma_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Isoma_th[i] = 0;

	_array_Core_1_kn = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_kn[i] = 0;

	_array_Core_1_kp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_kp[i] = 0;

	_array_Core_1_lastspike = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_lastspike[i] = 0;

	_array_Core_1_not_refractory = new char[256];
    
	for(int i=0; i<256; i++) _array_Core_1_not_refractory[i] = 0;

	_array_Core_1_soma_refP = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_soma_refP[i] = 0;

	_array_Core_1_subgroup__sub_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_Core_1_subgroup__sub_idx[i] = 0;

	_array_Core_1_Ut = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_1_Ut[i] = 0;

	_array_Core_1_Vnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_1_Vnmda[i] = 0;

	_array_Core_2__spikespace = new int32_t[257];
    
	for(int i=0; i<257; i++) _array_Core_2__spikespace[i] = 0;

	_array_Core_2_alpha = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_alpha[i] = 0;

	_array_Core_2_Campa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Campa[i] = 0;

	_array_Core_2_Cgaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Cgaba_a[i] = 0;

	_array_Core_2_Cgaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Cgaba_b[i] = 0;

	_array_Core_2_Cnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Cnmda[i] = 0;

	_array_Core_2_Csoma_ahp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Csoma_ahp[i] = 0;

	_array_Core_2_Csoma_mem = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Csoma_mem[i] = 0;

	_array_Core_2_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_2_i[i] = 0;

	_array_Core_2_I0 = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_I0[i] = 0;

	_array_Core_2_Iampa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Iampa[i] = 0;

	_array_Core_2_Iampa_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Iampa_tau[i] = 0;

	_array_Core_2_Iampa_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Iampa_w0[i] = 0;

	_array_Core_2_Igaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_a[i] = 0;

	_array_Core_2_Igaba_a_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_a_tau[i] = 0;

	_array_Core_2_Igaba_a_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_a_w0[i] = 0;

	_array_Core_2_Igaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_b[i] = 0;

	_array_Core_2_Igaba_b_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_b_tau[i] = 0;

	_array_Core_2_Igaba_b_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Igaba_b_w0[i] = 0;

	_array_Core_2_Inmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Inmda[i] = 0;

	_array_Core_2_Inmda_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Inmda_tau[i] = 0;

	_array_Core_2_Inmda_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Inmda_w0[i] = 0;

	_array_Core_2_Isoma_ahp = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Isoma_ahp[i] = 0;

	_array_Core_2_Isoma_ahp_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_ahp_tau[i] = 0;

	_array_Core_2_Isoma_ahp_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_ahp_th[i] = 0;

	_array_Core_2_Isoma_ahp_w = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Isoma_ahp_w[i] = 0;

	_array_Core_2_Isoma_const = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Isoma_const[i] = 0;

	_array_Core_2_Isoma_dpi_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_dpi_tau[i] = 0;

	_array_Core_2_Isoma_mem = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Isoma_mem[i] = 0;

	_array_Core_2_Isoma_pfb_gain = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_pfb_gain[i] = 0;

	_array_Core_2_Isoma_pfb_norm = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_pfb_norm[i] = 0;

	_array_Core_2_Isoma_pfb_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_pfb_th[i] = 0;

	_array_Core_2_Isoma_reset = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_reset[i] = 0;

	_array_Core_2_Isoma_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Isoma_th[i] = 0;

	_array_Core_2_kn = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_kn[i] = 0;

	_array_Core_2_kp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_kp[i] = 0;

	_array_Core_2_lastspike = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_lastspike[i] = 0;

	_array_Core_2_not_refractory = new char[256];
    
	for(int i=0; i<256; i++) _array_Core_2_not_refractory[i] = 0;

	_array_Core_2_soma_refP = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_soma_refP[i] = 0;

	_array_Core_2_Ut = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_2_Ut[i] = 0;

	_array_Core_2_Vnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_2_Vnmda[i] = 0;

	_array_Core_3__spikespace = new int32_t[257];
    
	for(int i=0; i<257; i++) _array_Core_3__spikespace[i] = 0;

	_array_Core_3_alpha = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_alpha[i] = 0;

	_array_Core_3_Campa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Campa[i] = 0;

	_array_Core_3_Cgaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Cgaba_a[i] = 0;

	_array_Core_3_Cgaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Cgaba_b[i] = 0;

	_array_Core_3_Cnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Cnmda[i] = 0;

	_array_Core_3_Csoma_ahp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Csoma_ahp[i] = 0;

	_array_Core_3_Csoma_mem = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Csoma_mem[i] = 0;

	_array_Core_3_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_3_i[i] = 0;

	_array_Core_3_I0 = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_I0[i] = 0;

	_array_Core_3_Iampa = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Iampa[i] = 0;

	_array_Core_3_Iampa_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Iampa_tau[i] = 0;

	_array_Core_3_Iampa_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Iampa_w0[i] = 0;

	_array_Core_3_Igaba_a = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_a[i] = 0;

	_array_Core_3_Igaba_a_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_a_tau[i] = 0;

	_array_Core_3_Igaba_a_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_a_w0[i] = 0;

	_array_Core_3_Igaba_b = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_b[i] = 0;

	_array_Core_3_Igaba_b_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_b_tau[i] = 0;

	_array_Core_3_Igaba_b_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Igaba_b_w0[i] = 0;

	_array_Core_3_Inmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Inmda[i] = 0;

	_array_Core_3_Inmda_tau = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Inmda_tau[i] = 0;

	_array_Core_3_Inmda_w0 = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Inmda_w0[i] = 0;

	_array_Core_3_Isoma_ahp = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Isoma_ahp[i] = 0;

	_array_Core_3_Isoma_ahp_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_ahp_tau[i] = 0;

	_array_Core_3_Isoma_ahp_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_ahp_th[i] = 0;

	_array_Core_3_Isoma_ahp_w = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Isoma_ahp_w[i] = 0;

	_array_Core_3_Isoma_const = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Isoma_const[i] = 0;

	_array_Core_3_Isoma_dpi_tau = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_dpi_tau[i] = 0;

	_array_Core_3_Isoma_mem = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Isoma_mem[i] = 0;

	_array_Core_3_Isoma_pfb_gain = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_pfb_gain[i] = 0;

	_array_Core_3_Isoma_pfb_norm = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_pfb_norm[i] = 0;

	_array_Core_3_Isoma_pfb_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_pfb_th[i] = 0;

	_array_Core_3_Isoma_reset = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_reset[i] = 0;

	_array_Core_3_Isoma_th = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Isoma_th[i] = 0;

	_array_Core_3_kn = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_kn[i] = 0;

	_array_Core_3_kp = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_kp[i] = 0;

	_array_Core_3_lastspike = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_lastspike[i] = 0;

	_array_Core_3_not_refractory = new char[256];
    
	for(int i=0; i<256; i++) _array_Core_3_not_refractory[i] = 0;

	_array_Core_3_soma_refP = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_soma_refP[i] = 0;

	_array_Core_3_Ut = new double[1];
    
	for(int i=0; i<1; i++) _array_Core_3_Ut[i] = 0;

	_array_Core_3_Vnmda = new double[256];
    
	for(int i=0; i<256; i++) _array_Core_3_Vnmda[i] = 0;

	_array_defaultclock_dt = new double[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;

	_array_defaultclock_t = new double[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;

	_array_defaultclock_timestep = new int64_t[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;

	_array_mon_neuron_input__source_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_input__source_idx[i] = 0;

	_array_mon_neuron_input_count = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_input_count[i] = 0;

	_array_mon_neuron_input_N = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_mon_neuron_input_N[i] = 0;

	_array_mon_neuron_output__source_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_output__source_idx[i] = 0;

	_array_mon_neuron_output_count = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_output_count[i] = 0;

	_array_mon_neuron_output_N = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_mon_neuron_output_N[i] = 0;

	_array_NMDA0_N = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_NMDA0_N[i] = 0;

	_array_poissongroup__spikespace = new int32_t[201];
    
	for(int i=0; i<201; i++) _array_poissongroup__spikespace[i] = 0;

	_array_poissongroup_i = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_poissongroup_i[i] = 0;


	// Arrays initialized to an "arange"
	_array_Core_0_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_0_i[i] = 0 + i;

	_array_Core_1_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_1_i[i] = 0 + i;

	_array_Core_1_subgroup__sub_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_Core_1_subgroup__sub_idx[i] = 0 + i;

	_array_Core_2_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_2_i[i] = 0 + i;

	_array_Core_3_i = new int32_t[256];
    
	for(int i=0; i<256; i++) _array_Core_3_i[i] = 0 + i;

	_array_mon_neuron_input__source_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_input__source_idx[i] = 0 + i;

	_array_mon_neuron_output__source_idx = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_mon_neuron_output__source_idx[i] = 0 + i;

	_array_poissongroup_i = new int32_t[200];
    
	for(int i=0; i<200; i++) _array_poissongroup_i[i] = 0 + i;


	// static arrays
	_timedarray_values = new double[1000000];

	// Random number generator states
	for (int i=0; i<1; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_timedarray_values;
	f_timedarray_values.open("static_arrays/_timedarray_values", ios::in | ios::binary);
	if(f_timedarray_values.is_open())
	{
		f_timedarray_values.read(reinterpret_cast<char*>(_timedarray_values), 1000000*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_values." << endl;
	}
}

void _write_arrays()
{
	using namespace brian;

	ofstream outfile__array_Core_0__spikespace;
	outfile__array_Core_0__spikespace.open("results/_array_Core_0__spikespace_8097146419878060048", ios::binary | ios::out);
	if(outfile__array_Core_0__spikespace.is_open())
	{
		outfile__array_Core_0__spikespace.write(reinterpret_cast<char*>(_array_Core_0__spikespace), 257*sizeof(_array_Core_0__spikespace[0]));
		outfile__array_Core_0__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0__spikespace." << endl;
	}
	ofstream outfile__array_Core_0_alpha;
	outfile__array_Core_0_alpha.open("results/_array_Core_0_alpha_6656631157718018654", ios::binary | ios::out);
	if(outfile__array_Core_0_alpha.is_open())
	{
		outfile__array_Core_0_alpha.write(reinterpret_cast<char*>(_array_Core_0_alpha), 1*sizeof(_array_Core_0_alpha[0]));
		outfile__array_Core_0_alpha.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_alpha." << endl;
	}
	ofstream outfile__array_Core_0_Campa;
	outfile__array_Core_0_Campa.open("results/_array_Core_0_Campa_4450902435179399783", ios::binary | ios::out);
	if(outfile__array_Core_0_Campa.is_open())
	{
		outfile__array_Core_0_Campa.write(reinterpret_cast<char*>(_array_Core_0_Campa), 256*sizeof(_array_Core_0_Campa[0]));
		outfile__array_Core_0_Campa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Campa." << endl;
	}
	ofstream outfile__array_Core_0_Cgaba_a;
	outfile__array_Core_0_Cgaba_a.open("results/_array_Core_0_Cgaba_a_-917725417654597309", ios::binary | ios::out);
	if(outfile__array_Core_0_Cgaba_a.is_open())
	{
		outfile__array_Core_0_Cgaba_a.write(reinterpret_cast<char*>(_array_Core_0_Cgaba_a), 256*sizeof(_array_Core_0_Cgaba_a[0]));
		outfile__array_Core_0_Cgaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Cgaba_a." << endl;
	}
	ofstream outfile__array_Core_0_Cgaba_b;
	outfile__array_Core_0_Cgaba_b.open("results/_array_Core_0_Cgaba_b_-3843185481736145924", ios::binary | ios::out);
	if(outfile__array_Core_0_Cgaba_b.is_open())
	{
		outfile__array_Core_0_Cgaba_b.write(reinterpret_cast<char*>(_array_Core_0_Cgaba_b), 256*sizeof(_array_Core_0_Cgaba_b[0]));
		outfile__array_Core_0_Cgaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Cgaba_b." << endl;
	}
	ofstream outfile__array_Core_0_Cnmda;
	outfile__array_Core_0_Cnmda.open("results/_array_Core_0_Cnmda_-9010445725137049026", ios::binary | ios::out);
	if(outfile__array_Core_0_Cnmda.is_open())
	{
		outfile__array_Core_0_Cnmda.write(reinterpret_cast<char*>(_array_Core_0_Cnmda), 256*sizeof(_array_Core_0_Cnmda[0]));
		outfile__array_Core_0_Cnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Cnmda." << endl;
	}
	ofstream outfile__array_Core_0_Csoma_ahp;
	outfile__array_Core_0_Csoma_ahp.open("results/_array_Core_0_Csoma_ahp_-3937814710898105789", ios::binary | ios::out);
	if(outfile__array_Core_0_Csoma_ahp.is_open())
	{
		outfile__array_Core_0_Csoma_ahp.write(reinterpret_cast<char*>(_array_Core_0_Csoma_ahp), 1*sizeof(_array_Core_0_Csoma_ahp[0]));
		outfile__array_Core_0_Csoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Csoma_ahp." << endl;
	}
	ofstream outfile__array_Core_0_Csoma_mem;
	outfile__array_Core_0_Csoma_mem.open("results/_array_Core_0_Csoma_mem_5783530513863066554", ios::binary | ios::out);
	if(outfile__array_Core_0_Csoma_mem.is_open())
	{
		outfile__array_Core_0_Csoma_mem.write(reinterpret_cast<char*>(_array_Core_0_Csoma_mem), 1*sizeof(_array_Core_0_Csoma_mem[0]));
		outfile__array_Core_0_Csoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Csoma_mem." << endl;
	}
	ofstream outfile__array_Core_0_i;
	outfile__array_Core_0_i.open("results/_array_Core_0_i_-5960520545996854005", ios::binary | ios::out);
	if(outfile__array_Core_0_i.is_open())
	{
		outfile__array_Core_0_i.write(reinterpret_cast<char*>(_array_Core_0_i), 256*sizeof(_array_Core_0_i[0]));
		outfile__array_Core_0_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_i." << endl;
	}
	ofstream outfile__array_Core_0_I0;
	outfile__array_Core_0_I0.open("results/_array_Core_0_I0_-3021259697338811357", ios::binary | ios::out);
	if(outfile__array_Core_0_I0.is_open())
	{
		outfile__array_Core_0_I0.write(reinterpret_cast<char*>(_array_Core_0_I0), 1*sizeof(_array_Core_0_I0[0]));
		outfile__array_Core_0_I0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_I0." << endl;
	}
	ofstream outfile__array_Core_0_Iampa;
	outfile__array_Core_0_Iampa.open("results/_array_Core_0_Iampa_-1457427188187744816", ios::binary | ios::out);
	if(outfile__array_Core_0_Iampa.is_open())
	{
		outfile__array_Core_0_Iampa.write(reinterpret_cast<char*>(_array_Core_0_Iampa), 256*sizeof(_array_Core_0_Iampa[0]));
		outfile__array_Core_0_Iampa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Iampa." << endl;
	}
	ofstream outfile__array_Core_0_Iampa_tau;
	outfile__array_Core_0_Iampa_tau.open("results/_array_Core_0_Iampa_tau_6037061684352902888", ios::binary | ios::out);
	if(outfile__array_Core_0_Iampa_tau.is_open())
	{
		outfile__array_Core_0_Iampa_tau.write(reinterpret_cast<char*>(_array_Core_0_Iampa_tau), 256*sizeof(_array_Core_0_Iampa_tau[0]));
		outfile__array_Core_0_Iampa_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Iampa_tau." << endl;
	}
	ofstream outfile__array_Core_0_Iampa_w0;
	outfile__array_Core_0_Iampa_w0.open("results/_array_Core_0_Iampa_w0_-3347142568882949289", ios::binary | ios::out);
	if(outfile__array_Core_0_Iampa_w0.is_open())
	{
		outfile__array_Core_0_Iampa_w0.write(reinterpret_cast<char*>(_array_Core_0_Iampa_w0), 256*sizeof(_array_Core_0_Iampa_w0[0]));
		outfile__array_Core_0_Iampa_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Iampa_w0." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_a;
	outfile__array_Core_0_Igaba_a.open("results/_array_Core_0_Igaba_a_5102956005768905199", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_a.is_open())
	{
		outfile__array_Core_0_Igaba_a.write(reinterpret_cast<char*>(_array_Core_0_Igaba_a), 256*sizeof(_array_Core_0_Igaba_a[0]));
		outfile__array_Core_0_Igaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_a." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_a_tau;
	outfile__array_Core_0_Igaba_a_tau.open("results/_array_Core_0_Igaba_a_tau_-2548481926726133868", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_a_tau.is_open())
	{
		outfile__array_Core_0_Igaba_a_tau.write(reinterpret_cast<char*>(_array_Core_0_Igaba_a_tau), 256*sizeof(_array_Core_0_Igaba_a_tau[0]));
		outfile__array_Core_0_Igaba_a_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_a_tau." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_a_w0;
	outfile__array_Core_0_Igaba_a_w0.open("results/_array_Core_0_Igaba_a_w0_-851216200074112354", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_a_w0.is_open())
	{
		outfile__array_Core_0_Igaba_a_w0.write(reinterpret_cast<char*>(_array_Core_0_Igaba_a_w0), 256*sizeof(_array_Core_0_Igaba_a_w0[0]));
		outfile__array_Core_0_Igaba_a_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_a_w0." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_b;
	outfile__array_Core_0_Igaba_b.open("results/_array_Core_0_Igaba_b_-959612227091446291", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_b.is_open())
	{
		outfile__array_Core_0_Igaba_b.write(reinterpret_cast<char*>(_array_Core_0_Igaba_b), 256*sizeof(_array_Core_0_Igaba_b[0]));
		outfile__array_Core_0_Igaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_b." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_b_tau;
	outfile__array_Core_0_Igaba_b_tau.open("results/_array_Core_0_Igaba_b_tau_-818511749708402440", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_b_tau.is_open())
	{
		outfile__array_Core_0_Igaba_b_tau.write(reinterpret_cast<char*>(_array_Core_0_Igaba_b_tau), 256*sizeof(_array_Core_0_Igaba_b_tau[0]));
		outfile__array_Core_0_Igaba_b_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_b_tau." << endl;
	}
	ofstream outfile__array_Core_0_Igaba_b_w0;
	outfile__array_Core_0_Igaba_b_w0.open("results/_array_Core_0_Igaba_b_w0_3463872194572203552", ios::binary | ios::out);
	if(outfile__array_Core_0_Igaba_b_w0.is_open())
	{
		outfile__array_Core_0_Igaba_b_w0.write(reinterpret_cast<char*>(_array_Core_0_Igaba_b_w0), 256*sizeof(_array_Core_0_Igaba_b_w0[0]));
		outfile__array_Core_0_Igaba_b_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Igaba_b_w0." << endl;
	}
	ofstream outfile__array_Core_0_Inmda;
	outfile__array_Core_0_Inmda.open("results/_array_Core_0_Inmda_6776487861035268665", ios::binary | ios::out);
	if(outfile__array_Core_0_Inmda.is_open())
	{
		outfile__array_Core_0_Inmda.write(reinterpret_cast<char*>(_array_Core_0_Inmda), 256*sizeof(_array_Core_0_Inmda[0]));
		outfile__array_Core_0_Inmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Inmda." << endl;
	}
	ofstream outfile__array_Core_0_Inmda_tau;
	outfile__array_Core_0_Inmda_tau.open("results/_array_Core_0_Inmda_tau_-2637992797932041448", ios::binary | ios::out);
	if(outfile__array_Core_0_Inmda_tau.is_open())
	{
		outfile__array_Core_0_Inmda_tau.write(reinterpret_cast<char*>(_array_Core_0_Inmda_tau), 256*sizeof(_array_Core_0_Inmda_tau[0]));
		outfile__array_Core_0_Inmda_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Inmda_tau." << endl;
	}
	ofstream outfile__array_Core_0_Inmda_w0;
	outfile__array_Core_0_Inmda_w0.open("results/_array_Core_0_Inmda_w0_-8190665152088380923", ios::binary | ios::out);
	if(outfile__array_Core_0_Inmda_w0.is_open())
	{
		outfile__array_Core_0_Inmda_w0.write(reinterpret_cast<char*>(_array_Core_0_Inmda_w0), 256*sizeof(_array_Core_0_Inmda_w0[0]));
		outfile__array_Core_0_Inmda_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Inmda_w0." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_ahp;
	outfile__array_Core_0_Isoma_ahp.open("results/_array_Core_0_Isoma_ahp_-5627077552545085753", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_ahp.is_open())
	{
		outfile__array_Core_0_Isoma_ahp.write(reinterpret_cast<char*>(_array_Core_0_Isoma_ahp), 256*sizeof(_array_Core_0_Isoma_ahp[0]));
		outfile__array_Core_0_Isoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_ahp." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_ahp_tau;
	outfile__array_Core_0_Isoma_ahp_tau.open("results/_array_Core_0_Isoma_ahp_tau_-5049599420396491452", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_ahp_tau.is_open())
	{
		outfile__array_Core_0_Isoma_ahp_tau.write(reinterpret_cast<char*>(_array_Core_0_Isoma_ahp_tau), 1*sizeof(_array_Core_0_Isoma_ahp_tau[0]));
		outfile__array_Core_0_Isoma_ahp_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_ahp_tau." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_ahp_th;
	outfile__array_Core_0_Isoma_ahp_th.open("results/_array_Core_0_Isoma_ahp_th_-8166506817854959084", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_ahp_th.is_open())
	{
		outfile__array_Core_0_Isoma_ahp_th.write(reinterpret_cast<char*>(_array_Core_0_Isoma_ahp_th), 1*sizeof(_array_Core_0_Isoma_ahp_th[0]));
		outfile__array_Core_0_Isoma_ahp_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_ahp_th." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_ahp_w;
	outfile__array_Core_0_Isoma_ahp_w.open("results/_array_Core_0_Isoma_ahp_w_-2967083318276009201", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_ahp_w.is_open())
	{
		outfile__array_Core_0_Isoma_ahp_w.write(reinterpret_cast<char*>(_array_Core_0_Isoma_ahp_w), 256*sizeof(_array_Core_0_Isoma_ahp_w[0]));
		outfile__array_Core_0_Isoma_ahp_w.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_ahp_w." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_const;
	outfile__array_Core_0_Isoma_const.open("results/_array_Core_0_Isoma_const_-7075051387610329213", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_const.is_open())
	{
		outfile__array_Core_0_Isoma_const.write(reinterpret_cast<char*>(_array_Core_0_Isoma_const), 256*sizeof(_array_Core_0_Isoma_const[0]));
		outfile__array_Core_0_Isoma_const.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_const." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_dpi_tau;
	outfile__array_Core_0_Isoma_dpi_tau.open("results/_array_Core_0_Isoma_dpi_tau_-5595998642996866899", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_dpi_tau.is_open())
	{
		outfile__array_Core_0_Isoma_dpi_tau.write(reinterpret_cast<char*>(_array_Core_0_Isoma_dpi_tau), 1*sizeof(_array_Core_0_Isoma_dpi_tau[0]));
		outfile__array_Core_0_Isoma_dpi_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_dpi_tau." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_mem;
	outfile__array_Core_0_Isoma_mem.open("results/_array_Core_0_Isoma_mem_3266650047496024670", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_mem.is_open())
	{
		outfile__array_Core_0_Isoma_mem.write(reinterpret_cast<char*>(_array_Core_0_Isoma_mem), 256*sizeof(_array_Core_0_Isoma_mem[0]));
		outfile__array_Core_0_Isoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_mem." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_pfb_gain;
	outfile__array_Core_0_Isoma_pfb_gain.open("results/_array_Core_0_Isoma_pfb_gain_-8008490571635100289", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_pfb_gain.is_open())
	{
		outfile__array_Core_0_Isoma_pfb_gain.write(reinterpret_cast<char*>(_array_Core_0_Isoma_pfb_gain), 1*sizeof(_array_Core_0_Isoma_pfb_gain[0]));
		outfile__array_Core_0_Isoma_pfb_gain.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_pfb_gain." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_pfb_norm;
	outfile__array_Core_0_Isoma_pfb_norm.open("results/_array_Core_0_Isoma_pfb_norm_-5154477065552101188", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_pfb_norm.is_open())
	{
		outfile__array_Core_0_Isoma_pfb_norm.write(reinterpret_cast<char*>(_array_Core_0_Isoma_pfb_norm), 1*sizeof(_array_Core_0_Isoma_pfb_norm[0]));
		outfile__array_Core_0_Isoma_pfb_norm.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_pfb_norm." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_pfb_th;
	outfile__array_Core_0_Isoma_pfb_th.open("results/_array_Core_0_Isoma_pfb_th_8158106460523666065", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_pfb_th.is_open())
	{
		outfile__array_Core_0_Isoma_pfb_th.write(reinterpret_cast<char*>(_array_Core_0_Isoma_pfb_th), 1*sizeof(_array_Core_0_Isoma_pfb_th[0]));
		outfile__array_Core_0_Isoma_pfb_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_pfb_th." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_reset;
	outfile__array_Core_0_Isoma_reset.open("results/_array_Core_0_Isoma_reset_1168648158986230994", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_reset.is_open())
	{
		outfile__array_Core_0_Isoma_reset.write(reinterpret_cast<char*>(_array_Core_0_Isoma_reset), 1*sizeof(_array_Core_0_Isoma_reset[0]));
		outfile__array_Core_0_Isoma_reset.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_reset." << endl;
	}
	ofstream outfile__array_Core_0_Isoma_th;
	outfile__array_Core_0_Isoma_th.open("results/_array_Core_0_Isoma_th_-7166173712634343664", ios::binary | ios::out);
	if(outfile__array_Core_0_Isoma_th.is_open())
	{
		outfile__array_Core_0_Isoma_th.write(reinterpret_cast<char*>(_array_Core_0_Isoma_th), 1*sizeof(_array_Core_0_Isoma_th[0]));
		outfile__array_Core_0_Isoma_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Isoma_th." << endl;
	}
	ofstream outfile__array_Core_0_kn;
	outfile__array_Core_0_kn.open("results/_array_Core_0_kn_-2585032584961745283", ios::binary | ios::out);
	if(outfile__array_Core_0_kn.is_open())
	{
		outfile__array_Core_0_kn.write(reinterpret_cast<char*>(_array_Core_0_kn), 1*sizeof(_array_Core_0_kn[0]));
		outfile__array_Core_0_kn.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_kn." << endl;
	}
	ofstream outfile__array_Core_0_kp;
	outfile__array_Core_0_kp.open("results/_array_Core_0_kp_6534127534784763213", ios::binary | ios::out);
	if(outfile__array_Core_0_kp.is_open())
	{
		outfile__array_Core_0_kp.write(reinterpret_cast<char*>(_array_Core_0_kp), 1*sizeof(_array_Core_0_kp[0]));
		outfile__array_Core_0_kp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_kp." << endl;
	}
	ofstream outfile__array_Core_0_lastspike;
	outfile__array_Core_0_lastspike.open("results/_array_Core_0_lastspike_280399881318240997", ios::binary | ios::out);
	if(outfile__array_Core_0_lastspike.is_open())
	{
		outfile__array_Core_0_lastspike.write(reinterpret_cast<char*>(_array_Core_0_lastspike), 256*sizeof(_array_Core_0_lastspike[0]));
		outfile__array_Core_0_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_lastspike." << endl;
	}
	ofstream outfile__array_Core_0_not_refractory;
	outfile__array_Core_0_not_refractory.open("results/_array_Core_0_not_refractory_7608920729004545751", ios::binary | ios::out);
	if(outfile__array_Core_0_not_refractory.is_open())
	{
		outfile__array_Core_0_not_refractory.write(reinterpret_cast<char*>(_array_Core_0_not_refractory), 256*sizeof(_array_Core_0_not_refractory[0]));
		outfile__array_Core_0_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_not_refractory." << endl;
	}
	ofstream outfile__array_Core_0_soma_refP;
	outfile__array_Core_0_soma_refP.open("results/_array_Core_0_soma_refP_-1057281137238136185", ios::binary | ios::out);
	if(outfile__array_Core_0_soma_refP.is_open())
	{
		outfile__array_Core_0_soma_refP.write(reinterpret_cast<char*>(_array_Core_0_soma_refP), 1*sizeof(_array_Core_0_soma_refP[0]));
		outfile__array_Core_0_soma_refP.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_soma_refP." << endl;
	}
	ofstream outfile__array_Core_0_Ut;
	outfile__array_Core_0_Ut.open("results/_array_Core_0_Ut_8332919010992360385", ios::binary | ios::out);
	if(outfile__array_Core_0_Ut.is_open())
	{
		outfile__array_Core_0_Ut.write(reinterpret_cast<char*>(_array_Core_0_Ut), 1*sizeof(_array_Core_0_Ut[0]));
		outfile__array_Core_0_Ut.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Ut." << endl;
	}
	ofstream outfile__array_Core_0_Vnmda;
	outfile__array_Core_0_Vnmda.open("results/_array_Core_0_Vnmda_-1660270172492938575", ios::binary | ios::out);
	if(outfile__array_Core_0_Vnmda.is_open())
	{
		outfile__array_Core_0_Vnmda.write(reinterpret_cast<char*>(_array_Core_0_Vnmda), 256*sizeof(_array_Core_0_Vnmda[0]));
		outfile__array_Core_0_Vnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_0_Vnmda." << endl;
	}
	ofstream outfile__array_Core_1__spikespace;
	outfile__array_Core_1__spikespace.open("results/_array_Core_1__spikespace_4745029304618448235", ios::binary | ios::out);
	if(outfile__array_Core_1__spikespace.is_open())
	{
		outfile__array_Core_1__spikespace.write(reinterpret_cast<char*>(_array_Core_1__spikespace), 257*sizeof(_array_Core_1__spikespace[0]));
		outfile__array_Core_1__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1__spikespace." << endl;
	}
	ofstream outfile__array_Core_1_alpha;
	outfile__array_Core_1_alpha.open("results/_array_Core_1_alpha_571658138794351274", ios::binary | ios::out);
	if(outfile__array_Core_1_alpha.is_open())
	{
		outfile__array_Core_1_alpha.write(reinterpret_cast<char*>(_array_Core_1_alpha), 1*sizeof(_array_Core_1_alpha[0]));
		outfile__array_Core_1_alpha.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_alpha." << endl;
	}
	ofstream outfile__array_Core_1_Campa;
	outfile__array_Core_1_Campa.open("results/_array_Core_1_Campa_8609382552111055249", ios::binary | ios::out);
	if(outfile__array_Core_1_Campa.is_open())
	{
		outfile__array_Core_1_Campa.write(reinterpret_cast<char*>(_array_Core_1_Campa), 256*sizeof(_array_Core_1_Campa[0]));
		outfile__array_Core_1_Campa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Campa." << endl;
	}
	ofstream outfile__array_Core_1_Cgaba_a;
	outfile__array_Core_1_Cgaba_a.open("results/_array_Core_1_Cgaba_a_-8995723779350095530", ios::binary | ios::out);
	if(outfile__array_Core_1_Cgaba_a.is_open())
	{
		outfile__array_Core_1_Cgaba_a.write(reinterpret_cast<char*>(_array_Core_1_Cgaba_a), 256*sizeof(_array_Core_1_Cgaba_a[0]));
		outfile__array_Core_1_Cgaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Cgaba_a." << endl;
	}
	ofstream outfile__array_Core_1_Cgaba_b;
	outfile__array_Core_1_Cgaba_b.open("results/_array_Core_1_Cgaba_b_3663097468023035989", ios::binary | ios::out);
	if(outfile__array_Core_1_Cgaba_b.is_open())
	{
		outfile__array_Core_1_Cgaba_b.write(reinterpret_cast<char*>(_array_Core_1_Cgaba_b), 256*sizeof(_array_Core_1_Cgaba_b[0]));
		outfile__array_Core_1_Cgaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Cgaba_b." << endl;
	}
	ofstream outfile__array_Core_1_Cnmda;
	outfile__array_Core_1_Cnmda.open("results/_array_Core_1_Cnmda_-681375116183163330", ios::binary | ios::out);
	if(outfile__array_Core_1_Cnmda.is_open())
	{
		outfile__array_Core_1_Cnmda.write(reinterpret_cast<char*>(_array_Core_1_Cnmda), 256*sizeof(_array_Core_1_Cnmda[0]));
		outfile__array_Core_1_Cnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Cnmda." << endl;
	}
	ofstream outfile__array_Core_1_Csoma_ahp;
	outfile__array_Core_1_Csoma_ahp.open("results/_array_Core_1_Csoma_ahp_-2699986503282591293", ios::binary | ios::out);
	if(outfile__array_Core_1_Csoma_ahp.is_open())
	{
		outfile__array_Core_1_Csoma_ahp.write(reinterpret_cast<char*>(_array_Core_1_Csoma_ahp), 1*sizeof(_array_Core_1_Csoma_ahp[0]));
		outfile__array_Core_1_Csoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Csoma_ahp." << endl;
	}
	ofstream outfile__array_Core_1_Csoma_mem;
	outfile__array_Core_1_Csoma_mem.open("results/_array_Core_1_Csoma_mem_648782377989776482", ios::binary | ios::out);
	if(outfile__array_Core_1_Csoma_mem.is_open())
	{
		outfile__array_Core_1_Csoma_mem.write(reinterpret_cast<char*>(_array_Core_1_Csoma_mem), 1*sizeof(_array_Core_1_Csoma_mem[0]));
		outfile__array_Core_1_Csoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Csoma_mem." << endl;
	}
	ofstream outfile__array_Core_1_i;
	outfile__array_Core_1_i.open("results/_array_Core_1_i_-3467675460354046822", ios::binary | ios::out);
	if(outfile__array_Core_1_i.is_open())
	{
		outfile__array_Core_1_i.write(reinterpret_cast<char*>(_array_Core_1_i), 256*sizeof(_array_Core_1_i[0]));
		outfile__array_Core_1_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_i." << endl;
	}
	ofstream outfile__array_Core_1_I0;
	outfile__array_Core_1_I0.open("results/_array_Core_1_I0_-6096087739673173138", ios::binary | ios::out);
	if(outfile__array_Core_1_I0.is_open())
	{
		outfile__array_Core_1_I0.write(reinterpret_cast<char*>(_array_Core_1_I0), 1*sizeof(_array_Core_1_I0[0]));
		outfile__array_Core_1_I0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_I0." << endl;
	}
	ofstream outfile__array_Core_1_Iampa;
	outfile__array_Core_1_Iampa.open("results/_array_Core_1_Iampa_-3193502329945830275", ios::binary | ios::out);
	if(outfile__array_Core_1_Iampa.is_open())
	{
		outfile__array_Core_1_Iampa.write(reinterpret_cast<char*>(_array_Core_1_Iampa), 256*sizeof(_array_Core_1_Iampa[0]));
		outfile__array_Core_1_Iampa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Iampa." << endl;
	}
	ofstream outfile__array_Core_1_Iampa_tau;
	outfile__array_Core_1_Iampa_tau.open("results/_array_Core_1_Iampa_tau_7016003092968859404", ios::binary | ios::out);
	if(outfile__array_Core_1_Iampa_tau.is_open())
	{
		outfile__array_Core_1_Iampa_tau.write(reinterpret_cast<char*>(_array_Core_1_Iampa_tau), 256*sizeof(_array_Core_1_Iampa_tau[0]));
		outfile__array_Core_1_Iampa_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Iampa_tau." << endl;
	}
	ofstream outfile__array_Core_1_Iampa_w0;
	outfile__array_Core_1_Iampa_w0.open("results/_array_Core_1_Iampa_w0_-4840741910656773315", ios::binary | ios::out);
	if(outfile__array_Core_1_Iampa_w0.is_open())
	{
		outfile__array_Core_1_Iampa_w0.write(reinterpret_cast<char*>(_array_Core_1_Iampa_w0), 256*sizeof(_array_Core_1_Iampa_w0[0]));
		outfile__array_Core_1_Iampa_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Iampa_w0." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_a;
	outfile__array_Core_1_Igaba_a.open("results/_array_Core_1_Igaba_a_6484421537062458286", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_a.is_open())
	{
		outfile__array_Core_1_Igaba_a.write(reinterpret_cast<char*>(_array_Core_1_Igaba_a), 256*sizeof(_array_Core_1_Igaba_a[0]));
		outfile__array_Core_1_Igaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_a." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_a_tau;
	outfile__array_Core_1_Igaba_a_tau.open("results/_array_Core_1_Igaba_a_tau_-527859343927284139", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_a_tau.is_open())
	{
		outfile__array_Core_1_Igaba_a_tau.write(reinterpret_cast<char*>(_array_Core_1_Igaba_a_tau), 256*sizeof(_array_Core_1_Igaba_a_tau[0]));
		outfile__array_Core_1_Igaba_a_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_a_tau." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_a_w0;
	outfile__array_Core_1_Igaba_a_w0.open("results/_array_Core_1_Igaba_a_w0_3111911350743917276", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_a_w0.is_open())
	{
		outfile__array_Core_1_Igaba_a_w0.write(reinterpret_cast<char*>(_array_Core_1_Igaba_a_w0), 256*sizeof(_array_Core_1_Igaba_a_w0[0]));
		outfile__array_Core_1_Igaba_a_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_a_w0." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_b;
	outfile__array_Core_1_Igaba_b.open("results/_array_Core_1_Igaba_b_4181302973271716044", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_b.is_open())
	{
		outfile__array_Core_1_Igaba_b.write(reinterpret_cast<char*>(_array_Core_1_Igaba_b), 256*sizeof(_array_Core_1_Igaba_b[0]));
		outfile__array_Core_1_Igaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_b." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_b_tau;
	outfile__array_Core_1_Igaba_b_tau.open("results/_array_Core_1_Igaba_b_tau_1255687534998843380", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_b_tau.is_open())
	{
		outfile__array_Core_1_Igaba_b_tau.write(reinterpret_cast<char*>(_array_Core_1_Igaba_b_tau), 256*sizeof(_array_Core_1_Igaba_b_tau[0]));
		outfile__array_Core_1_Igaba_b_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_b_tau." << endl;
	}
	ofstream outfile__array_Core_1_Igaba_b_w0;
	outfile__array_Core_1_Igaba_b_w0.open("results/_array_Core_1_Igaba_b_w0_234580988452980457", ios::binary | ios::out);
	if(outfile__array_Core_1_Igaba_b_w0.is_open())
	{
		outfile__array_Core_1_Igaba_b_w0.write(reinterpret_cast<char*>(_array_Core_1_Igaba_b_w0), 256*sizeof(_array_Core_1_Igaba_b_w0[0]));
		outfile__array_Core_1_Igaba_b_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Igaba_b_w0." << endl;
	}
	ofstream outfile__array_Core_1_Inmda;
	outfile__array_Core_1_Inmda.open("results/_array_Core_1_Inmda_5165625148638551004", ios::binary | ios::out);
	if(outfile__array_Core_1_Inmda.is_open())
	{
		outfile__array_Core_1_Inmda.write(reinterpret_cast<char*>(_array_Core_1_Inmda), 256*sizeof(_array_Core_1_Inmda[0]));
		outfile__array_Core_1_Inmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Inmda." << endl;
	}
	ofstream outfile__array_Core_1_Inmda_tau;
	outfile__array_Core_1_Inmda_tau.open("results/_array_Core_1_Inmda_tau_7570685690612394729", ios::binary | ios::out);
	if(outfile__array_Core_1_Inmda_tau.is_open())
	{
		outfile__array_Core_1_Inmda_tau.write(reinterpret_cast<char*>(_array_Core_1_Inmda_tau), 256*sizeof(_array_Core_1_Inmda_tau[0]));
		outfile__array_Core_1_Inmda_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Inmda_tau." << endl;
	}
	ofstream outfile__array_Core_1_Inmda_w0;
	outfile__array_Core_1_Inmda_w0.open("results/_array_Core_1_Inmda_w0_-1267867838482046779", ios::binary | ios::out);
	if(outfile__array_Core_1_Inmda_w0.is_open())
	{
		outfile__array_Core_1_Inmda_w0.write(reinterpret_cast<char*>(_array_Core_1_Inmda_w0), 256*sizeof(_array_Core_1_Inmda_w0[0]));
		outfile__array_Core_1_Inmda_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Inmda_w0." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_ahp;
	outfile__array_Core_1_Isoma_ahp.open("results/_array_Core_1_Isoma_ahp_6912298048512682808", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_ahp.is_open())
	{
		outfile__array_Core_1_Isoma_ahp.write(reinterpret_cast<char*>(_array_Core_1_Isoma_ahp), 256*sizeof(_array_Core_1_Isoma_ahp[0]));
		outfile__array_Core_1_Isoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_ahp." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_ahp_tau;
	outfile__array_Core_1_Isoma_ahp_tau.open("results/_array_Core_1_Isoma_ahp_tau_7223716057247902436", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_ahp_tau.is_open())
	{
		outfile__array_Core_1_Isoma_ahp_tau.write(reinterpret_cast<char*>(_array_Core_1_Isoma_ahp_tau), 1*sizeof(_array_Core_1_Isoma_ahp_tau[0]));
		outfile__array_Core_1_Isoma_ahp_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_ahp_tau." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_ahp_th;
	outfile__array_Core_1_Isoma_ahp_th.open("results/_array_Core_1_Isoma_ahp_th_7404222482541630090", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_ahp_th.is_open())
	{
		outfile__array_Core_1_Isoma_ahp_th.write(reinterpret_cast<char*>(_array_Core_1_Isoma_ahp_th), 1*sizeof(_array_Core_1_Isoma_ahp_th[0]));
		outfile__array_Core_1_Isoma_ahp_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_ahp_th." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_ahp_w;
	outfile__array_Core_1_Isoma_ahp_w.open("results/_array_Core_1_Isoma_ahp_w_7072443153751518999", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_ahp_w.is_open())
	{
		outfile__array_Core_1_Isoma_ahp_w.write(reinterpret_cast<char*>(_array_Core_1_Isoma_ahp_w), 256*sizeof(_array_Core_1_Isoma_ahp_w[0]));
		outfile__array_Core_1_Isoma_ahp_w.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_ahp_w." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_const;
	outfile__array_Core_1_Isoma_const.open("results/_array_Core_1_Isoma_const_-3741468716725083720", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_const.is_open())
	{
		outfile__array_Core_1_Isoma_const.write(reinterpret_cast<char*>(_array_Core_1_Isoma_const), 256*sizeof(_array_Core_1_Isoma_const[0]));
		outfile__array_Core_1_Isoma_const.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_const." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_dpi_tau;
	outfile__array_Core_1_Isoma_dpi_tau.open("results/_array_Core_1_Isoma_dpi_tau_-3795651691745873183", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_dpi_tau.is_open())
	{
		outfile__array_Core_1_Isoma_dpi_tau.write(reinterpret_cast<char*>(_array_Core_1_Isoma_dpi_tau), 1*sizeof(_array_Core_1_Isoma_dpi_tau[0]));
		outfile__array_Core_1_Isoma_dpi_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_dpi_tau." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_mem;
	outfile__array_Core_1_Isoma_mem.open("results/_array_Core_1_Isoma_mem_6763980833865717731", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_mem.is_open())
	{
		outfile__array_Core_1_Isoma_mem.write(reinterpret_cast<char*>(_array_Core_1_Isoma_mem), 256*sizeof(_array_Core_1_Isoma_mem[0]));
		outfile__array_Core_1_Isoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_mem." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_pfb_gain;
	outfile__array_Core_1_Isoma_pfb_gain.open("results/_array_Core_1_Isoma_pfb_gain_-6323092088045192551", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_pfb_gain.is_open())
	{
		outfile__array_Core_1_Isoma_pfb_gain.write(reinterpret_cast<char*>(_array_Core_1_Isoma_pfb_gain), 1*sizeof(_array_Core_1_Isoma_pfb_gain[0]));
		outfile__array_Core_1_Isoma_pfb_gain.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_pfb_gain." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_pfb_norm;
	outfile__array_Core_1_Isoma_pfb_norm.open("results/_array_Core_1_Isoma_pfb_norm_4943974822919639169", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_pfb_norm.is_open())
	{
		outfile__array_Core_1_Isoma_pfb_norm.write(reinterpret_cast<char*>(_array_Core_1_Isoma_pfb_norm), 1*sizeof(_array_Core_1_Isoma_pfb_norm[0]));
		outfile__array_Core_1_Isoma_pfb_norm.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_pfb_norm." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_pfb_th;
	outfile__array_Core_1_Isoma_pfb_th.open("results/_array_Core_1_Isoma_pfb_th_-1940280114659747011", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_pfb_th.is_open())
	{
		outfile__array_Core_1_Isoma_pfb_th.write(reinterpret_cast<char*>(_array_Core_1_Isoma_pfb_th), 1*sizeof(_array_Core_1_Isoma_pfb_th[0]));
		outfile__array_Core_1_Isoma_pfb_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_pfb_th." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_reset;
	outfile__array_Core_1_Isoma_reset.open("results/_array_Core_1_Isoma_reset_-5972044550959237026", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_reset.is_open())
	{
		outfile__array_Core_1_Isoma_reset.write(reinterpret_cast<char*>(_array_Core_1_Isoma_reset), 1*sizeof(_array_Core_1_Isoma_reset[0]));
		outfile__array_Core_1_Isoma_reset.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_reset." << endl;
	}
	ofstream outfile__array_Core_1_Isoma_th;
	outfile__array_Core_1_Isoma_th.open("results/_array_Core_1_Isoma_th_3562425046288098391", ios::binary | ios::out);
	if(outfile__array_Core_1_Isoma_th.is_open())
	{
		outfile__array_Core_1_Isoma_th.write(reinterpret_cast<char*>(_array_Core_1_Isoma_th), 1*sizeof(_array_Core_1_Isoma_th[0]));
		outfile__array_Core_1_Isoma_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Isoma_th." << endl;
	}
	ofstream outfile__array_Core_1_kn;
	outfile__array_Core_1_kn.open("results/_array_Core_1_kn_-7632872205077892285", ios::binary | ios::out);
	if(outfile__array_Core_1_kn.is_open())
	{
		outfile__array_Core_1_kn.write(reinterpret_cast<char*>(_array_Core_1_kn), 1*sizeof(_array_Core_1_kn[0]));
		outfile__array_Core_1_kn.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_kn." << endl;
	}
	ofstream outfile__array_Core_1_kp;
	outfile__array_Core_1_kp.open("results/_array_Core_1_kp_1377722054779167041", ios::binary | ios::out);
	if(outfile__array_Core_1_kp.is_open())
	{
		outfile__array_Core_1_kp.write(reinterpret_cast<char*>(_array_Core_1_kp), 1*sizeof(_array_Core_1_kp[0]));
		outfile__array_Core_1_kp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_kp." << endl;
	}
	ofstream outfile__array_Core_1_lastspike;
	outfile__array_Core_1_lastspike.open("results/_array_Core_1_lastspike_-364132331603415851", ios::binary | ios::out);
	if(outfile__array_Core_1_lastspike.is_open())
	{
		outfile__array_Core_1_lastspike.write(reinterpret_cast<char*>(_array_Core_1_lastspike), 256*sizeof(_array_Core_1_lastspike[0]));
		outfile__array_Core_1_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_lastspike." << endl;
	}
	ofstream outfile__array_Core_1_not_refractory;
	outfile__array_Core_1_not_refractory.open("results/_array_Core_1_not_refractory_-8385826990595957215", ios::binary | ios::out);
	if(outfile__array_Core_1_not_refractory.is_open())
	{
		outfile__array_Core_1_not_refractory.write(reinterpret_cast<char*>(_array_Core_1_not_refractory), 256*sizeof(_array_Core_1_not_refractory[0]));
		outfile__array_Core_1_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_not_refractory." << endl;
	}
	ofstream outfile__array_Core_1_soma_refP;
	outfile__array_Core_1_soma_refP.open("results/_array_Core_1_soma_refP_9064680826640654083", ios::binary | ios::out);
	if(outfile__array_Core_1_soma_refP.is_open())
	{
		outfile__array_Core_1_soma_refP.write(reinterpret_cast<char*>(_array_Core_1_soma_refP), 1*sizeof(_array_Core_1_soma_refP[0]));
		outfile__array_Core_1_soma_refP.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_soma_refP." << endl;
	}
	ofstream outfile__array_Core_1_subgroup__sub_idx;
	outfile__array_Core_1_subgroup__sub_idx.open("results/_array_Core_1_subgroup__sub_idx_3899378301309207361", ios::binary | ios::out);
	if(outfile__array_Core_1_subgroup__sub_idx.is_open())
	{
		outfile__array_Core_1_subgroup__sub_idx.write(reinterpret_cast<char*>(_array_Core_1_subgroup__sub_idx), 200*sizeof(_array_Core_1_subgroup__sub_idx[0]));
		outfile__array_Core_1_subgroup__sub_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_subgroup__sub_idx." << endl;
	}
	ofstream outfile__array_Core_1_Ut;
	outfile__array_Core_1_Ut.open("results/_array_Core_1_Ut_8561731026037120197", ios::binary | ios::out);
	if(outfile__array_Core_1_Ut.is_open())
	{
		outfile__array_Core_1_Ut.write(reinterpret_cast<char*>(_array_Core_1_Ut), 1*sizeof(_array_Core_1_Ut[0]));
		outfile__array_Core_1_Ut.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Ut." << endl;
	}
	ofstream outfile__array_Core_1_Vnmda;
	outfile__array_Core_1_Vnmda.open("results/_array_Core_1_Vnmda_3067544200709820555", ios::binary | ios::out);
	if(outfile__array_Core_1_Vnmda.is_open())
	{
		outfile__array_Core_1_Vnmda.write(reinterpret_cast<char*>(_array_Core_1_Vnmda), 256*sizeof(_array_Core_1_Vnmda[0]));
		outfile__array_Core_1_Vnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_1_Vnmda." << endl;
	}
	ofstream outfile__array_Core_2__spikespace;
	outfile__array_Core_2__spikespace.open("results/_array_Core_2__spikespace_2714288011452248463", ios::binary | ios::out);
	if(outfile__array_Core_2__spikespace.is_open())
	{
		outfile__array_Core_2__spikespace.write(reinterpret_cast<char*>(_array_Core_2__spikespace), 257*sizeof(_array_Core_2__spikespace[0]));
		outfile__array_Core_2__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2__spikespace." << endl;
	}
	ofstream outfile__array_Core_2_alpha;
	outfile__array_Core_2_alpha.open("results/_array_Core_2_alpha_4829133918550077155", ios::binary | ios::out);
	if(outfile__array_Core_2_alpha.is_open())
	{
		outfile__array_Core_2_alpha.write(reinterpret_cast<char*>(_array_Core_2_alpha), 1*sizeof(_array_Core_2_alpha[0]));
		outfile__array_Core_2_alpha.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_alpha." << endl;
	}
	ofstream outfile__array_Core_2_Campa;
	outfile__array_Core_2_Campa.open("results/_array_Core_2_Campa_7047386180542495969", ios::binary | ios::out);
	if(outfile__array_Core_2_Campa.is_open())
	{
		outfile__array_Core_2_Campa.write(reinterpret_cast<char*>(_array_Core_2_Campa), 256*sizeof(_array_Core_2_Campa[0]));
		outfile__array_Core_2_Campa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Campa." << endl;
	}
	ofstream outfile__array_Core_2_Cgaba_a;
	outfile__array_Core_2_Cgaba_a.open("results/_array_Core_2_Cgaba_a_-8905278182133238360", ios::binary | ios::out);
	if(outfile__array_Core_2_Cgaba_a.is_open())
	{
		outfile__array_Core_2_Cgaba_a.write(reinterpret_cast<char*>(_array_Core_2_Cgaba_a), 256*sizeof(_array_Core_2_Cgaba_a[0]));
		outfile__array_Core_2_Cgaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Cgaba_a." << endl;
	}
	ofstream outfile__array_Core_2_Cgaba_b;
	outfile__array_Core_2_Cgaba_b.open("results/_array_Core_2_Cgaba_b_-1971229590902507628", ios::binary | ios::out);
	if(outfile__array_Core_2_Cgaba_b.is_open())
	{
		outfile__array_Core_2_Cgaba_b.write(reinterpret_cast<char*>(_array_Core_2_Cgaba_b), 256*sizeof(_array_Core_2_Cgaba_b[0]));
		outfile__array_Core_2_Cgaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Cgaba_b." << endl;
	}
	ofstream outfile__array_Core_2_Cnmda;
	outfile__array_Core_2_Cnmda.open("results/_array_Core_2_Cnmda_-1986683408901046977", ios::binary | ios::out);
	if(outfile__array_Core_2_Cnmda.is_open())
	{
		outfile__array_Core_2_Cnmda.write(reinterpret_cast<char*>(_array_Core_2_Cnmda), 256*sizeof(_array_Core_2_Cnmda[0]));
		outfile__array_Core_2_Cnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Cnmda." << endl;
	}
	ofstream outfile__array_Core_2_Csoma_ahp;
	outfile__array_Core_2_Csoma_ahp.open("results/_array_Core_2_Csoma_ahp_1138092045845657559", ios::binary | ios::out);
	if(outfile__array_Core_2_Csoma_ahp.is_open())
	{
		outfile__array_Core_2_Csoma_ahp.write(reinterpret_cast<char*>(_array_Core_2_Csoma_ahp), 1*sizeof(_array_Core_2_Csoma_ahp[0]));
		outfile__array_Core_2_Csoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Csoma_ahp." << endl;
	}
	ofstream outfile__array_Core_2_Csoma_mem;
	outfile__array_Core_2_Csoma_mem.open("results/_array_Core_2_Csoma_mem_-6925124155335664906", ios::binary | ios::out);
	if(outfile__array_Core_2_Csoma_mem.is_open())
	{
		outfile__array_Core_2_Csoma_mem.write(reinterpret_cast<char*>(_array_Core_2_Csoma_mem), 1*sizeof(_array_Core_2_Csoma_mem[0]));
		outfile__array_Core_2_Csoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Csoma_mem." << endl;
	}
	ofstream outfile__array_Core_2_i;
	outfile__array_Core_2_i.open("results/_array_Core_2_i_2968909171995655679", ios::binary | ios::out);
	if(outfile__array_Core_2_i.is_open())
	{
		outfile__array_Core_2_i.write(reinterpret_cast<char*>(_array_Core_2_i), 256*sizeof(_array_Core_2_i[0]));
		outfile__array_Core_2_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_i." << endl;
	}
	ofstream outfile__array_Core_2_I0;
	outfile__array_Core_2_I0.open("results/_array_Core_2_I0_682878840923591516", ios::binary | ios::out);
	if(outfile__array_Core_2_I0.is_open())
	{
		outfile__array_Core_2_I0.write(reinterpret_cast<char*>(_array_Core_2_I0), 1*sizeof(_array_Core_2_I0[0]));
		outfile__array_Core_2_I0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_I0." << endl;
	}
	ofstream outfile__array_Core_2_Iampa;
	outfile__array_Core_2_Iampa.open("results/_array_Core_2_Iampa_-974519474301591764", ios::binary | ios::out);
	if(outfile__array_Core_2_Iampa.is_open())
	{
		outfile__array_Core_2_Iampa.write(reinterpret_cast<char*>(_array_Core_2_Iampa), 256*sizeof(_array_Core_2_Iampa[0]));
		outfile__array_Core_2_Iampa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Iampa." << endl;
	}
	ofstream outfile__array_Core_2_Iampa_tau;
	outfile__array_Core_2_Iampa_tau.open("results/_array_Core_2_Iampa_tau_-544179754447697738", ios::binary | ios::out);
	if(outfile__array_Core_2_Iampa_tau.is_open())
	{
		outfile__array_Core_2_Iampa_tau.write(reinterpret_cast<char*>(_array_Core_2_Iampa_tau), 256*sizeof(_array_Core_2_Iampa_tau[0]));
		outfile__array_Core_2_Iampa_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Iampa_tau." << endl;
	}
	ofstream outfile__array_Core_2_Iampa_w0;
	outfile__array_Core_2_Iampa_w0.open("results/_array_Core_2_Iampa_w0_-4669962256799230612", ios::binary | ios::out);
	if(outfile__array_Core_2_Iampa_w0.is_open())
	{
		outfile__array_Core_2_Iampa_w0.write(reinterpret_cast<char*>(_array_Core_2_Iampa_w0), 256*sizeof(_array_Core_2_Iampa_w0[0]));
		outfile__array_Core_2_Iampa_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Iampa_w0." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_a;
	outfile__array_Core_2_Igaba_a.open("results/_array_Core_2_Igaba_a_-8178172425224408915", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_a.is_open())
	{
		outfile__array_Core_2_Igaba_a.write(reinterpret_cast<char*>(_array_Core_2_Igaba_a), 256*sizeof(_array_Core_2_Igaba_a[0]));
		outfile__array_Core_2_Igaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_a." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_a_tau;
	outfile__array_Core_2_Igaba_a_tau.open("results/_array_Core_2_Igaba_a_tau_-1929011691668370533", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_a_tau.is_open())
	{
		outfile__array_Core_2_Igaba_a_tau.write(reinterpret_cast<char*>(_array_Core_2_Igaba_a_tau), 256*sizeof(_array_Core_2_Igaba_a_tau[0]));
		outfile__array_Core_2_Igaba_a_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_a_tau." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_a_w0;
	outfile__array_Core_2_Igaba_a_w0.open("results/_array_Core_2_Igaba_a_w0_-208738957315584050", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_a_w0.is_open())
	{
		outfile__array_Core_2_Igaba_a_w0.write(reinterpret_cast<char*>(_array_Core_2_Igaba_a_w0), 256*sizeof(_array_Core_2_Igaba_a_w0[0]));
		outfile__array_Core_2_Igaba_a_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_a_w0." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_b;
	outfile__array_Core_2_Igaba_b.open("results/_array_Core_2_Igaba_b_8301307780769450769", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_b.is_open())
	{
		outfile__array_Core_2_Igaba_b.write(reinterpret_cast<char*>(_array_Core_2_Igaba_b), 256*sizeof(_array_Core_2_Igaba_b[0]));
		outfile__array_Core_2_Igaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_b." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_b_tau;
	outfile__array_Core_2_Igaba_b_tau.open("results/_array_Core_2_Igaba_b_tau_1374943908093485507", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_b_tau.is_open())
	{
		outfile__array_Core_2_Igaba_b_tau.write(reinterpret_cast<char*>(_array_Core_2_Igaba_b_tau), 256*sizeof(_array_Core_2_Igaba_b_tau[0]));
		outfile__array_Core_2_Igaba_b_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_b_tau." << endl;
	}
	ofstream outfile__array_Core_2_Igaba_b_w0;
	outfile__array_Core_2_Igaba_b_w0.open("results/_array_Core_2_Igaba_b_w0_-5816505228291761172", ios::binary | ios::out);
	if(outfile__array_Core_2_Igaba_b_w0.is_open())
	{
		outfile__array_Core_2_Igaba_b_w0.write(reinterpret_cast<char*>(_array_Core_2_Igaba_b_w0), 256*sizeof(_array_Core_2_Igaba_b_w0[0]));
		outfile__array_Core_2_Igaba_b_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Igaba_b_w0." << endl;
	}
	ofstream outfile__array_Core_2_Inmda;
	outfile__array_Core_2_Inmda.open("results/_array_Core_2_Inmda_7539677657246542543", ios::binary | ios::out);
	if(outfile__array_Core_2_Inmda.is_open())
	{
		outfile__array_Core_2_Inmda.write(reinterpret_cast<char*>(_array_Core_2_Inmda), 256*sizeof(_array_Core_2_Inmda[0]));
		outfile__array_Core_2_Inmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Inmda." << endl;
	}
	ofstream outfile__array_Core_2_Inmda_tau;
	outfile__array_Core_2_Inmda_tau.open("results/_array_Core_2_Inmda_tau_5841666628376730255", ios::binary | ios::out);
	if(outfile__array_Core_2_Inmda_tau.is_open())
	{
		outfile__array_Core_2_Inmda_tau.write(reinterpret_cast<char*>(_array_Core_2_Inmda_tau), 256*sizeof(_array_Core_2_Inmda_tau[0]));
		outfile__array_Core_2_Inmda_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Inmda_tau." << endl;
	}
	ofstream outfile__array_Core_2_Inmda_w0;
	outfile__array_Core_2_Inmda_w0.open("results/_array_Core_2_Inmda_w0_-3925476476951872238", ios::binary | ios::out);
	if(outfile__array_Core_2_Inmda_w0.is_open())
	{
		outfile__array_Core_2_Inmda_w0.write(reinterpret_cast<char*>(_array_Core_2_Inmda_w0), 256*sizeof(_array_Core_2_Inmda_w0[0]));
		outfile__array_Core_2_Inmda_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Inmda_w0." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_ahp;
	outfile__array_Core_2_Isoma_ahp.open("results/_array_Core_2_Isoma_ahp_5357668929509119229", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_ahp.is_open())
	{
		outfile__array_Core_2_Isoma_ahp.write(reinterpret_cast<char*>(_array_Core_2_Isoma_ahp), 256*sizeof(_array_Core_2_Isoma_ahp[0]));
		outfile__array_Core_2_Isoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_ahp." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_ahp_tau;
	outfile__array_Core_2_Isoma_ahp_tau.open("results/_array_Core_2_Isoma_ahp_tau_-5655943616640550978", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_ahp_tau.is_open())
	{
		outfile__array_Core_2_Isoma_ahp_tau.write(reinterpret_cast<char*>(_array_Core_2_Isoma_ahp_tau), 1*sizeof(_array_Core_2_Isoma_ahp_tau[0]));
		outfile__array_Core_2_Isoma_ahp_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_ahp_tau." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_ahp_th;
	outfile__array_Core_2_Isoma_ahp_th.open("results/_array_Core_2_Isoma_ahp_th_1323424557956543554", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_ahp_th.is_open())
	{
		outfile__array_Core_2_Isoma_ahp_th.write(reinterpret_cast<char*>(_array_Core_2_Isoma_ahp_th), 1*sizeof(_array_Core_2_Isoma_ahp_th[0]));
		outfile__array_Core_2_Isoma_ahp_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_ahp_th." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_ahp_w;
	outfile__array_Core_2_Isoma_ahp_w.open("results/_array_Core_2_Isoma_ahp_w_-684021785476012181", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_ahp_w.is_open())
	{
		outfile__array_Core_2_Isoma_ahp_w.write(reinterpret_cast<char*>(_array_Core_2_Isoma_ahp_w), 256*sizeof(_array_Core_2_Isoma_ahp_w[0]));
		outfile__array_Core_2_Isoma_ahp_w.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_ahp_w." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_const;
	outfile__array_Core_2_Isoma_const.open("results/_array_Core_2_Isoma_const_4039152237381443532", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_const.is_open())
	{
		outfile__array_Core_2_Isoma_const.write(reinterpret_cast<char*>(_array_Core_2_Isoma_const), 256*sizeof(_array_Core_2_Isoma_const[0]));
		outfile__array_Core_2_Isoma_const.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_const." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_dpi_tau;
	outfile__array_Core_2_Isoma_dpi_tau.open("results/_array_Core_2_Isoma_dpi_tau_3046647324501871057", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_dpi_tau.is_open())
	{
		outfile__array_Core_2_Isoma_dpi_tau.write(reinterpret_cast<char*>(_array_Core_2_Isoma_dpi_tau), 1*sizeof(_array_Core_2_Isoma_dpi_tau[0]));
		outfile__array_Core_2_Isoma_dpi_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_dpi_tau." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_mem;
	outfile__array_Core_2_Isoma_mem.open("results/_array_Core_2_Isoma_mem_3157564184355940062", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_mem.is_open())
	{
		outfile__array_Core_2_Isoma_mem.write(reinterpret_cast<char*>(_array_Core_2_Isoma_mem), 256*sizeof(_array_Core_2_Isoma_mem[0]));
		outfile__array_Core_2_Isoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_mem." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_pfb_gain;
	outfile__array_Core_2_Isoma_pfb_gain.open("results/_array_Core_2_Isoma_pfb_gain_-7420569219266410542", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_pfb_gain.is_open())
	{
		outfile__array_Core_2_Isoma_pfb_gain.write(reinterpret_cast<char*>(_array_Core_2_Isoma_pfb_gain), 1*sizeof(_array_Core_2_Isoma_pfb_gain[0]));
		outfile__array_Core_2_Isoma_pfb_gain.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_pfb_gain." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_pfb_norm;
	outfile__array_Core_2_Isoma_pfb_norm.open("results/_array_Core_2_Isoma_pfb_norm_5549272964492540823", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_pfb_norm.is_open())
	{
		outfile__array_Core_2_Isoma_pfb_norm.write(reinterpret_cast<char*>(_array_Core_2_Isoma_pfb_norm), 1*sizeof(_array_Core_2_Isoma_pfb_norm[0]));
		outfile__array_Core_2_Isoma_pfb_norm.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_pfb_norm." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_pfb_th;
	outfile__array_Core_2_Isoma_pfb_th.open("results/_array_Core_2_Isoma_pfb_th_-4891874736898815128", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_pfb_th.is_open())
	{
		outfile__array_Core_2_Isoma_pfb_th.write(reinterpret_cast<char*>(_array_Core_2_Isoma_pfb_th), 1*sizeof(_array_Core_2_Isoma_pfb_th[0]));
		outfile__array_Core_2_Isoma_pfb_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_pfb_th." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_reset;
	outfile__array_Core_2_Isoma_reset.open("results/_array_Core_2_Isoma_reset_9053691545112417306", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_reset.is_open())
	{
		outfile__array_Core_2_Isoma_reset.write(reinterpret_cast<char*>(_array_Core_2_Isoma_reset), 1*sizeof(_array_Core_2_Isoma_reset[0]));
		outfile__array_Core_2_Isoma_reset.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_reset." << endl;
	}
	ofstream outfile__array_Core_2_Isoma_th;
	outfile__array_Core_2_Isoma_th.open("results/_array_Core_2_Isoma_th_4014757945149860417", ios::binary | ios::out);
	if(outfile__array_Core_2_Isoma_th.is_open())
	{
		outfile__array_Core_2_Isoma_th.write(reinterpret_cast<char*>(_array_Core_2_Isoma_th), 1*sizeof(_array_Core_2_Isoma_th[0]));
		outfile__array_Core_2_Isoma_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Isoma_th." << endl;
	}
	ofstream outfile__array_Core_2_kn;
	outfile__array_Core_2_kn.open("results/_array_Core_2_kn_5461670827021580351", ios::binary | ios::out);
	if(outfile__array_Core_2_kn.is_open())
	{
		outfile__array_Core_2_kn.write(reinterpret_cast<char*>(_array_Core_2_kn), 1*sizeof(_array_Core_2_kn[0]));
		outfile__array_Core_2_kn.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_kn." << endl;
	}
	ofstream outfile__array_Core_2_kp;
	outfile__array_Core_2_kp.open("results/_array_Core_2_kp_1030090911550580157", ios::binary | ios::out);
	if(outfile__array_Core_2_kp.is_open())
	{
		outfile__array_Core_2_kp.write(reinterpret_cast<char*>(_array_Core_2_kp), 1*sizeof(_array_Core_2_kp[0]));
		outfile__array_Core_2_kp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_kp." << endl;
	}
	ofstream outfile__array_Core_2_lastspike;
	outfile__array_Core_2_lastspike.open("results/_array_Core_2_lastspike_2829611964094689767", ios::binary | ios::out);
	if(outfile__array_Core_2_lastspike.is_open())
	{
		outfile__array_Core_2_lastspike.write(reinterpret_cast<char*>(_array_Core_2_lastspike), 256*sizeof(_array_Core_2_lastspike[0]));
		outfile__array_Core_2_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_lastspike." << endl;
	}
	ofstream outfile__array_Core_2_not_refractory;
	outfile__array_Core_2_not_refractory.open("results/_array_Core_2_not_refractory_-6684990326119076777", ios::binary | ios::out);
	if(outfile__array_Core_2_not_refractory.is_open())
	{
		outfile__array_Core_2_not_refractory.write(reinterpret_cast<char*>(_array_Core_2_not_refractory), 256*sizeof(_array_Core_2_not_refractory[0]));
		outfile__array_Core_2_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_not_refractory." << endl;
	}
	ofstream outfile__array_Core_2_soma_refP;
	outfile__array_Core_2_soma_refP.open("results/_array_Core_2_soma_refP_1669772621464462177", ios::binary | ios::out);
	if(outfile__array_Core_2_soma_refP.is_open())
	{
		outfile__array_Core_2_soma_refP.write(reinterpret_cast<char*>(_array_Core_2_soma_refP), 1*sizeof(_array_Core_2_soma_refP[0]));
		outfile__array_Core_2_soma_refP.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_soma_refP." << endl;
	}
	ofstream outfile__array_Core_2_Ut;
	outfile__array_Core_2_Ut.open("results/_array_Core_2_Ut_-1343622603439472376", ios::binary | ios::out);
	if(outfile__array_Core_2_Ut.is_open())
	{
		outfile__array_Core_2_Ut.write(reinterpret_cast<char*>(_array_Core_2_Ut), 1*sizeof(_array_Core_2_Ut[0]));
		outfile__array_Core_2_Ut.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Ut." << endl;
	}
	ofstream outfile__array_Core_2_Vnmda;
	outfile__array_Core_2_Vnmda.open("results/_array_Core_2_Vnmda_8312045330919369892", ios::binary | ios::out);
	if(outfile__array_Core_2_Vnmda.is_open())
	{
		outfile__array_Core_2_Vnmda.write(reinterpret_cast<char*>(_array_Core_2_Vnmda), 256*sizeof(_array_Core_2_Vnmda[0]));
		outfile__array_Core_2_Vnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_2_Vnmda." << endl;
	}
	ofstream outfile__array_Core_3__spikespace;
	outfile__array_Core_3__spikespace.open("results/_array_Core_3__spikespace_925090370136711148", ios::binary | ios::out);
	if(outfile__array_Core_3__spikespace.is_open())
	{
		outfile__array_Core_3__spikespace.write(reinterpret_cast<char*>(_array_Core_3__spikespace), 257*sizeof(_array_Core_3__spikespace[0]));
		outfile__array_Core_3__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3__spikespace." << endl;
	}
	ofstream outfile__array_Core_3_alpha;
	outfile__array_Core_3_alpha.open("results/_array_Core_3_alpha_7023912066492965080", ios::binary | ios::out);
	if(outfile__array_Core_3_alpha.is_open())
	{
		outfile__array_Core_3_alpha.write(reinterpret_cast<char*>(_array_Core_3_alpha), 1*sizeof(_array_Core_3_alpha[0]));
		outfile__array_Core_3_alpha.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_alpha." << endl;
	}
	ofstream outfile__array_Core_3_Campa;
	outfile__array_Core_3_Campa.open("results/_array_Core_3_Campa_5609373316574156489", ios::binary | ios::out);
	if(outfile__array_Core_3_Campa.is_open())
	{
		outfile__array_Core_3_Campa.write(reinterpret_cast<char*>(_array_Core_3_Campa), 256*sizeof(_array_Core_3_Campa[0]));
		outfile__array_Core_3_Campa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Campa." << endl;
	}
	ofstream outfile__array_Core_3_Cgaba_a;
	outfile__array_Core_3_Cgaba_a.open("results/_array_Core_3_Cgaba_a_-5882218436572377360", ios::binary | ios::out);
	if(outfile__array_Core_3_Cgaba_a.is_open())
	{
		outfile__array_Core_3_Cgaba_a.write(reinterpret_cast<char*>(_array_Core_3_Cgaba_a), 256*sizeof(_array_Core_3_Cgaba_a[0]));
		outfile__array_Core_3_Cgaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Cgaba_a." << endl;
	}
	ofstream outfile__array_Core_3_Cgaba_b;
	outfile__array_Core_3_Cgaba_b.open("results/_array_Core_3_Cgaba_b_3618821169187363468", ios::binary | ios::out);
	if(outfile__array_Core_3_Cgaba_b.is_open())
	{
		outfile__array_Core_3_Cgaba_b.write(reinterpret_cast<char*>(_array_Core_3_Cgaba_b), 256*sizeof(_array_Core_3_Cgaba_b[0]));
		outfile__array_Core_3_Cgaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Cgaba_b." << endl;
	}
	ofstream outfile__array_Core_3_Cnmda;
	outfile__array_Core_3_Cnmda.open("results/_array_Core_3_Cnmda_8116837111563939204", ios::binary | ios::out);
	if(outfile__array_Core_3_Cnmda.is_open())
	{
		outfile__array_Core_3_Cnmda.write(reinterpret_cast<char*>(_array_Core_3_Cnmda), 256*sizeof(_array_Core_3_Cnmda[0]));
		outfile__array_Core_3_Cnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Cnmda." << endl;
	}
	ofstream outfile__array_Core_3_Csoma_ahp;
	outfile__array_Core_3_Csoma_ahp.open("results/_array_Core_3_Csoma_ahp_-3793483691755332528", ios::binary | ios::out);
	if(outfile__array_Core_3_Csoma_ahp.is_open())
	{
		outfile__array_Core_3_Csoma_ahp.write(reinterpret_cast<char*>(_array_Core_3_Csoma_ahp), 1*sizeof(_array_Core_3_Csoma_ahp[0]));
		outfile__array_Core_3_Csoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Csoma_ahp." << endl;
	}
	ofstream outfile__array_Core_3_Csoma_mem;
	outfile__array_Core_3_Csoma_mem.open("results/_array_Core_3_Csoma_mem_-5889416037957276374", ios::binary | ios::out);
	if(outfile__array_Core_3_Csoma_mem.is_open())
	{
		outfile__array_Core_3_Csoma_mem.write(reinterpret_cast<char*>(_array_Core_3_Csoma_mem), 1*sizeof(_array_Core_3_Csoma_mem[0]));
		outfile__array_Core_3_Csoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Csoma_mem." << endl;
	}
	ofstream outfile__array_Core_3_i;
	outfile__array_Core_3_i.open("results/_array_Core_3_i_8591714197859633164", ios::binary | ios::out);
	if(outfile__array_Core_3_i.is_open())
	{
		outfile__array_Core_3_i.write(reinterpret_cast<char*>(_array_Core_3_i), 256*sizeof(_array_Core_3_i[0]));
		outfile__array_Core_3_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_i." << endl;
	}
	ofstream outfile__array_Core_3_I0;
	outfile__array_Core_3_I0.open("results/_array_Core_3_I0_3075125978108322843", ios::binary | ios::out);
	if(outfile__array_Core_3_I0.is_open())
	{
		outfile__array_Core_3_I0.write(reinterpret_cast<char*>(_array_Core_3_I0), 1*sizeof(_array_Core_3_I0[0]));
		outfile__array_Core_3_I0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_I0." << endl;
	}
	ofstream outfile__array_Core_3_Iampa;
	outfile__array_Core_3_Iampa.open("results/_array_Core_3_Iampa_19337635397207523", ios::binary | ios::out);
	if(outfile__array_Core_3_Iampa.is_open())
	{
		outfile__array_Core_3_Iampa.write(reinterpret_cast<char*>(_array_Core_3_Iampa), 256*sizeof(_array_Core_3_Iampa[0]));
		outfile__array_Core_3_Iampa.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Iampa." << endl;
	}
	ofstream outfile__array_Core_3_Iampa_tau;
	outfile__array_Core_3_Iampa_tau.open("results/_array_Core_3_Iampa_tau_-5291959108577628328", ios::binary | ios::out);
	if(outfile__array_Core_3_Iampa_tau.is_open())
	{
		outfile__array_Core_3_Iampa_tau.write(reinterpret_cast<char*>(_array_Core_3_Iampa_tau), 256*sizeof(_array_Core_3_Iampa_tau[0]));
		outfile__array_Core_3_Iampa_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Iampa_tau." << endl;
	}
	ofstream outfile__array_Core_3_Iampa_w0;
	outfile__array_Core_3_Iampa_w0.open("results/_array_Core_3_Iampa_w0_1106509828622934673", ios::binary | ios::out);
	if(outfile__array_Core_3_Iampa_w0.is_open())
	{
		outfile__array_Core_3_Iampa_w0.write(reinterpret_cast<char*>(_array_Core_3_Iampa_w0), 256*sizeof(_array_Core_3_Iampa_w0[0]));
		outfile__array_Core_3_Iampa_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Iampa_w0." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_a;
	outfile__array_Core_3_Igaba_a.open("results/_array_Core_3_Igaba_a_-2449380911012052708", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_a.is_open())
	{
		outfile__array_Core_3_Igaba_a.write(reinterpret_cast<char*>(_array_Core_3_Igaba_a), 256*sizeof(_array_Core_3_Igaba_a[0]));
		outfile__array_Core_3_Igaba_a.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_a." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_a_tau;
	outfile__array_Core_3_Igaba_a_tau.open("results/_array_Core_3_Igaba_a_tau_4553950090293187123", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_a_tau.is_open())
	{
		outfile__array_Core_3_Igaba_a_tau.write(reinterpret_cast<char*>(_array_Core_3_Igaba_a_tau), 256*sizeof(_array_Core_3_Igaba_a_tau[0]));
		outfile__array_Core_3_Igaba_a_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_a_tau." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_a_w0;
	outfile__array_Core_3_Igaba_a_w0.open("results/_array_Core_3_Igaba_a_w0_8547859044992465929", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_a_w0.is_open())
	{
		outfile__array_Core_3_Igaba_a_w0.write(reinterpret_cast<char*>(_array_Core_3_Igaba_a_w0), 256*sizeof(_array_Core_3_Igaba_a_w0[0]));
		outfile__array_Core_3_Igaba_a_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_a_w0." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_b;
	outfile__array_Core_3_Igaba_b.open("results/_array_Core_3_Igaba_b_-2140163809950727578", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_b.is_open())
	{
		outfile__array_Core_3_Igaba_b.write(reinterpret_cast<char*>(_array_Core_3_Igaba_b), 256*sizeof(_array_Core_3_Igaba_b[0]));
		outfile__array_Core_3_Igaba_b.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_b." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_b_tau;
	outfile__array_Core_3_Igaba_b_tau.open("results/_array_Core_3_Igaba_b_tau_-4753321281669164666", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_b_tau.is_open())
	{
		outfile__array_Core_3_Igaba_b_tau.write(reinterpret_cast<char*>(_array_Core_3_Igaba_b_tau), 256*sizeof(_array_Core_3_Igaba_b_tau[0]));
		outfile__array_Core_3_Igaba_b_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_b_tau." << endl;
	}
	ofstream outfile__array_Core_3_Igaba_b_w0;
	outfile__array_Core_3_Igaba_b_w0.open("results/_array_Core_3_Igaba_b_w0_-3055420324327214710", ios::binary | ios::out);
	if(outfile__array_Core_3_Igaba_b_w0.is_open())
	{
		outfile__array_Core_3_Igaba_b_w0.write(reinterpret_cast<char*>(_array_Core_3_Igaba_b_w0), 256*sizeof(_array_Core_3_Igaba_b_w0[0]));
		outfile__array_Core_3_Igaba_b_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Igaba_b_w0." << endl;
	}
	ofstream outfile__array_Core_3_Inmda;
	outfile__array_Core_3_Inmda.open("results/_array_Core_3_Inmda_-3903150805585098248", ios::binary | ios::out);
	if(outfile__array_Core_3_Inmda.is_open())
	{
		outfile__array_Core_3_Inmda.write(reinterpret_cast<char*>(_array_Core_3_Inmda), 256*sizeof(_array_Core_3_Inmda[0]));
		outfile__array_Core_3_Inmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Inmda." << endl;
	}
	ofstream outfile__array_Core_3_Inmda_tau;
	outfile__array_Core_3_Inmda_tau.open("results/_array_Core_3_Inmda_tau_4142429705960385365", ios::binary | ios::out);
	if(outfile__array_Core_3_Inmda_tau.is_open())
	{
		outfile__array_Core_3_Inmda_tau.write(reinterpret_cast<char*>(_array_Core_3_Inmda_tau), 256*sizeof(_array_Core_3_Inmda_tau[0]));
		outfile__array_Core_3_Inmda_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Inmda_tau." << endl;
	}
	ofstream outfile__array_Core_3_Inmda_w0;
	outfile__array_Core_3_Inmda_w0.open("results/_array_Core_3_Inmda_w0_8805286254491715353", ios::binary | ios::out);
	if(outfile__array_Core_3_Inmda_w0.is_open())
	{
		outfile__array_Core_3_Inmda_w0.write(reinterpret_cast<char*>(_array_Core_3_Inmda_w0), 256*sizeof(_array_Core_3_Inmda_w0[0]));
		outfile__array_Core_3_Inmda_w0.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Inmda_w0." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_ahp;
	outfile__array_Core_3_Isoma_ahp.open("results/_array_Core_3_Isoma_ahp_-5129403062028320144", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_ahp.is_open())
	{
		outfile__array_Core_3_Isoma_ahp.write(reinterpret_cast<char*>(_array_Core_3_Isoma_ahp), 256*sizeof(_array_Core_3_Isoma_ahp[0]));
		outfile__array_Core_3_Isoma_ahp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_ahp." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_ahp_tau;
	outfile__array_Core_3_Isoma_ahp_tau.open("results/_array_Core_3_Isoma_ahp_tau_6770652563715760761", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_ahp_tau.is_open())
	{
		outfile__array_Core_3_Isoma_ahp_tau.write(reinterpret_cast<char*>(_array_Core_3_Isoma_ahp_tau), 1*sizeof(_array_Core_3_Isoma_ahp_tau[0]));
		outfile__array_Core_3_Isoma_ahp_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_ahp_tau." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_ahp_th;
	outfile__array_Core_3_Isoma_ahp_th.open("results/_array_Core_3_Isoma_ahp_th_1847420288654044548", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_ahp_th.is_open())
	{
		outfile__array_Core_3_Isoma_ahp_th.write(reinterpret_cast<char*>(_array_Core_3_Isoma_ahp_th), 1*sizeof(_array_Core_3_Isoma_ahp_th[0]));
		outfile__array_Core_3_Isoma_ahp_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_ahp_th." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_ahp_w;
	outfile__array_Core_3_Isoma_ahp_w.open("results/_array_Core_3_Isoma_ahp_w_-1453795793844147635", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_ahp_w.is_open())
	{
		outfile__array_Core_3_Isoma_ahp_w.write(reinterpret_cast<char*>(_array_Core_3_Isoma_ahp_w), 256*sizeof(_array_Core_3_Isoma_ahp_w[0]));
		outfile__array_Core_3_Isoma_ahp_w.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_ahp_w." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_const;
	outfile__array_Core_3_Isoma_const.open("results/_array_Core_3_Isoma_const_7084774506719423407", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_const.is_open())
	{
		outfile__array_Core_3_Isoma_const.write(reinterpret_cast<char*>(_array_Core_3_Isoma_const), 256*sizeof(_array_Core_3_Isoma_const[0]));
		outfile__array_Core_3_Isoma_const.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_const." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_dpi_tau;
	outfile__array_Core_3_Isoma_dpi_tau.open("results/_array_Core_3_Isoma_dpi_tau_-7023190833673187737", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_dpi_tau.is_open())
	{
		outfile__array_Core_3_Isoma_dpi_tau.write(reinterpret_cast<char*>(_array_Core_3_Isoma_dpi_tau), 1*sizeof(_array_Core_3_Isoma_dpi_tau[0]));
		outfile__array_Core_3_Isoma_dpi_tau.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_dpi_tau." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_mem;
	outfile__array_Core_3_Isoma_mem.open("results/_array_Core_3_Isoma_mem_6235376650510639338", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_mem.is_open())
	{
		outfile__array_Core_3_Isoma_mem.write(reinterpret_cast<char*>(_array_Core_3_Isoma_mem), 256*sizeof(_array_Core_3_Isoma_mem[0]));
		outfile__array_Core_3_Isoma_mem.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_mem." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_pfb_gain;
	outfile__array_Core_3_Isoma_pfb_gain.open("results/_array_Core_3_Isoma_pfb_gain_5534789968890428236", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_pfb_gain.is_open())
	{
		outfile__array_Core_3_Isoma_pfb_gain.write(reinterpret_cast<char*>(_array_Core_3_Isoma_pfb_gain), 1*sizeof(_array_Core_3_Isoma_pfb_gain[0]));
		outfile__array_Core_3_Isoma_pfb_gain.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_pfb_gain." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_pfb_norm;
	outfile__array_Core_3_Isoma_pfb_norm.open("results/_array_Core_3_Isoma_pfb_norm_-7643615461468460731", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_pfb_norm.is_open())
	{
		outfile__array_Core_3_Isoma_pfb_norm.write(reinterpret_cast<char*>(_array_Core_3_Isoma_pfb_norm), 1*sizeof(_array_Core_3_Isoma_pfb_norm[0]));
		outfile__array_Core_3_Isoma_pfb_norm.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_pfb_norm." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_pfb_th;
	outfile__array_Core_3_Isoma_pfb_th.open("results/_array_Core_3_Isoma_pfb_th_5465390347841482473", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_pfb_th.is_open())
	{
		outfile__array_Core_3_Isoma_pfb_th.write(reinterpret_cast<char*>(_array_Core_3_Isoma_pfb_th), 1*sizeof(_array_Core_3_Isoma_pfb_th[0]));
		outfile__array_Core_3_Isoma_pfb_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_pfb_th." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_reset;
	outfile__array_Core_3_Isoma_reset.open("results/_array_Core_3_Isoma_reset_-7838902228483364023", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_reset.is_open())
	{
		outfile__array_Core_3_Isoma_reset.write(reinterpret_cast<char*>(_array_Core_3_Isoma_reset), 1*sizeof(_array_Core_3_Isoma_reset[0]));
		outfile__array_Core_3_Isoma_reset.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_reset." << endl;
	}
	ofstream outfile__array_Core_3_Isoma_th;
	outfile__array_Core_3_Isoma_th.open("results/_array_Core_3_Isoma_th_5036898745500982380", ios::binary | ios::out);
	if(outfile__array_Core_3_Isoma_th.is_open())
	{
		outfile__array_Core_3_Isoma_th.write(reinterpret_cast<char*>(_array_Core_3_Isoma_th), 1*sizeof(_array_Core_3_Isoma_th[0]));
		outfile__array_Core_3_Isoma_th.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Isoma_th." << endl;
	}
	ofstream outfile__array_Core_3_kn;
	outfile__array_Core_3_kn.open("results/_array_Core_3_kn_-5596949268185192018", ios::binary | ios::out);
	if(outfile__array_Core_3_kn.is_open())
	{
		outfile__array_Core_3_kn.write(reinterpret_cast<char*>(_array_Core_3_kn), 1*sizeof(_array_Core_3_kn[0]));
		outfile__array_Core_3_kn.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_kn." << endl;
	}
	ofstream outfile__array_Core_3_kp;
	outfile__array_Core_3_kp.open("results/_array_Core_3_kp_531274463574687305", ios::binary | ios::out);
	if(outfile__array_Core_3_kp.is_open())
	{
		outfile__array_Core_3_kp.write(reinterpret_cast<char*>(_array_Core_3_kp), 1*sizeof(_array_Core_3_kp[0]));
		outfile__array_Core_3_kp.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_kp." << endl;
	}
	ofstream outfile__array_Core_3_lastspike;
	outfile__array_Core_3_lastspike.open("results/_array_Core_3_lastspike_5828215839358648469", ios::binary | ios::out);
	if(outfile__array_Core_3_lastspike.is_open())
	{
		outfile__array_Core_3_lastspike.write(reinterpret_cast<char*>(_array_Core_3_lastspike), 256*sizeof(_array_Core_3_lastspike[0]));
		outfile__array_Core_3_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_lastspike." << endl;
	}
	ofstream outfile__array_Core_3_not_refractory;
	outfile__array_Core_3_not_refractory.open("results/_array_Core_3_not_refractory_4285103038334197389", ios::binary | ios::out);
	if(outfile__array_Core_3_not_refractory.is_open())
	{
		outfile__array_Core_3_not_refractory.write(reinterpret_cast<char*>(_array_Core_3_not_refractory), 256*sizeof(_array_Core_3_not_refractory[0]));
		outfile__array_Core_3_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_not_refractory." << endl;
	}
	ofstream outfile__array_Core_3_soma_refP;
	outfile__array_Core_3_soma_refP.open("results/_array_Core_3_soma_refP_6738957755782379340", ios::binary | ios::out);
	if(outfile__array_Core_3_soma_refP.is_open())
	{
		outfile__array_Core_3_soma_refP.write(reinterpret_cast<char*>(_array_Core_3_soma_refP), 1*sizeof(_array_Core_3_soma_refP[0]));
		outfile__array_Core_3_soma_refP.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_soma_refP." << endl;
	}
	ofstream outfile__array_Core_3_Ut;
	outfile__array_Core_3_Ut.open("results/_array_Core_3_Ut_-4131228439442331510", ios::binary | ios::out);
	if(outfile__array_Core_3_Ut.is_open())
	{
		outfile__array_Core_3_Ut.write(reinterpret_cast<char*>(_array_Core_3_Ut), 1*sizeof(_array_Core_3_Ut[0]));
		outfile__array_Core_3_Ut.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Ut." << endl;
	}
	ofstream outfile__array_Core_3_Vnmda;
	outfile__array_Core_3_Vnmda.open("results/_array_Core_3_Vnmda_-6892628026604575664", ios::binary | ios::out);
	if(outfile__array_Core_3_Vnmda.is_open())
	{
		outfile__array_Core_3_Vnmda.write(reinterpret_cast<char*>(_array_Core_3_Vnmda), 256*sizeof(_array_Core_3_Vnmda[0]));
		outfile__array_Core_3_Vnmda.close();
	} else
	{
		std::cout << "Error writing output file for _array_Core_3_Vnmda." << endl;
	}
	ofstream outfile__array_defaultclock_dt;
	outfile__array_defaultclock_dt.open("results/_array_defaultclock_dt_-8325147453693480284", ios::binary | ios::out);
	if(outfile__array_defaultclock_dt.is_open())
	{
		outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
		outfile__array_defaultclock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
	}
	ofstream outfile__array_defaultclock_t;
	outfile__array_defaultclock_t.open("results/_array_defaultclock_t_3766643217606192898", ios::binary | ios::out);
	if(outfile__array_defaultclock_t.is_open())
	{
		outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
		outfile__array_defaultclock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_t." << endl;
	}
	ofstream outfile__array_defaultclock_timestep;
	outfile__array_defaultclock_timestep.open("results/_array_defaultclock_timestep_3459174420874988273", ios::binary | ios::out);
	if(outfile__array_defaultclock_timestep.is_open())
	{
		outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
		outfile__array_defaultclock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
	}
	ofstream outfile__array_mon_neuron_input__source_idx;
	outfile__array_mon_neuron_input__source_idx.open("results/_array_mon_neuron_input__source_idx_-1747735418370144032", ios::binary | ios::out);
	if(outfile__array_mon_neuron_input__source_idx.is_open())
	{
		outfile__array_mon_neuron_input__source_idx.write(reinterpret_cast<char*>(_array_mon_neuron_input__source_idx), 200*sizeof(_array_mon_neuron_input__source_idx[0]));
		outfile__array_mon_neuron_input__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_input__source_idx." << endl;
	}
	ofstream outfile__array_mon_neuron_input_count;
	outfile__array_mon_neuron_input_count.open("results/_array_mon_neuron_input_count_135422840089586940", ios::binary | ios::out);
	if(outfile__array_mon_neuron_input_count.is_open())
	{
		outfile__array_mon_neuron_input_count.write(reinterpret_cast<char*>(_array_mon_neuron_input_count), 200*sizeof(_array_mon_neuron_input_count[0]));
		outfile__array_mon_neuron_input_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_input_count." << endl;
	}
	ofstream outfile__array_mon_neuron_input_N;
	outfile__array_mon_neuron_input_N.open("results/_array_mon_neuron_input_N_3839155734707027753", ios::binary | ios::out);
	if(outfile__array_mon_neuron_input_N.is_open())
	{
		outfile__array_mon_neuron_input_N.write(reinterpret_cast<char*>(_array_mon_neuron_input_N), 1*sizeof(_array_mon_neuron_input_N[0]));
		outfile__array_mon_neuron_input_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_input_N." << endl;
	}
	ofstream outfile__array_mon_neuron_output__source_idx;
	outfile__array_mon_neuron_output__source_idx.open("results/_array_mon_neuron_output__source_idx_6806756496035578946", ios::binary | ios::out);
	if(outfile__array_mon_neuron_output__source_idx.is_open())
	{
		outfile__array_mon_neuron_output__source_idx.write(reinterpret_cast<char*>(_array_mon_neuron_output__source_idx), 200*sizeof(_array_mon_neuron_output__source_idx[0]));
		outfile__array_mon_neuron_output__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_output__source_idx." << endl;
	}
	ofstream outfile__array_mon_neuron_output_count;
	outfile__array_mon_neuron_output_count.open("results/_array_mon_neuron_output_count_542873530678187427", ios::binary | ios::out);
	if(outfile__array_mon_neuron_output_count.is_open())
	{
		outfile__array_mon_neuron_output_count.write(reinterpret_cast<char*>(_array_mon_neuron_output_count), 200*sizeof(_array_mon_neuron_output_count[0]));
		outfile__array_mon_neuron_output_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_output_count." << endl;
	}
	ofstream outfile__array_mon_neuron_output_N;
	outfile__array_mon_neuron_output_N.open("results/_array_mon_neuron_output_N_-1063235664358423436", ios::binary | ios::out);
	if(outfile__array_mon_neuron_output_N.is_open())
	{
		outfile__array_mon_neuron_output_N.write(reinterpret_cast<char*>(_array_mon_neuron_output_N), 1*sizeof(_array_mon_neuron_output_N[0]));
		outfile__array_mon_neuron_output_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_mon_neuron_output_N." << endl;
	}
	ofstream outfile__array_NMDA0_N;
	outfile__array_NMDA0_N.open("results/_array_NMDA0_N_-5284015481342164905", ios::binary | ios::out);
	if(outfile__array_NMDA0_N.is_open())
	{
		outfile__array_NMDA0_N.write(reinterpret_cast<char*>(_array_NMDA0_N), 1*sizeof(_array_NMDA0_N[0]));
		outfile__array_NMDA0_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_NMDA0_N." << endl;
	}
	ofstream outfile__array_poissongroup__spikespace;
	outfile__array_poissongroup__spikespace.open("results/_array_poissongroup__spikespace_-2645669618090137782", ios::binary | ios::out);
	if(outfile__array_poissongroup__spikespace.is_open())
	{
		outfile__array_poissongroup__spikespace.write(reinterpret_cast<char*>(_array_poissongroup__spikespace), 201*sizeof(_array_poissongroup__spikespace[0]));
		outfile__array_poissongroup__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_poissongroup__spikespace." << endl;
	}
	ofstream outfile__array_poissongroup_i;
	outfile__array_poissongroup_i.open("results/_array_poissongroup_i_-4574790696377929542", ios::binary | ios::out);
	if(outfile__array_poissongroup_i.is_open())
	{
		outfile__array_poissongroup_i.write(reinterpret_cast<char*>(_array_poissongroup_i), 200*sizeof(_array_poissongroup_i[0]));
		outfile__array_poissongroup_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_poissongroup_i." << endl;
	}

	ofstream outfile__dynamic_array_mon_neuron_input_i;
	outfile__dynamic_array_mon_neuron_input_i.open("results/_dynamic_array_mon_neuron_input_i_6675097172020364218", ios::binary | ios::out);
	if(outfile__dynamic_array_mon_neuron_input_i.is_open())
	{
        if (! _dynamic_array_mon_neuron_input_i.empty() )
        {
			outfile__dynamic_array_mon_neuron_input_i.write(reinterpret_cast<char*>(&_dynamic_array_mon_neuron_input_i[0]), _dynamic_array_mon_neuron_input_i.size()*sizeof(_dynamic_array_mon_neuron_input_i[0]));
		    outfile__dynamic_array_mon_neuron_input_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_mon_neuron_input_i." << endl;
	}
	ofstream outfile__dynamic_array_mon_neuron_input_t;
	outfile__dynamic_array_mon_neuron_input_t.open("results/_dynamic_array_mon_neuron_input_t_565560797149160085", ios::binary | ios::out);
	if(outfile__dynamic_array_mon_neuron_input_t.is_open())
	{
        if (! _dynamic_array_mon_neuron_input_t.empty() )
        {
			outfile__dynamic_array_mon_neuron_input_t.write(reinterpret_cast<char*>(&_dynamic_array_mon_neuron_input_t[0]), _dynamic_array_mon_neuron_input_t.size()*sizeof(_dynamic_array_mon_neuron_input_t[0]));
		    outfile__dynamic_array_mon_neuron_input_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_mon_neuron_input_t." << endl;
	}
	ofstream outfile__dynamic_array_mon_neuron_output_i;
	outfile__dynamic_array_mon_neuron_output_i.open("results/_dynamic_array_mon_neuron_output_i_8122066792676967782", ios::binary | ios::out);
	if(outfile__dynamic_array_mon_neuron_output_i.is_open())
	{
        if (! _dynamic_array_mon_neuron_output_i.empty() )
        {
			outfile__dynamic_array_mon_neuron_output_i.write(reinterpret_cast<char*>(&_dynamic_array_mon_neuron_output_i[0]), _dynamic_array_mon_neuron_output_i.size()*sizeof(_dynamic_array_mon_neuron_output_i[0]));
		    outfile__dynamic_array_mon_neuron_output_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_mon_neuron_output_i." << endl;
	}
	ofstream outfile__dynamic_array_mon_neuron_output_t;
	outfile__dynamic_array_mon_neuron_output_t.open("results/_dynamic_array_mon_neuron_output_t_-1126661221566807506", ios::binary | ios::out);
	if(outfile__dynamic_array_mon_neuron_output_t.is_open())
	{
        if (! _dynamic_array_mon_neuron_output_t.empty() )
        {
			outfile__dynamic_array_mon_neuron_output_t.write(reinterpret_cast<char*>(&_dynamic_array_mon_neuron_output_t[0]), _dynamic_array_mon_neuron_output_t.size()*sizeof(_dynamic_array_mon_neuron_output_t[0]));
		    outfile__dynamic_array_mon_neuron_output_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_mon_neuron_output_t." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0__synaptic_post;
	outfile__dynamic_array_NMDA0__synaptic_post.open("results/_dynamic_array_NMDA0__synaptic_post_122494412421897090", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0__synaptic_post.is_open())
	{
        if (! _dynamic_array_NMDA0__synaptic_post.empty() )
        {
			outfile__dynamic_array_NMDA0__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0__synaptic_post[0]), _dynamic_array_NMDA0__synaptic_post.size()*sizeof(_dynamic_array_NMDA0__synaptic_post[0]));
		    outfile__dynamic_array_NMDA0__synaptic_post.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0__synaptic_post." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0__synaptic_pre;
	outfile__dynamic_array_NMDA0__synaptic_pre.open("results/_dynamic_array_NMDA0__synaptic_pre_3868825687735666305", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0__synaptic_pre.is_open())
	{
        if (! _dynamic_array_NMDA0__synaptic_pre.empty() )
        {
			outfile__dynamic_array_NMDA0__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0__synaptic_pre[0]), _dynamic_array_NMDA0__synaptic_pre.size()*sizeof(_dynamic_array_NMDA0__synaptic_pre[0]));
		    outfile__dynamic_array_NMDA0__synaptic_pre.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0__synaptic_pre." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0_delay;
	outfile__dynamic_array_NMDA0_delay.open("results/_dynamic_array_NMDA0_delay_2099685013340918456", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0_delay.is_open())
	{
        if (! _dynamic_array_NMDA0_delay.empty() )
        {
			outfile__dynamic_array_NMDA0_delay.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0_delay[0]), _dynamic_array_NMDA0_delay.size()*sizeof(_dynamic_array_NMDA0_delay[0]));
		    outfile__dynamic_array_NMDA0_delay.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0_delay." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0_delay_1;
	outfile__dynamic_array_NMDA0_delay_1.open("results/_dynamic_array_NMDA0_delay_1_539988439311615425", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0_delay_1.is_open())
	{
        if (! _dynamic_array_NMDA0_delay_1.empty() )
        {
			outfile__dynamic_array_NMDA0_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0_delay_1[0]), _dynamic_array_NMDA0_delay_1.size()*sizeof(_dynamic_array_NMDA0_delay_1[0]));
		    outfile__dynamic_array_NMDA0_delay_1.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0_delay_1." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0_N_incoming;
	outfile__dynamic_array_NMDA0_N_incoming.open("results/_dynamic_array_NMDA0_N_incoming_-622646284535003387", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0_N_incoming.is_open())
	{
        if (! _dynamic_array_NMDA0_N_incoming.empty() )
        {
			outfile__dynamic_array_NMDA0_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0_N_incoming[0]), _dynamic_array_NMDA0_N_incoming.size()*sizeof(_dynamic_array_NMDA0_N_incoming[0]));
		    outfile__dynamic_array_NMDA0_N_incoming.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0_N_incoming." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0_N_outgoing;
	outfile__dynamic_array_NMDA0_N_outgoing.open("results/_dynamic_array_NMDA0_N_outgoing_5185441601732951842", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0_N_outgoing.is_open())
	{
        if (! _dynamic_array_NMDA0_N_outgoing.empty() )
        {
			outfile__dynamic_array_NMDA0_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0_N_outgoing[0]), _dynamic_array_NMDA0_N_outgoing.size()*sizeof(_dynamic_array_NMDA0_N_outgoing[0]));
		    outfile__dynamic_array_NMDA0_N_outgoing.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0_N_outgoing." << endl;
	}
	ofstream outfile__dynamic_array_NMDA0_weight;
	outfile__dynamic_array_NMDA0_weight.open("results/_dynamic_array_NMDA0_weight_-3239081387640548461", ios::binary | ios::out);
	if(outfile__dynamic_array_NMDA0_weight.is_open())
	{
        if (! _dynamic_array_NMDA0_weight.empty() )
        {
			outfile__dynamic_array_NMDA0_weight.write(reinterpret_cast<char*>(&_dynamic_array_NMDA0_weight[0]), _dynamic_array_NMDA0_weight.size()*sizeof(_dynamic_array_NMDA0_weight[0]));
		    outfile__dynamic_array_NMDA0_weight.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_NMDA0_weight." << endl;
	}

	// Write last run info to disk
	ofstream outfile_last_run_info;
	outfile_last_run_info.open("results/last_run_info.txt", ios::out);
	if(outfile_last_run_info.is_open())
	{
		outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
		outfile_last_run_info.close();
	} else
	{
	    std::cout << "Error writing last run info to file." << std::endl;
	}
}

void _dealloc_arrays()
{
	using namespace brian;


	// static arrays
	if(_timedarray_values!=0)
	{
		delete [] _timedarray_values;
		_timedarray_values = 0;
	}
}

