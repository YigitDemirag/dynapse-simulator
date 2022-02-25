#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>

#include "run.h"
#include "brianlib/common_math.h"
#include "randomkit.h"

#include "code_objects/AMPA0_post_codeobject.h"
#include "code_objects/AMPA0_post_push_spikes.h"
#include "code_objects/before_run_AMPA0_post_push_spikes.h"
#include "code_objects/AMPA0_pre_codeobject.h"
#include "code_objects/AMPA0_pre_push_spikes.h"
#include "code_objects/before_run_AMPA0_pre_push_spikes.h"
#include "code_objects/AMPA0_synapses_create_generator_codeobject.h"
#include "code_objects/Core_0_resetter_1_codeobject.h"
#include "code_objects/Core_0_stateupdater_1_codeobject.h"
#include "code_objects/Core_0_thresholder_1_codeobject.h"
#include "code_objects/Core_1_resetter_1_codeobject.h"
#include "code_objects/Core_1_stateupdater_1_codeobject.h"
#include "code_objects/Core_1_thresholder_1_codeobject.h"
#include "code_objects/Core_2_resetter_1_codeobject.h"
#include "code_objects/Core_2_stateupdater_1_codeobject.h"
#include "code_objects/Core_2_thresholder_1_codeobject.h"
#include "code_objects/Core_3_resetter_1_codeobject.h"
#include "code_objects/Core_3_stateupdater_1_codeobject.h"
#include "code_objects/Core_3_thresholder_1_codeobject.h"
#include "code_objects/InpSpikeGenerator_codeobject.h"
#include "code_objects/mon_neuron_input_codeobject.h"
#include "code_objects/mon_neuron_output_codeobject.h"
#include "code_objects/statemonitor_1_codeobject.h"
#include "code_objects/statemonitor_2_codeobject.h"


#include <iostream>
#include <fstream>
#include <string>




int main(int argc, char **argv)
{
        

	brian_start();
        

	{
		using namespace brian;

		
                
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 1.9999999999999998e-05;
        _dynamic_array_InpSpikeGenerator_spike_number.resize(500);
        
                        
                        for(int i=0; i<_dynamic_array_InpSpikeGenerator_spike_number.size(); i++)
                        {
                            _dynamic_array_InpSpikeGenerator_spike_number[i] = _static_array__dynamic_array_InpSpikeGenerator_spike_number[i];
                        }
                        
        _dynamic_array_InpSpikeGenerator_neuron_index.resize(500);
        
                        
                        for(int i=0; i<_dynamic_array_InpSpikeGenerator_neuron_index.size(); i++)
                        {
                            _dynamic_array_InpSpikeGenerator_neuron_index[i] = _static_array__dynamic_array_InpSpikeGenerator_neuron_index[i];
                        }
                        
        _dynamic_array_InpSpikeGenerator_spike_time.resize(500);
        
                        
                        for(int i=0; i<_dynamic_array_InpSpikeGenerator_spike_time.size(); i++)
                        {
                            _dynamic_array_InpSpikeGenerator_spike_time[i] = _static_array__dynamic_array_InpSpikeGenerator_spike_time[i];
                        }
                        
        _dynamic_array_InpSpikeGenerator__timebins.resize(500);
        _array_InpSpikeGenerator__lastindex[0] = 0;
        _array_InpSpikeGenerator_period[0] = 0.0;
        
                        
                        for(int i=0; i<_num__array_Core_0_lastspike; i++)
                        {
                            _array_Core_0_lastspike[i] = - 10000.0;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_not_refractory; i++)
                        {
                            _array_Core_0_not_refractory[i] = true;
                        }
                        
        _array_Core_0_kn[0] = 0.75;
        _array_Core_0_kp[0] = 0.66;
        _array_Core_0_Ut[0] = 0.025;
        _array_Core_0_I0[0] = 5e-13;
        _array_Core_0_alpha[0] = 4;
        _array_Core_0_Csoma_mem[0] = 1.5e-12;
        _array_Core_0_Isoma_dpi_tau[0] = 3e-12;
        _array_Core_0_Isoma_th[0] = 1e-09;
        _array_Core_0_Isoma_reset[0] = 6e-13;
        
                        
                        for(int i=0; i<_num__array_Core_0_Isoma_const; i++)
                        {
                            _array_Core_0_Isoma_const[i] = 5e-13;
                        }
                        
        _array_Core_0_soma_refP[0] = 0.005;
        _array_Core_0_Csoma_ahp[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_0_Isoma_ahp; i++)
                        {
                            _array_Core_0_Isoma_ahp[i] = 5e-13;
                        }
                        
        _array_Core_0_Isoma_ahp_tau[0] = 1e-12;
        _array_Core_0_Isoma_ahp_th[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_0_Isoma_ahp_w; i++)
                        {
                            _array_Core_0_Isoma_ahp_w[i] = 2e-12;
                        }
                        
        _array_Core_0_Isoma_pfb_gain[0] = 5e-11;
        _array_Core_0_Isoma_pfb_th[0] = 5e-10;
        _array_Core_0_Isoma_pfb_norm[0] = 1e-11;
        
                        
                        for(int i=0; i<_num__array_Core_0_Cnmda; i++)
                        {
                            _array_Core_0_Cnmda[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Inmda_tau; i++)
                        {
                            _array_Core_0_Inmda_tau[i] = 1e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Inmda_w0; i++)
                        {
                            _array_Core_0_Inmda_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Inmda; i++)
                        {
                            _array_Core_0_Inmda[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Vnmda; i++)
                        {
                            _array_Core_0_Vnmda[i] = 0.01;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Campa; i++)
                        {
                            _array_Core_0_Campa[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Iampa_tau; i++)
                        {
                            _array_Core_0_Iampa_tau[i] = 6e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Iampa_w0; i++)
                        {
                            _array_Core_0_Iampa_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Iampa; i++)
                        {
                            _array_Core_0_Iampa[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Cgaba_b; i++)
                        {
                            _array_Core_0_Cgaba_b[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_b_tau; i++)
                        {
                            _array_Core_0_Igaba_b_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_b_w0; i++)
                        {
                            _array_Core_0_Igaba_b_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_b; i++)
                        {
                            _array_Core_0_Igaba_b[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Cgaba_a; i++)
                        {
                            _array_Core_0_Cgaba_a[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_a_tau; i++)
                        {
                            _array_Core_0_Igaba_a_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_a_w0; i++)
                        {
                            _array_Core_0_Igaba_a_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_0_Igaba_a; i++)
                        {
                            _array_Core_0_Igaba_a[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_lastspike; i++)
                        {
                            _array_Core_1_lastspike[i] = - 10000.0;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_not_refractory; i++)
                        {
                            _array_Core_1_not_refractory[i] = true;
                        }
                        
        _array_Core_1_kn[0] = 0.75;
        _array_Core_1_kp[0] = 0.66;
        _array_Core_1_Ut[0] = 0.025;
        _array_Core_1_I0[0] = 5e-13;
        _array_Core_1_alpha[0] = 4;
        _array_Core_1_Csoma_mem[0] = 1.5e-12;
        _array_Core_1_Isoma_dpi_tau[0] = 3e-12;
        _array_Core_1_Isoma_th[0] = 1e-09;
        _array_Core_1_Isoma_reset[0] = 6e-13;
        
                        
                        for(int i=0; i<_num__array_Core_1_Isoma_const; i++)
                        {
                            _array_Core_1_Isoma_const[i] = 5e-13;
                        }
                        
        _array_Core_1_soma_refP[0] = 0.005;
        _array_Core_1_Csoma_ahp[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_1_Isoma_ahp; i++)
                        {
                            _array_Core_1_Isoma_ahp[i] = 5e-13;
                        }
                        
        _array_Core_1_Isoma_ahp_tau[0] = 1e-12;
        _array_Core_1_Isoma_ahp_th[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_1_Isoma_ahp_w; i++)
                        {
                            _array_Core_1_Isoma_ahp_w[i] = 2e-12;
                        }
                        
        _array_Core_1_Isoma_pfb_gain[0] = 5e-11;
        _array_Core_1_Isoma_pfb_th[0] = 5e-10;
        _array_Core_1_Isoma_pfb_norm[0] = 1e-11;
        
                        
                        for(int i=0; i<_num__array_Core_1_Cnmda; i++)
                        {
                            _array_Core_1_Cnmda[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Inmda_tau; i++)
                        {
                            _array_Core_1_Inmda_tau[i] = 1e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Inmda_w0; i++)
                        {
                            _array_Core_1_Inmda_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Inmda; i++)
                        {
                            _array_Core_1_Inmda[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Vnmda; i++)
                        {
                            _array_Core_1_Vnmda[i] = 0.01;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Campa; i++)
                        {
                            _array_Core_1_Campa[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Iampa_tau; i++)
                        {
                            _array_Core_1_Iampa_tau[i] = 6e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Iampa_w0; i++)
                        {
                            _array_Core_1_Iampa_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Iampa; i++)
                        {
                            _array_Core_1_Iampa[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Cgaba_b; i++)
                        {
                            _array_Core_1_Cgaba_b[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_b_tau; i++)
                        {
                            _array_Core_1_Igaba_b_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_b_w0; i++)
                        {
                            _array_Core_1_Igaba_b_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_b; i++)
                        {
                            _array_Core_1_Igaba_b[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Cgaba_a; i++)
                        {
                            _array_Core_1_Cgaba_a[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_a_tau; i++)
                        {
                            _array_Core_1_Igaba_a_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_a_w0; i++)
                        {
                            _array_Core_1_Igaba_a_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_1_Igaba_a; i++)
                        {
                            _array_Core_1_Igaba_a[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_lastspike; i++)
                        {
                            _array_Core_2_lastspike[i] = - 10000.0;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_not_refractory; i++)
                        {
                            _array_Core_2_not_refractory[i] = true;
                        }
                        
        _array_Core_2_kn[0] = 0.75;
        _array_Core_2_kp[0] = 0.66;
        _array_Core_2_Ut[0] = 0.025;
        _array_Core_2_I0[0] = 5e-13;
        _array_Core_2_alpha[0] = 4;
        _array_Core_2_Csoma_mem[0] = 1.5e-12;
        _array_Core_2_Isoma_dpi_tau[0] = 3e-12;
        _array_Core_2_Isoma_th[0] = 1e-09;
        _array_Core_2_Isoma_reset[0] = 6e-13;
        
                        
                        for(int i=0; i<_num__array_Core_2_Isoma_const; i++)
                        {
                            _array_Core_2_Isoma_const[i] = 5e-13;
                        }
                        
        _array_Core_2_soma_refP[0] = 0.005;
        _array_Core_2_Csoma_ahp[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_2_Isoma_ahp; i++)
                        {
                            _array_Core_2_Isoma_ahp[i] = 5e-13;
                        }
                        
        _array_Core_2_Isoma_ahp_tau[0] = 1e-12;
        _array_Core_2_Isoma_ahp_th[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_2_Isoma_ahp_w; i++)
                        {
                            _array_Core_2_Isoma_ahp_w[i] = 2e-12;
                        }
                        
        _array_Core_2_Isoma_pfb_gain[0] = 5e-11;
        _array_Core_2_Isoma_pfb_th[0] = 5e-10;
        _array_Core_2_Isoma_pfb_norm[0] = 1e-11;
        
                        
                        for(int i=0; i<_num__array_Core_2_Cnmda; i++)
                        {
                            _array_Core_2_Cnmda[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Inmda_tau; i++)
                        {
                            _array_Core_2_Inmda_tau[i] = 1e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Inmda_w0; i++)
                        {
                            _array_Core_2_Inmda_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Inmda; i++)
                        {
                            _array_Core_2_Inmda[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Vnmda; i++)
                        {
                            _array_Core_2_Vnmda[i] = 0.01;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Campa; i++)
                        {
                            _array_Core_2_Campa[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Iampa_tau; i++)
                        {
                            _array_Core_2_Iampa_tau[i] = 6e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Iampa_w0; i++)
                        {
                            _array_Core_2_Iampa_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Iampa; i++)
                        {
                            _array_Core_2_Iampa[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Cgaba_b; i++)
                        {
                            _array_Core_2_Cgaba_b[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_b_tau; i++)
                        {
                            _array_Core_2_Igaba_b_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_b_w0; i++)
                        {
                            _array_Core_2_Igaba_b_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_b; i++)
                        {
                            _array_Core_2_Igaba_b[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Cgaba_a; i++)
                        {
                            _array_Core_2_Cgaba_a[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_a_tau; i++)
                        {
                            _array_Core_2_Igaba_a_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_a_w0; i++)
                        {
                            _array_Core_2_Igaba_a_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_2_Igaba_a; i++)
                        {
                            _array_Core_2_Igaba_a[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_lastspike; i++)
                        {
                            _array_Core_3_lastspike[i] = - 10000.0;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_not_refractory; i++)
                        {
                            _array_Core_3_not_refractory[i] = true;
                        }
                        
        _array_Core_3_kn[0] = 0.75;
        _array_Core_3_kp[0] = 0.66;
        _array_Core_3_Ut[0] = 0.025;
        _array_Core_3_I0[0] = 5e-13;
        _array_Core_3_alpha[0] = 4;
        _array_Core_3_Csoma_mem[0] = 1.5e-12;
        _array_Core_3_Isoma_dpi_tau[0] = 3e-12;
        _array_Core_3_Isoma_th[0] = 1e-09;
        _array_Core_3_Isoma_reset[0] = 6e-13;
        
                        
                        for(int i=0; i<_num__array_Core_3_Isoma_const; i++)
                        {
                            _array_Core_3_Isoma_const[i] = 5e-13;
                        }
                        
        _array_Core_3_soma_refP[0] = 0.005;
        _array_Core_3_Csoma_ahp[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_3_Isoma_ahp; i++)
                        {
                            _array_Core_3_Isoma_ahp[i] = 5e-13;
                        }
                        
        _array_Core_3_Isoma_ahp_tau[0] = 1e-12;
        _array_Core_3_Isoma_ahp_th[0] = 1e-12;
        
                        
                        for(int i=0; i<_num__array_Core_3_Isoma_ahp_w; i++)
                        {
                            _array_Core_3_Isoma_ahp_w[i] = 2e-12;
                        }
                        
        _array_Core_3_Isoma_pfb_gain[0] = 5e-11;
        _array_Core_3_Isoma_pfb_th[0] = 5e-10;
        _array_Core_3_Isoma_pfb_norm[0] = 1e-11;
        
                        
                        for(int i=0; i<_num__array_Core_3_Cnmda; i++)
                        {
                            _array_Core_3_Cnmda[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Inmda_tau; i++)
                        {
                            _array_Core_3_Inmda_tau[i] = 1e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Inmda_w0; i++)
                        {
                            _array_Core_3_Inmda_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Inmda; i++)
                        {
                            _array_Core_3_Inmda[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Vnmda; i++)
                        {
                            _array_Core_3_Vnmda[i] = 0.01;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Campa; i++)
                        {
                            _array_Core_3_Campa[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Iampa_tau; i++)
                        {
                            _array_Core_3_Iampa_tau[i] = 6e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Iampa_w0; i++)
                        {
                            _array_Core_3_Iampa_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Iampa; i++)
                        {
                            _array_Core_3_Iampa[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Cgaba_b; i++)
                        {
                            _array_Core_3_Cgaba_b[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_b_tau; i++)
                        {
                            _array_Core_3_Igaba_b_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_b_w0; i++)
                        {
                            _array_Core_3_Igaba_b_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_b; i++)
                        {
                            _array_Core_3_Igaba_b[i] = 5e-13;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Cgaba_a; i++)
                        {
                            _array_Core_3_Cgaba_a[i] = 1.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_a_tau; i++)
                        {
                            _array_Core_3_Igaba_a_tau[i] = 2.5e-12;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_a_w0; i++)
                        {
                            _array_Core_3_Igaba_a_w0[i] = 5e-11;
                        }
                        
        
                        
                        for(int i=0; i<_num__array_Core_3_Igaba_a; i++)
                        {
                            _array_Core_3_Igaba_a[i] = 5e-13;
                        }
                        
        _run_AMPA0_synapses_create_generator_codeobject();
        
                        
                        for(int i=0; i<_dynamic_array_AMPA0_weight.size(); i++)
                        {
                            _dynamic_array_AMPA0_weight[i] = 500;
                        }
                        
        _array_Core_1_Isoma_ahp_tau[0] = 4e-14;
        _array_Core_1_Isoma_ahp_th[0] = 5e-14;
        
                        
                        for(int i=0; i<_num__array_Core_1_Iampa_tau; i++)
                        {
                            _array_Core_1_Iampa_tau[i] = 1e-11;
                        }
                        
        _array_statemonitor_1__indices[0] = 0;
        _array_statemonitor_2__indices[0] = 0;
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_t[0] = 0.0;
        _array_InpSpikeGenerator__lastindex[0] = 0;
        
                        
                        for(int i=0; i<_dynamic_array_InpSpikeGenerator__timebins.size(); i++)
                        {
                            _dynamic_array_InpSpikeGenerator__timebins[i] = _static_array__dynamic_array_InpSpikeGenerator__timebins[i];
                        }
                        
        _array_InpSpikeGenerator__period_bins[0] = 0.0;
        _before_run_AMPA0_pre_push_spikes();
        _before_run_AMPA0_post_push_spikes();
        network_1.clear();
        network_1.add(&defaultclock, _run_statemonitor_1_codeobject);
        network_1.add(&defaultclock, _run_statemonitor_2_codeobject);
        network_1.add(&defaultclock, _run_Core_0_stateupdater_1_codeobject);
        network_1.add(&defaultclock, _run_Core_1_stateupdater_1_codeobject);
        network_1.add(&defaultclock, _run_Core_2_stateupdater_1_codeobject);
        network_1.add(&defaultclock, _run_Core_3_stateupdater_1_codeobject);
        network_1.add(&defaultclock, _run_Core_0_thresholder_1_codeobject);
        network_1.add(&defaultclock, _run_Core_1_thresholder_1_codeobject);
        network_1.add(&defaultclock, _run_Core_2_thresholder_1_codeobject);
        network_1.add(&defaultclock, _run_Core_3_thresholder_1_codeobject);
        network_1.add(&defaultclock, _run_InpSpikeGenerator_codeobject);
        network_1.add(&defaultclock, _run_mon_neuron_input_codeobject);
        network_1.add(&defaultclock, _run_mon_neuron_output_codeobject);
        network_1.add(&defaultclock, _run_AMPA0_pre_push_spikes);
        network_1.add(&defaultclock, _run_AMPA0_pre_codeobject);
        network_1.add(&defaultclock, _run_AMPA0_post_push_spikes);
        network_1.add(&defaultclock, _run_AMPA0_post_codeobject);
        network_1.add(&defaultclock, _run_Core_0_resetter_1_codeobject);
        network_1.add(&defaultclock, _run_Core_1_resetter_1_codeobject);
        network_1.add(&defaultclock, _run_Core_2_resetter_1_codeobject);
        network_1.add(&defaultclock, _run_Core_3_resetter_1_codeobject);
        network_1.run(5.0, NULL, 10.0);
        #ifdef DEBUG
        _debugmsg_mon_neuron_input_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_mon_neuron_output_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_AMPA0_pre_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_AMPA0_post_codeobject();
        #endif

	}
        

	brian_end();
        

	return 0;
}