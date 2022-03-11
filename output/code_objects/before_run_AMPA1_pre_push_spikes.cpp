#include "code_objects/before_run_AMPA1_pre_push_spikes.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>

////// SUPPORT CODE ///////
namespace {
        
    template < typename T1, typename T2 > struct _higher_type;
    template < > struct _higher_type<int,int> { typedef int type; };
    template < > struct _higher_type<int,long> { typedef long type; };
    template < > struct _higher_type<int,long long> { typedef long long type; };
    template < > struct _higher_type<int,float> { typedef float type; };
    template < > struct _higher_type<int,double> { typedef double type; };
    template < > struct _higher_type<int,long double> { typedef long double type; };
    template < > struct _higher_type<long,int> { typedef long type; };
    template < > struct _higher_type<long,long> { typedef long type; };
    template < > struct _higher_type<long,long long> { typedef long long type; };
    template < > struct _higher_type<long,float> { typedef float type; };
    template < > struct _higher_type<long,double> { typedef double type; };
    template < > struct _higher_type<long,long double> { typedef long double type; };
    template < > struct _higher_type<long long,int> { typedef long long type; };
    template < > struct _higher_type<long long,long> { typedef long long type; };
    template < > struct _higher_type<long long,long long> { typedef long long type; };
    template < > struct _higher_type<long long,float> { typedef float type; };
    template < > struct _higher_type<long long,double> { typedef double type; };
    template < > struct _higher_type<long long,long double> { typedef long double type; };
    template < > struct _higher_type<float,int> { typedef float type; };
    template < > struct _higher_type<float,long> { typedef float type; };
    template < > struct _higher_type<float,long long> { typedef float type; };
    template < > struct _higher_type<float,float> { typedef float type; };
    template < > struct _higher_type<float,double> { typedef double type; };
    template < > struct _higher_type<float,long double> { typedef long double type; };
    template < > struct _higher_type<double,int> { typedef double type; };
    template < > struct _higher_type<double,long> { typedef double type; };
    template < > struct _higher_type<double,long long> { typedef double type; };
    template < > struct _higher_type<double,float> { typedef double type; };
    template < > struct _higher_type<double,double> { typedef double type; };
    template < > struct _higher_type<double,long double> { typedef long double type; };
    template < > struct _higher_type<long double,int> { typedef long double type; };
    template < > struct _higher_type<long double,long> { typedef long double type; };
    template < > struct _higher_type<long double,long long> { typedef long double type; };
    template < > struct _higher_type<long double,float> { typedef long double type; };
    template < > struct _higher_type<long double,double> { typedef long double type; };
    template < > struct _higher_type<long double,long double> { typedef long double type; };
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_mod(T1 x, T2 y)
    {{
        return x-y*floor(1.0*x/y);
    }}
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_floordiv(T1 x, T2 y)
    {{
        return floor(1.0*x/y);
    }}
    #ifdef _MSC_VER
    #define _brian_pow(x, y) (pow((double)(x), (y)))
    #else
    #define _brian_pow(x, y) (pow((x), (y)))
    #endif

}

void _before_run_AMPA1_pre_push_spikes()
{
    using namespace brian;
    ///// CONSTANTS ///////////
    const size_t _num_source_dt = 1;
const int64_t _n_sources = 1;
double* const _array_AMPA1_delay = _dynamic_array_AMPA1_delay.empty()? 0 : &_dynamic_array_AMPA1_delay[0];
const size_t _numdelay = _dynamic_array_AMPA1_delay.size();
const int64_t _n_targets = 1;
const size_t _num_spikespace = 2;
    ///// POINTERS ////////////
        
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_AMPA1_delay = _array_AMPA1_delay;
    int32_t* __restrict  _ptr_array_InpSpikeGenerator__spikespace = _array_InpSpikeGenerator__spikespace;

    std::vector<double> &real_delays = _dynamic_array_AMPA1_delay;
    double* real_delays_data = real_delays.empty() ? 0 : &(real_delays[0]);
    int32_t* sources = AMPA1_pre.sources.empty() ? 0 : &(AMPA1_pre.sources[0]);
    const size_t n_delays = real_delays.size();
    const size_t n_synapses = AMPA1_pre.sources.size();
    AMPA1_pre.prepare(_n_sources,
                           _n_targets,
                           real_delays_data, n_delays, sources,
                           n_synapses,
                           _ptr_array_defaultclock_dt[0]);

}
