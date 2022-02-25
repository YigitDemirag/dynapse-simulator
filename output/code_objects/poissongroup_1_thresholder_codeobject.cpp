#include "code_objects/poissongroup_1_thresholder_codeobject.h"
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
        
    static double* _namespace_timedarray_1_values;
    static inline double _timedarray_1(const double t, const int i)
    {
        const double epsilon = 0.000001000000000000 / 1;
        if (i < 0 || i >= 200)
            return NAN;
        int timestep = (int)((t/epsilon + 0.5)/1);
        if(timestep < 0)
           timestep = 0;
        else if(timestep >= 5000)
            timestep = 5000-1;
        return _namespace_timedarray_1_values[timestep*200 + i];
    }
    double _rand(const int _vectorisation_idx) {
        return rk_double(brian::_mersenne_twister_states[0]);
    }
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

////// HASH DEFINES ///////



void _run_poissongroup_1_thresholder_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numi = 200;
const size_t _numt = 1;
const size_t _numdt = 1;
const size_t _num_spikespace = 201;
    ///// POINTERS ////////////
        
    int32_t* __restrict  _ptr_array_poissongroup_1_i = _array_poissongroup_1_i;
    double*   _ptr_array_defaultclock_t = _array_defaultclock_t;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    int32_t* __restrict  _ptr_array_poissongroup_1__spikespace = _array_poissongroup_1__spikespace;
    _namespace_timedarray_1_values = _timedarray_1_values;



    //// MAIN CODE ////////////
    // scalar code
    const int _vectorisation_idx = -1;
        
    const double t = _ptr_array_defaultclock_t[0];
    const double dt = _ptr_array_defaultclock_dt[0];



    long _count = 0;
    for(size_t _idx=0; _idx<200; _idx++)
    {
        const size_t _vectorisation_idx = _idx;
                
        const int32_t i = _ptr_array_poissongroup_1_i[_idx];
        const double rates = _timedarray_1(t, i);
        const char _cond = _rand(_vectorisation_idx) < (dt * rates);

        if(_cond) {
            _ptr_array_poissongroup_1__spikespace[_count++] = _idx;
        }
    }
    _ptr_array_poissongroup_1__spikespace[200] = _count;

}


