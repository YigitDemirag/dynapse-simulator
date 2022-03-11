#include "code_objects/statemonitor_2_codeobject.h"
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

////// HASH DEFINES ///////



void _run_statemonitor_2_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _num_clock_t = 1;
const size_t _num_indices = 1;
const size_t _num_source_Isoma_mem = 256;
double* const _array_statemonitor_2_t = _dynamic_array_statemonitor_2_t.empty()? 0 : &_dynamic_array_statemonitor_2_t[0];
const size_t _numt = _dynamic_array_statemonitor_2_t.size();
const size_t _numN = 1;
const size_t _numnot_refractory = 256;
const size_t _num_sub_idx = 1;
    ///// POINTERS ////////////
        
    double*   _ptr_array_defaultclock_t = _array_defaultclock_t;
    int32_t*   _ptr_array_statemonitor_2__indices = _array_statemonitor_2__indices;
    double* __restrict  _ptr_array_Core_1_Isoma_mem = _array_Core_1_Isoma_mem;
    double* __restrict  _ptr_array_statemonitor_2_t = _array_statemonitor_2_t;
    int32_t*   _ptr_array_statemonitor_2_N = _array_statemonitor_2_N;
    char* __restrict  _ptr_array_Core_1_not_refractory = _array_Core_1_not_refractory;
    int32_t*   _ptr_array_Core_1_subgroup__sub_idx = _array_Core_1_subgroup__sub_idx;


    _dynamic_array_statemonitor_2_t.push_back(_ptr_array_defaultclock_t[0]);

    const size_t _new_size = _dynamic_array_statemonitor_2_t.size();
    // Resize the dynamic arrays
    _dynamic_array_statemonitor_2_Isoma_mem.resize(_new_size, _num_indices);

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    
    for (int _i = 0; _i < (int)_num_indices; _i++)
    {
        // vector code
        const size_t _idx = _ptr_array_statemonitor_2__indices[_i];
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _sub_idx = _ptr_array_Core_1_subgroup__sub_idx[_idx];
        const double _source_Isoma_mem = _ptr_array_Core_1_Isoma_mem[_sub_idx];
        const double _to_record_Isoma_mem = _source_Isoma_mem;


        _dynamic_array_statemonitor_2_Isoma_mem(_new_size-1, _i) = _to_record_Isoma_mem;
    }

    _ptr_array_statemonitor_2_N[0] = _new_size;


}


