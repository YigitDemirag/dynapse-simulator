#include "code_objects/Core_2_thresholder_1_codeobject.h"
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



void _run_Core_2_thresholder_1_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numIsoma_mem = 256;
const size_t _numlastspike = 256;
const size_t _num_spikespace = 257;
const size_t _numnot_refractory = 256;
const size_t _numIsoma_th = 1;
const size_t _numt = 1;
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_Core_2_Isoma_mem = _array_Core_2_Isoma_mem;
    double* __restrict  _ptr_array_Core_2_lastspike = _array_Core_2_lastspike;
    int32_t* __restrict  _ptr_array_Core_2__spikespace = _array_Core_2__spikespace;
    char* __restrict  _ptr_array_Core_2_not_refractory = _array_Core_2_not_refractory;
    double*   _ptr_array_Core_2_Isoma_th = _array_Core_2_Isoma_th;
    double*   _ptr_array_defaultclock_t = _array_defaultclock_t;



    //// MAIN CODE ////////////
    // scalar code
    const int _vectorisation_idx = -1;
        
    const double Isoma_th = _ptr_array_Core_2_Isoma_th[0];



    long _count = 0;
    for(size_t _idx=0; _idx<256; _idx++)
    {
        const size_t _vectorisation_idx = _idx;
                
        const double Isoma_mem = _ptr_array_Core_2_Isoma_mem[_idx];
        const char not_refractory = _ptr_array_Core_2_not_refractory[_idx];
        char _cond;
        if(!not_refractory)
            _cond = (Isoma_mem > Isoma_th) && false;
        else 
            _cond = (Isoma_mem > Isoma_th) && true;

        if(_cond) {
            _ptr_array_Core_2__spikespace[_count++] = _idx;
            _ptr_array_Core_2_not_refractory[_idx] = false;
            _ptr_array_Core_2_lastspike[_idx] = _ptr_array_defaultclock_t[0];
        }
    }
    _ptr_array_Core_2__spikespace[256] = _count;

}


