#pragma once

//#ifdef _WIN32
//    #include <builtin_types.h>
//#else
//    #define __device__
//    #define __host__
//    #define __forceinline__
//    #define __thread__
//    #define __constant__
//#endif

#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

#ifdef __CUDACC__
    #define FSQRT	sqrtf
    #define FLENGTH length
    #define FMAX	fmaxf
    #define FMIN	fminf
#elif defined(__METAL_VERSION__)
    #define __thread__ thread
    #define FSQRT   metal::sqrt
    #define FLENGTH glm::length
    #define FMAX    metal::max
    #define FMIN    metal::min
#elif TARGET_CPU_ARM64

    #define __thread__ 
    #define FSQRT   glm::sqrt
    #define FLENGTH glm::length
    #define FMAX    glm::max
    #define FMIN    glm::min

#else

    #include <xmmintrin.h>
    #include <glm/glm.hpp>

    #ifndef __attribute_const__
        #if defined(__clang__)
            #define __attribute_const__  __attribute__((__const__))
        #elif defined(_MSC_VER)
            #define __attribute_const__  __declspec(noalias)
        #else
            #define __attribute_const__
        #endif
    #endif

    __attribute_const__ static inline float sqrtss(float number)
    {
        float temp;
        __m128 reg1 = _mm_load_ss(&number);
        __m128 reg2 = _mm_sqrt_ss(reg1);
        _mm_store_ss(&temp, reg2);
        return temp;
    }

    template <typename V>
    static inline float _length(const V& v)
    {
        return sqrtss(glm::dot(v, v));
    }

    #define FSQRT	sqrtss
    #define FLENGTH _length
    #define FMAX	glm::max
    #define FMIN	glm::min
#endif
