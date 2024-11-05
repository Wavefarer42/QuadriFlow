//
// Unbound Alpha - IA.h
// Copyright (c) 2015 Unbound Technologies, Inc. All rights reserved.
// 

#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/vec1.hpp>
#include <array>
#include "MathAliases.h"

namespace glm
{
	template <typename T>
	struct trange
	{
		T lb;
		T ub;

		
			trange() {}

		
			trange(__thread__ const T& lb, __thread__ const T& ub) :lb(lb), ub(ub) {}

		
			trange(__thread__ const T& v) : lb(v), ub(v) {}

		
			trange(__thread__ const trange<T>& r) : lb(r.lb), ub(r.ub) {}
	};

	template <typename T>
	struct trange3
	{
		trange<T> x, y, z;

		
			trange3() {}

		
			trange3(__thread__ const trange<T>& x, __thread__ const trange<T>& y, __thread__ const trange<T>& z) :
			x(x),
			y(y),
			z(z) {}

		
			trange3(const float lb, const float ub) :
			x(trange<T>(lb, ub)),
			y(trange<T>(lb, ub)),
			z(trange<T>(lb, ub)) {}

		
			trange3(__thread__ const tvec3<T>& lb, __thread__ const tvec3<T>& ub) :
			x(trange<T>(lb.x, ub.x)),
			y(trange<T>(lb.y, ub.y)),
			z(trange<T>(lb.z, ub.z)) {}

		
			trange3(__thread__ const tvec3<T>& v) :
			x(trange<T>(v.x)),
			y(trange<T>(v.y)),
			z(trange<T>(v.z)) {}

		
			trange3(__thread__ const trange3<T>& r) :
			x(r.x),
			y(r.y),
			z(r.z) {}
	};

	template <typename T>
	struct trange2
	{
		trange<T> x, y;

		
			trange2() {}

		
			trange2(__thread__ const trange<T>& x, __thread__ const trange<T>& y) :
			x(x),
			y(y){}

		
			trange2(__thread__ const tvec2<T>& lb, __thread__ const tvec2<T>& ub) :
			x(trange<T>(lb.x, ub.x)),
			y(trange<T>(lb.y, ub.y)) {}

		
			trange2(__thread__ const tvec2<T>& v) :
			x(trange<T>(v.x)),
			y(trange<T>(v.y)) {}

		
			trange2(__thread__ const trange2<T>& r) :
			x(r.x),
			y(r.y) {}

        
            trange2(__thread__ const trange3<T>& r) :
            x(r.x),
            y(r.y) {}

	};

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// functions
	//

	template<typename T>
	
	bool contains(__thread__  trange<T> const& r, __thread__ T const& v)
	{
        return ((v >= r.lb) && (v <= r.ub));
	}

	template<typename T>
	
	bool contains(__thread__  trange<T> const& a, __thread__  trange<T> const& b) // this should be called 'overlaps'
	{ 
		return !((a.ub < b.lb) || (a.lb > b.ub));
	}

    template<typename T>
    
    bool contains(__thread__ __thread__ trange3<T> const& a, __thread__ __thread__ trange3<T> const& b)
    {
        // true, iff b fully contained by a
        return (
            a.x.lb <= b.x.lb &&
            a.y.lb <= b.y.lb &&
            a.z.lb <= b.z.lb &&
            a.x.ub >= b.x.ub &&
            a.y.ub >= b.y.ub &&
            a.z.ub >= b.z.ub
        );
    }

    template<typename T>
    
    bool surfaceArea(__thread__ __thread__ trange3<T> const& a)
    {
        const float x = a.x.ub - a.x.lb;
        const float y = a.y.ub - a.y.lb;
        const float z = a.z.ub - a.z.lb;

        return 2.0f * (x * y + x * z + y * z);
    }


	template<typename T>
	
	trange<T> min(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return trange<T>(FMIN(a.lb, b.lb), FMIN(a.ub, b.ub));
	}

    template<typename T>
    
    trange2<T> min(__thread__ trange2<T> const& a, __thread__ trange2<T> const& b)
    {
        return trange2<T>(min(a.x, b.x), min(a.y, b.y));
    }

	template<typename T>
	
	trange<T> min(__thread__ trange<T> const& a, T b)
	{
		return trange<T>(FMIN(a.lb, b), FMIN(a.ub, b));
	}

	template<typename T>
	
	trange<T> min(T b, __thread__ trange<T> const& a)
	{
		return trange<T>(FMIN(a.lb, b), FMIN(a.ub, b));
	}

	template<typename T>
	
	trange<T> max(__thread__  trange<T> const& a, __thread__  trange<T> const& b)
	{
		return trange<T>(FMAX(a.lb, b.lb), FMAX(a.ub, b.ub));
	}

	template<typename T>
	
	trange<T> max(__thread__  trange<T> const& a, T b)
	{
		return trange<T>(FMAX(a.lb, b), FMAX(a.ub, b));
	}

	template<typename T>
	
	trange<T> max(T b, __thread__  trange<T> const& a)
	{
		return trange<T>(FMAX(a.lb, b), FMAX(a.ub, b));
	}


    template<typename T>
    trange<T> extend(__thread__  trange<T> const& a, float e)
    {
        return trange<T>(a.lb-e, a.ub+e);
    }

    template <typename T>
    trange3<T> extend(__thread__ trange3<T> const& a, float e)
    {
        return trange3<T>(
            extend(a.x, e),
            extend(a.y, e),
            extend(a.z, e)
        );
    }

    template<typename T>
    trange<T> scale(__thread__  trange<T> const& a, float e)
    {
        return trange<T>(a.lb*e, a.ub*e);
    }

    template <typename T>
    trange3<T> scale(__thread__ trange3<T> const& a, float e)
    {
        return trange3<T>(
            scale(a.x, e),
            scale(a.y, e),
            scale(a.z, e)
        );
    }


    template<typename T>
    trange<T> merge(__thread__  trange<T> const& a, __thread__  trange<T> const& b)
    {
        return trange<T>(FMIN(a.lb, b.lb), FMAX(a.ub, b.ub));
    }

    template <typename T>
    trange3<T> merge(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
    {
        return trange3<T>(
            merge(a.x, b.x),
            merge(a.y, b.y),
            merge(a.z, b.z)
        );
    }

    template<typename T>
    trange<T> merge(__thread__  trange<T> const& a, T b)
    {
        return trange<T>(FMIN(a.lb, b), FMAX(a.ub, b));
    }

    template <typename T>
    trange3<T> merge(__thread__ trange3<T> const& a, __thread__ tvec3<T> const& b)
    {
        return trange3<T>(
            merge(a.x, b.x),
            merge(a.y, b.y),
            merge(a.z, b.z)
        );
    }

    
    template<typename T>
    bool empty(__thread__  trange<T> const& a)
    {
        return (a.lb >= a.ub);
    }

    template <typename T>
    
    bool empty(__thread__ trange3<T> const& a)
    {
        return (empty(a.x) || empty(a.y) || empty(a.z));
    }

    template <typename T>
    
    bool intersects(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
    {
        return (a.x.lb <= b.x.ub && a.x.ub >= b.x.lb) &&
               (a.y.lb <= b.y.ub && a.y.ub >= b.y.lb) &&
               (a.z.lb <= b.z.ub && a.z.ub >= b.z.lb);
    }

    template <typename T>
    
    trange<T> intersect(__thread__ const trange<T>& a, __thread__ const trange<T>& b)
    {
        return trange<T>(glm::max(a.lb, b.lb), glm::min(a.ub, b.ub));
    }

    template <typename T>
    
    trange3<T> intersect(__thread__ const trange3<T>& a, __thread__ const trange3<T>& b)
    {
        return trange3<T>(
            intersect(a.x, b.x),
            intersect(a.y, b.y),
            intersect(a.z, b.z)
        );
    }

    template <typename T>
    trange<T> sign(__thread__ trange<T> const& a)
    {
        return trange<T>(sign(a.lb), sign(a.ub));
    }


	/////////////////////
	// Addition
	/////////////////////

	template <typename T>
	
	trange<T> operator+(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return trange<T>(a.lb + b.lb, a.ub + b.ub);
	}

	template <typename T>
	
	trange<T> operator+(__thread__ trange<T> const& a, __thread__ T const& b)
	{
		return trange<T>(a.lb + b, a.ub + b);
	}

	template <typename T>
	
	trange<T> operator+(__thread__ T const& b, __thread__ trange<T> const& a)
	{
		return trange<T>(a.lb + b, a.ub + b);
	}

	template <typename T>
	
	trange3<T> operator+(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	template <typename T>
	
	trange3<T> operator+(__thread__ tvec3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	template <typename T>
	
	trange3<T> operator+(__thread__ trange3<T> const& b, __thread__ tvec3<T> const& a)
	{
		return trange3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	/////////////////////
	// Subtraction
	/////////////////////
	template <typename T>
	
	trange<T> operator-(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return trange<T>(a.lb - b.ub, a.ub - b.lb);
	}

	template <typename T>
	
	trange<T> operator-(__thread__ trange<T> const& a, __thread__ T const& b)
	{
		return trange<T>(a.lb - b, a.ub - b);
	}


	template <typename T>
	
	trange<T> operator-(__thread__ T const& b, __thread__ trange<T> const& a)
	{
		return trange<T>(b - a.ub, b - a.lb);
	}

	template <typename T>
	
	trange3<T> operator-(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
	}

    template <typename T>
	
	trange2<T> operator-(__thread__ trange2<T> const& a, __thread__ trange2<T> const& b)
	{
		return trange2<T>(a.x - b.x, a.y - b.y);
	}

	template <typename T>
	
	trange3<T> operator-(__thread__ tvec3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	template <typename T>
	
	trange3<T> operator-(__thread__ trange3<T> const& a, __thread__ tvec3<T> const& b)
	{
		return trange3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	/////////////////////
	// Multiplication
	/////////////////////
	template <typename T>
	
	trange<T> operator*(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		const T f0(a.lb * b.lb);
		const T f1(a.lb * b.ub);
		const T f2(a.ub * b.lb);
		const T f3(a.ub * b.ub);
		return trange<T>(
			FMIN(FMIN(f0, f1), FMIN(f2, f3)),
			FMAX(FMAX(f0, f1), FMAX(f2, f3))
		);
	}

	template <typename T>
	
	trange<T> operator*(__thread__ trange<T> const& a, __thread__ T const& c)
	{
		const T f0(a.lb * c);
		const T f2(a.ub * c);

		return trange<T>(
			FMIN(f0, f2),
			FMAX(f0, f2)
		);
	}

	template <typename T>
	
	trange<T> operator*(__thread__ T const& c, __thread__ trange<T> const& a)
	{
		const T f0(a.lb * c);
		const T f2(a.ub * c);

		return trange<T>(
			FMIN(f0, f2),
			FMAX(f0, f2)
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			a.x * b.x,
			a.y * b.y,
			a.z * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ T const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			a * b.x,
			a * b.y,
			a * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ trange3<T> const& b, __thread__ T const& a)
	{
		return trange3<T>(
			a * b.x,
			a * b.y,
			a * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ tvec3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			a.x * b.x,
			a.y * b.y,
			a.z * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ trange3<T> const& b, __thread__ tvec3<T> const& a)
	{
		return trange3<T>(
			a.x * b.x,
			a.y * b.y,
			a.z * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ trange<T> const& a, __thread__ tvec3<T> const& b)
	{
		return trange3<T>(
			a * b.x,
			a * b.y,
			a * b.z
		);
	}

	template <typename T>
	
	trange3<T> operator*(__thread__ tvec3<T> const& b, __thread__ trange<T> const& a)
	{
		return trange3<T>(
			a * b.x,
			a * b.y,
			a * b.z
		);
	}

    // range2 multiplication
    template <typename T>
	
	trange2<T> operator*(__thread__ trange2<T> const& a, __thread__ trange2<T> const& b)
	{
		return trange2<T>(
			a.x * b.x,
			a.y * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ T const& a, __thread__ trange2<T> const& b)
	{
		return trange2<T>(
			a * b.x,
			a * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ trange2<T> const& b, __thread__ T const& a)
	{
		return trange2<T>(
			a * b.x,
			a * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ tvec2<T> const& a, __thread__ trange2<T> const& b)
	{
		return trange2<T>(
			a.x * b.x,
			a.y * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ trange2<T> const& b, __thread__ tvec2<T> const& a)
	{
		return trange2<T>(
			a.x * b.x,
			a.y * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ trange<T> const& a, __thread__ tvec2<T> const& b)
	{
		return trange2<T>(
			a * b.x,
			a * b.y
		);
	}

	template <typename T>
	
	trange2<T> operator*(__thread__ tvec2<T> const& b, __thread__ trange<T> const& a)
	{
		return trange2<T>(
			a * b.x,
			a * b.y
		);
	}

    template <typename T>
	
	trange2<T> operator*(__thread__ trange2<T> const& b, __thread__ trange<T> const& a)
	{
		return trange2<T>(
			a * b.x,
			a * b.y
		);
	}

    template <typename T>
	
	trange2<T> operator*(__thread__ trange<T> const& b, __thread__ trange2<T> const& a)
	{
		return trange2<T>(
			a.x * b,
			a.y * b
		);
	}

	/////////////////////
	// Division
	/////////////////////
	template <typename T>
	
	trange<T> operator/(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		const T f0(a.lb / b.lb);
		const T f1(a.lb / b.ub);
		const T f2(a.ub / b.lb);
		const T f3(a.ub / b.ub);
		return trange<T>(
			FMIN(FMIN(f0, f1), FMIN(f2, f3)),
			FMAX(FMAX(f0, f1), FMAX(f2, f3))
		);
	}

    template <typename T>
    
    trange<T> operator/(__thread__ T const& a, __thread__ trange<T> const& b)
    {
        const T f0(a / b.lb);
        const T f1(a / b.ub);
        return trange<T>(
            FMIN(f0, f1),
            FMAX(f0, f1)
        );
    }

	template <typename T>
	
	trange3<T> operator/(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			a.x / b.x,
			a.y / b.y,
			a.z / b.z
		);
	}

    template <typename T>
    
    trange3<T> operator/(__thread__ T const& a, __thread__ trange3<T> const& b)
    {
        return trange3<T>(
            a / b.x,
            a / b.y,
            a / b.z
        );
    }

    
	/////////////////////
	// Transcententals
	/////////////////////
	template <typename T>
	
	trange<T> sqrt(__thread__ trange<T> const &a)
	{
		return trange<T>(FSQRT(a.lb), FSQRT(a.ub));
	}

	/////////////////////
	// Linear Algebra
	/////////////////////
	template <typename T>
	 inline
	trange<T> pow2(__thread__ trange<T> const& a)
	{
		const trange<T> a2(a.lb*a.lb, a.ub*a.ub);
		return (a.lb >= 0.0f) ? a2 : (a.ub < 0.0f) ? trange<T>(a2.ub, a2.lb) : trange<T>(0.0f, FMAX(a2.lb, a2.ub));
	}

	template <typename T>
	 inline
	trange3<T> pow2(__thread__ trange3<T> const& a)
	{
		const trange<T> x(pow2(a.x));
		const trange<T> y(pow2(a.y));
		const trange<T> z(pow2(a.z));
		return trange3<T>(x, y, z);
	}

	template <typename T>
	 inline
	trange2<T> pow2(__thread__ trange2<T> const& a)
	{
		const trange<T> x(pow2(a.x));
		const trange<T> y(pow2(a.y));
		return trange2<T>(x, y);
	}

	template <typename T>
	 inline
	trange<T> length2(__thread__ trange3<T> const& a)
	{
		const trange3<T> c(pow2(a));
		return c.x + c.y + c.z;
	}

	template <typename T>
	 inline
	trange<T> length(__thread__ trange3<T> const& a)
	{
		const trange<T> l2(length2(a));
		return sqrt(l2);
	}

	template <typename T>
	 inline
	trange<T> length2(__thread__ trange2<T> const& a)
	{
		const trange2<T> c(pow2(a));
		return c.x + c.y;
	}

	template <typename T>
	 
	trange<T> length(__thread__ trange2<T> const& a)
	{
		const trange<T> l2(length2(a));
		return sqrt(l2);
	}

	template <typename T>
	
	trange<T> dot(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		const trange3<T> c(a*b);
		return c.x + c.y + c.z;
	}

    template <typename T>
	
	trange<T> dot(__thread__ trange2<T> const& a, __thread__ tvec2<T> const& b)
	{
		const trange2<T> c(a*b);
		return c.x + c.y;
	}

    template <typename T>
    
    trange<T> dot(__thread__ trange2<T> const& a, __thread__ trange2<T> const& b)
    {
        const trange2<T> c(a*b);
        return c.x + c.y;
    }

	template <typename T>
	
	trange3<T> cross(__thread__ trange3<T> const& a, __thread__ tvec3<T> const& b)
	{
		return trange3<T>(
			b.z * a.y - b.y * a.z,
			b.x * a.z - b.z * a.x,
			b.y * a.x - b.x * a.y
		);
	}

	template <typename T>
	
	trange3<T> cross(__thread__ tvec3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			b.z * a.y - b.y * a.z,
			b.x * a.z - b.z * a.x,
			b.y * a.x - b.x * a.y
		);
	}

	template <typename T>
	
	trange3<T> cross(__thread__ trange3<T> const& a, __thread__ trange3<T> const& b)
	{
		return trange3<T>(
			b.z * a.y - b.y * a.z,
			b.x * a.z - b.z * a.x,
			b.y * a.x - b.x * a.y
		);
	}

	/////////////////////
	// Comparison
	/////////////////////
	template <typename T>
	
	bool operator!=(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return (a.lb != b.lb) || (a.ub != b.ub);
	}

	template <typename T>
	
	bvec2 operator<(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return bvec2((a.lb < b.lb), (a.ub < b.ub));
	}

	template <typename T>
	
	bvec2 operator<=(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return bvec2((a.lb <= b.lb), (a.ub <= b.ub));
	}

	template <typename T>
	
	bvec2 operator>(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return bvec2((a.lb > b.lb), (a.ub > b.ub));
	}

	template <typename T>
	
	bvec2 operator>=(__thread__ trange<T> const& a, __thread__ trange<T> const& b)
	{
		return bvec2((a.lb >= b.lb), (a.ub >= b.ub));
	}

	/////////////////////
	// Other Arithmetic
	/////////////////////
	template <typename T>
	
	trange<T> operator-(__thread__ trange<T> const& a)
	{
		return trange<T>(-a.ub, -a.lb);
	}

	template <typename T>
	
	trange<T> abs(__thread__ trange<T> const& a)
	{
		if (a.lb >= static_cast<T>(0))
			return a;
		else if (a.ub <= static_cast<T>(0))
			return trange<T>(-a.ub, -a.lb);
		return trange<T>(static_cast<T>(0), FMAX(-a.lb, a.ub));
	}

    template <typename T>
    
    trange2<T> abs(__thread__ trange2<T> const& a)
    {
        return trange2<T>(
            abs(a.x),
            abs(a.y)
        );
    }

	template <typename T>
	
	trange3<T> abs(__thread__ trange3<T> const& a)
	{
		return trange3<T>(
			abs(a.x),
			abs(a.y),
			abs(a.z)
		);
	}

	template <typename T>
	
	trange<T> smin(__thread__ trange<T> const& a, __thread__ trange<T> const& b, T r)
	{
        const trange<T> ux(max(r - a, 0.0f));
        const trange<T> uy(max(r - b, 0.0f));
        const trange<T> len(sqrt(pow2(ux) + pow2(uy)));
        return max(r, min(a, b)) - len;

        //const trange<T> e(max(r - abs(a - b), static_cast<T>(0)));
        //return min(a, b) - pow2(e) * (static_cast<T>(0.25) / r);
	}

	template <typename T>
	
	trange<T> smax(__thread__ trange<T> const& a, __thread__ trange<T> const& b, T r)
	{
        const trange<T> ux(max(r + a, 0.0f));
        const trange<T> uy(max(r + b, 0.0f));
        const trange<T> len(sqrt(pow2(ux) + pow2(uy)));
        return min(-r, max(a, b)) + len;

		//const trange<T> e(max(r - abs(a - b), static_cast<T>(0)));
		//return max(a, b) + pow2(e) * (static_cast<T>(0.25) / r);
	}

	/////////////////////
	// Quaternion
	/////////////////////
	template <typename T>
	
	trange3<T> operator* (__thread__ tquat<T, highp> const& q, __thread__ trange3<T> const& v)
	{
		const trange3<T> uv(cross(tvec3<T>(q.x, q.y, q.z), v));
		const trange3<T> uuv(cross(tvec3<T>(q.x, q.y, q.z), uv));
		return v + static_cast<T>(2) * (q.w * uv + uuv);
	}

    template <typename T>
    
    trange<T> mix(__thread__ trange<T> const& x, __thread__ trange<T> const& y, T __thread__ const& a)
    {
        return x + a * (y - x);
    }

	/////////////////////
	// Signed Distance Functions
	/////////////////////
	template <typename T>
	
	trange<T> sdBox(__thread__ trange3<T> const& p, vec3 b)
	{
		const trange3<T> d(abs(p) - b);
		const trange<T> a(min(max(d.x, max(d.y, d.z)), trange<T>(static_cast<T>(0))));
		const trange3<T> d3(
			max(d.x, trange<T>(static_cast<T>(0))),
			max(d.y, trange<T>(static_cast<T>(0))),
			max(d.z, trange<T>(static_cast<T>(0)))
		);
		return length(d3) + a;
	}

	template <typename T>
	
	trange<T> sdSphere(__thread__ trange3<T> const& p, float r)
	{
		const trange<T> l(length(p));
		const trange<T> sd(l - r);
		return sd;
	}

#if 1
	template <typename T>
	
	trange<T> sdSuperprim(__thread__ trange3<T> const& p, __thread__ vec3 const& dim, const float thickness, const float cornerRadiusX, const float cornerRadiusY)
	{
		const float mdxy(max(dim.x, dim.y));
		const float rx(cornerRadiusX*mdxy);
		const float ry(cornerRadiusY*dim.z);

		const trange3<T> d(abs(p) - dim);
		const trange2<T> dr(
			max(d.x + rx, 0.0f),
			max(d.y + rx, 0.0f)
		);
		const float sw(thickness*mdxy);
		const trange<T> q(abs(length(dr) + min(-rx, max(d.x, d.y)) + sw) - sw);
		const trange2<T> qq(max(q + ry, 0.0f), max(d.z + ry, 0.0f));
		return length(qq) + min(-ry, max(q, d.z));
	}
#else
    // incomplete IA uberprim
    template <typename T>
    
        trange<T> sdSuperprim(__thread__ trange3<T> const& p, vec3 const& dim, const float conef, const float cornerRadiusX, const float cornerRadiusY)
    {
        const float mdxy = FMAX(dim.x, dim.y);
        const float rx = cornerRadiusX * mdxy;
        const float ry = cornerRadiusY * dim.z;

        const vec3 dr(
            dim.x - rx,
            dim.y - rx,
            dim.z - ry
        );
        
        const float cone = -(1.0f - conef)*FMIN(dim.x, dim.y);
        vec2 ba = vec2(cone, -2.0*dr.z);
        float sz2 = ba.y;
        ba /= dot(ba, ba);

        const trange3<T> d(
            abs(p.x) - dr.x,
            abs(p.y) - dr.y,
            abs(p.z) - dr.z
        );
        const trange2<T> pd(max(d.x, 0.0f), max(d.y, 0.0f));
        const trange<T> q = length(pd) + min(0.0f, max(d.x, d.y)) - rx + ry;

        trange2<T> pa(q, p.z - dr.z);
        trange2<T> diag = pa - trange2<T>(cone, sz2) * min(max(dot(pa, ba), 0.0f), 1.0f);
        trange2<T> h0(max(q - cone, 0.0f), p.z + dr.z);
        trange2<T> h1(max(q, 0.0f), p.z - dr.z);

        return sqrt(min(dot(diag, diag), min(dot(h0, h0), dot(h1, h1)))) /* * sign(max(dot(pa, trange2<T>(-ba.y, ba.x)), d.z))*/ - ry;
    }
#endif
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // common typedefs
    //
    typedef trange3<int32_t> irange3;
    typedef trange3<uint32_t> urange3;
    typedef trange3<float> range3;
    typedef trange2<float> range2;
    typedef trange2<int32_t> irange2;
    typedef trange2<uint32_t> urange2;
    typedef trange<float> range;
    typedef trange<uint32_t> urange;
    typedef trange<int32_t> irange;
    
    
    
    static inline glm::range2 solve_quadratic(const __thread__ float& a, const __thread__ range& b, const __thread__ range& c)
    {
        using namespace glm;
        // ax^2 + bx + c = 0, a non-zero
        const range q = b*b - 4.0f*a*c;
        if (q.lb < 0.0f || q.ub < 0.0f)
        {
            return range2(1000.0f, 1000.0f);
        }
        
        const float aDoubleDiv = 1.0f/(2.0f * a);
        range r0 = -b * aDoubleDiv;
        range r1 = sqrt(q) * aDoubleDiv;
        return range2(r0 - r1, r0 + r1);
    }

    // Media Molecule's max norm sphere: https://www.shadertoy.com/view/Mt2XWG
    
    inline range sdMaxNormSphere(glm::range3 __thread__ const& p, glm::vec3 __thread__ const& r)
    {
        using namespace glm;
        
        // move ellipse so that target point is at origin, centre in positive space
        // f(v) = (v.x - c.x)^2/r.x^2 + (v.y - c.y)^2/r.y^2 + (v.z - c.z)^2/r.z^2 - 1
        const range3 c = abs(p);
        const range3 c2 = c*c;
        const vec3 r2 = r*r;
        range d(1000.0f);

        const vec3 r2Div(1.0f/r2.x, 1.0f/r2.y, 1.0f/r2.z);

        // gather terms of quadratic
        const vec3 qa(r2Div);
        const range3 qb = -2.0f*c * r2Div;
        const range3 qc = c2 * r2Div;
        const range qcs = (qc.x + qc.y + qc.z) - 1.0f;

        // check corners:
        // solve f(v)=0 for v.x=v.y=v.z=t
        const range2 t0 = abs(solve_quadratic(qa.x+qa.y+qa.z, qb.x+qb.y+qb.z, qcs));
        d = min(min(d, t0.x), t0.y);

        // interior of convex shape always hits corners first, so early out
        if (qcs.ub <= 0.0f)
        {
            return -d;
        }

        // check edges:
        // df/dx=0 => v.x=c.x, solve f(v)=0 for v.x=c.x, v.y=v.z=t
        // then do the same for y and z cases
        range2 t = abs(solve_quadratic(qa.y + qa.z, qb.y + qb.z, qc.y + qc.z - 1.0f));
        d = min(d, max(min(t.x, t.y), c.x));

        t = abs(solve_quadratic(qa.x + qa.z, qb.x + qb.z, qc.x + qc.z - 1.0f));
        d = min(d, max(min(t.x, t.y), c.y));

        t = abs(solve_quadratic(qa.x + qa.y, qb.x + qb.y, qc.x + qc.y - 1.0f));
        d = min(d, max(min(t.x, t.y), c.z));

        // check faces:
        // df/dx=df/dy=0 => v.xy=c.xy, so f(v)=0 => |v.z - c.z|=r.z
        d = min(d, max(max(c.x, c.y), abs(c.z - r.z)));
        d = min(d, max(max(c.x, abs(c.y - r.y)), c.z));
        d = min(d, max(max(abs(c.x - r.x), c.y), c.z));

        // done
        return d;
    }

    
    static range sdEllipsoid(__thread__ glm::range3 const& p, __thread__ glm::vec3 const& dim, float __thread__ const& thickness, float __thread__ const& cornerRadiusY, float __thread__ const& slice)
    {
        using namespace glm;
        const float round = cornerRadiusY - thickness;
                
        const range dPerfectSphere = length(p) - dim.x + round;

        // mediamolecule's max norm sphere, needed for squashed spheroids
        const range dMaxNormSphere = sdMaxNormSphere(p, dim) + round;

        const float maxScale = max(dim.x, max(dim.y, dim.z));
        const float minScale = min(dim.x, min(dim.y, dim.z));
        const float squash = clamp((maxScale - minScale) * 1000.0f, 0.0f, 1.0f);

        const range d = mix(dPerfectSphere, dMaxNormSphere, squash);

        // slicing
        const float sliceY = (-dim.y * ((slice * 2.0f) - 1.0f));
        const range slicedY = -p.y - sliceY;
        const range2 A = range2(max(d, 0.0f), max(slicedY, 0.0f));
        return length(A) + min(max(d, slicedY), 0.0f) - round;
    }
    
    
    //
    // Some functions to do aabb style calculations
    //
    
    template <typename T>
    T center(__thread__ trange<T> const& a)
    {
        return (a.lb + a.ub) / static_cast<T>(2);
    }

    template <typename T>
    tvec3<T> center(__thread__ trange3<T> const& a)
    {
        return tvec3<T>(
            (a.x.lb + a.x.ub) / static_cast<T>(2),
            (a.y.lb + a.y.ub) / static_cast<T>(2),
            (a.z.lb + a.z.ub) / static_cast<T>(2)
        );
    }
    
    template <typename T>
    T width(__thread__ trange<T> const& a)
    {
        return (a.ub - a.lb);
    }

    template <typename T>
    tvec3<T> width(__thread__ trange3<T> const& a)
    {
        return tvec3<T>(
            (a.x.ub - a.x.lb),
            (a.y.ub - a.y.lb),
            (a.z.ub - a.z.lb)
        );
    }

    template <typename T>
    bool inside(__thread__ trange<T> const& a, T p)
    {
        return (p >= a.lb && p <= a.ub);
    }

	template <typename T>
	bool inside(__thread__ trange<T> const& a, __thread__ trange<T> const& p)
	{
		return (p.lb >= a.lb && p.ub <= a.ub);
	}

    template <typename T>
    bool inside(__thread__ trange3<T> const& a, __thread__ const glm::tvec3<T>& p)
    {
        return inside(a.x, p.x) && inside(a.y, p.y) && inside(a.z, p.z);
    }

	template <typename T>
	bool inside(__thread__ trange3<T> const& a, __thread__ const trange3<T>& p)
	{
		return inside(a.x, p.x) && inside(a.y, p.y) && inside(a.z, p.z);
	}

    template <typename T>
    bool inside(__thread__ trange2<T> const& a, __thread__ const glm::tvec2<T>& p)
    {
        return inside(a.x, p.x) && inside(a.y, p.y);
    }


	template <typename T>
	bool overlaps(__thread__ trange3<T> const& a, __thread__ const trange3<T>& b)
	{
		const glm::tvec3<T> d1(
			b.x.lb - a.x.ub,
			b.y.lb - a.y.ub,
			b.z.lb - a.z.ub
		);
		
		const glm::tvec3<T> d2(
			a.x.lb - b.x.ub,
			a.y.lb - b.y.ub,
			a.z.lb - b.z.ub
		);

		if (d1.x > 0.0f || d1.y > 0.0f || d1.z > 0.0f)
			return false;

		if (d2.x > 0.0f || d2.y > 0.0f || d2.z > 0.0f)
			return false;

		return true;
	}

    
    template <typename T>
    std::array<glm::tvec3<T>, 8> corners(__thread__ trange3<T> const& a)
    {
        return {
            tvec3<T>(a.x.lb, a.y.lb, a.z.lb),
            tvec3<T>(a.x.ub, a.y.lb, a.z.lb),
            tvec3<T>(a.x.lb, a.y.ub, a.z.lb),
            tvec3<T>(a.x.ub, a.y.ub, a.z.lb),
            tvec3<T>(a.x.lb, a.y.lb, a.z.ub),
            tvec3<T>(a.x.ub, a.y.lb, a.z.ub),
            tvec3<T>(a.x.lb, a.y.ub, a.z.ub),
            tvec3<T>(a.x.ub, a.y.ub, a.z.ub)
        };
    }


}


