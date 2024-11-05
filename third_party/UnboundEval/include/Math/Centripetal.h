//
// Unbound - Centripetal.h
// Copyright (c) 2017 Unbound Technologies, Inc. All rights reserved.
// 

#pragma once
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

namespace Math
{
	struct CubicPoly
    {
		float c0,c1,c2,c3;
		CubicPoly(){}
		CubicPoly( float x0, float x1, float t0, float t1 )
        {
			c0 = x0;
			c1 = t0;
			c2 = -3.0f*x0 + 3.0f*x1 - 2.0f*t0 - t1;
			c3 =  2.0f*x0 - 2.0f*x1 +      t0 + t1;
		}

		float eval( float t ) const
        {
			float t2 = t  * t;
			float t3 = t2 * t;
			return c0 + c1*t + c2*t2 + c3*t3;
		}
	};

    CubicPoly makeNonUniformCatmullRom(float x0, float x1, float x2, float x3, float dt0, float dt1, float dt2);

	struct Centripedal3
    {
		CubicPoly x;
		CubicPoly y;
		CubicPoly z;

        Centripedal3() = default;
        Centripedal3(const Centripedal3&) = default;
        Centripedal3(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
        glm::vec3 eval(float t) const;
        float getLength() const;
	};
};
