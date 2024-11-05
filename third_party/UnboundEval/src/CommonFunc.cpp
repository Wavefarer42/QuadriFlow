//
// Unbound - CommonFunc.cpp
// Copyright (c) 2014-current Unbound Technologies, Inc. All rights reserved.
// 
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/vec1.hpp>
#include <string>
#include <algorithm>
#include <cctype>

#include "Util/CommonFunc.h"
#include "CSG/IA.h"

namespace Util
{
	glm::ivec3 saveNextMultipleOfN(const glm::ivec3& a, const int n)
	{
		glm::ivec3 aa = a;
		if (aa.x <= 0){ aa.x = n; }
		if (aa.y <= 0){ aa.y = n; }
		if (aa.z <= 0){ aa.z = n; }
		return glm::ceil(glm::vec3(aa) / float(n)) * float(n);
	}

	int saveNextMultipleOfN(const int a, const int n)
	{
		if (a == 0){
			return n;
		}
		return int(glm::ceil(float(a) / float(n)) * float(n));
	}

	glm::ivec3 nextMultipleOfN(const glm::ivec3& a, const int n)
	{
		return glm::ceil(glm::vec3(a) / float(n)) * float(n);
	}

	glm::ivec3 nextMultipleOfN(const glm::ivec3& a, const glm::ivec3& n)
	{
		return glm::ceil(glm::vec3(a) / glm::vec3(n)) * glm::vec3(n);
	}

	glm::ivec2 nextMultipleOfN(const glm::ivec2& a, const int n)
	{
		return glm::ceil(glm::vec2(a) / float(n)) * float(n);
	}

	int nextMultipleOfN(const int a, const int n)
	{
		return glm::ceil(float(a) / float(n)) * float(n);
	}

	unsigned int roundUpToNextPowerOfTwo(unsigned int x)
	{
		x--;
		x |= x >> 1;  // handle  2 bit numbers
		x |= x >> 2;  // handle  4 bit numbers
		x |= x >> 4;  // handle  8 bit numbers
		x |= x >> 8;  // handle 16 bit numbers
		x |= x >> 16; // handle 32 bit numbers
		x++;
		return x;
	}

    size_t roundUpToNextPowerOfTwo(size_t x)
    {
        assert(x > 0);
        x--;
        x |= x >> 1;  // handle  2 bit numbers
        x |= x >> 2;  // handle  4 bit numbers
        x |= x >> 4;  // handle  8 bit numbers
        x |= x >> 8;  // handle 16 bit numbers
        x |= x >> 16; // handle 32 bit numbers
        x |= x >> 32; // handle 64 bit numbers
        x++;
        return x;
    }


	glm::vec3 hsv2rgb(float h, float s, float v)
	{
		glm::vec4 K = glm::vec4(1.0f, 2.0f / 3.0f, 1.0f / 3.0f, 3.0f);
		glm::vec3 p = glm::abs(glm::fract(glm::vec3(h) + glm::vec3(K)) * 6.0f - glm::vec3(K.w));
		return v * glm::mix(glm::vec3(K.x), glm::clamp(p - glm::vec3(K.x), 0.0f, 1.0f), s);
	}

	glm::vec3 rgb2hsv(float r, float g, float b)
	{
		using namespace glm;
		const vec4 K(0.0f, -1.0f / 3.0f, 2.0f / 3.0f, -1.0f);
		vec4 p = mix(vec4(b, g, K.w, K.z), vec4(g, b, K.x, K.y), vec4(step(vec1(b), vec1(g)).x));
		vec4 q = mix(vec4(p.x, p.y, p.w, r), vec4(r, p.y, p.z, p.x), vec4(step(vec1(p.x), vec1(r)).x));

		float d = q.x - min(q.w, q.y);
		float e = 1.0e-10f;
		return vec3(abs(q.z + (q.w - q.y) / (6.0f * d + e)), d / (q.x + e), q.x);
	}

	glm::vec3 spherical(float theta, float phi)
	{
		return glm::vec3(glm::cos(phi) * glm::sin(theta), glm::sin(phi) * glm::sin(theta), glm::cos(theta));
	}

	float angularDifference(const glm::quat& a, const glm::quat& b)
	{
		glm::quat rot = a * glm::conjugate(b);
		return 2.0f * glm::acos(rot.w);
	}

	glm::quat lookAtOrientation(glm::vec3 dirNorm, glm::vec3 desiredUp)
	{
		using namespace glm;

        if (glm::abs(glm::abs(glm::dot(dirNorm, desiredUp)) - 1.0f) < FLT_EPSILON)
        {
            desiredUp = glm::vec3(dirNorm.y, -dirNorm.x, dirNorm.z);
        }
		const vec3 s(normalize(cross(dirNorm, desiredUp)));
		const vec3 u(cross(s, dirNorm));

		mat3 m(1);
		m[0][0] = s.x;
		m[1][0] = s.y;
		m[2][0] = s.z;
		m[0][1] = u.x;
		m[1][1] = u.y;
		m[2][1] = u.z;
		m[0][2] = -dirNorm.x;
		m[1][2] = -dirNorm.y;
		m[2][2] = -dirNorm.z;
		return conjugate(normalize(quat_cast(m)));
	}

	bool startsWith(const std::string& str, const std::string& prefix)
	{
		return (prefix.length() <= str.length()) && std::equal(prefix.begin(), prefix.end(), str.begin());
	}

	bool strContains( const std::string& str, const std::string& to_find )
	{
        std::string s1( str );
        std::transform(s1.begin(), s1.end(), s1.begin(), ::tolower);
        if( s1.find(to_find) != std::string::npos ){ return true; }
		return false;
	}

	// trim from start (in place)
	void ltrim(std::string &s) 
	{
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
			return !std::isspace(ch);
		}));
	}

	// trim from end (in place)
	void rtrim(std::string &s) 
	{
		s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
			return !std::isspace(ch);
		}).base(), s.end());
	}

	// trim from both ends (in place)
	void trim(std::string &s) 
	{
		ltrim(s);
		rtrim(s);
	}

	glm::vec2 intersectRayAABB(const glm::vec3& p, const glm::vec3& d, const glm::vec3& lb, const glm::vec3& ub)
	{
		using namespace glm;

		bvec3 parallel = lessThan(abs(d), vec3(FLT_EPSILON));

		vec3 d_safe = vec3(
			parallel.x ? 1.0f : d.x,
			parallel.y ? 1.0f : d.y,
			parallel.z ? 1.0f : d.z
		);

		const glm::vec3 inv_d = 1.0f / d_safe;
		vec3 t1 = (lb - p) * inv_d;
		vec3 t2 = (ub - p) * inv_d;
		vec3 vmin(
			parallel.x ? -FLT_MAX : min(t1.x, t2.x),
			parallel.y ? -FLT_MAX : min(t1.y, t2.y),
			parallel.z ? -FLT_MAX : min(t1.z, t2.z)
		);
		vec3 vmax(
			parallel.x ? FLT_MAX : max(t1.x, t2.x),
			parallel.y ? FLT_MAX : max(t1.y, t2.y),
			parallel.z ? FLT_MAX : max(t1.z, t2.z)
		);

		float tmin = max(max(vmin.x, vmin.y), vmin.z);
		float tmax = min(min(vmax.x, vmax.y), vmax.z);

		if (parallel.x && !(p.x >= lb.x && p.x <= ub.x)) return vec2(FLT_MAX, -FLT_MAX);
		else if (parallel.y && !(p.y >= lb.y && p.y <= ub.y)) return vec2(FLT_MAX, -FLT_MAX);
		else if (parallel.z && !(p.z >= lb.z && p.z <= ub.z)) return vec2(FLT_MAX, -FLT_MAX);
		else return vec2(tmin, tmax);
	}

    glm::vec2 intersectRayAABB(const glm::vec3& p, const glm::vec3& d, const glm::range3& aabb)
    {
        using namespace glm;

        const glm::vec3 inv_d = 1.0f / d;
        float t1 = (aabb.x.lb - p.x)*inv_d.x;
        float t2 = (aabb.x.ub - p.x)*inv_d.x;
        float t3 = (aabb.y.lb - p.y)*inv_d.y;
        float t4 = (aabb.y.ub - p.y)*inv_d.y;
        float t5 = (aabb.z.lb - p.z)*inv_d.z;
        float t6 = (aabb.z.ub - p.z)*inv_d.z;

        float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
        float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
        return vec2(tmin, tmax);
    }

    
    glm::vec2 intersectRayAABBInvDist(const glm::vec3& p, const glm::vec3& inv_d, const glm::vec3& lb, const glm::vec3& ub)
    {
        using namespace glm;

        float t1 = (lb.x - p.x)*inv_d.x;
        float t2 = (ub.x - p.x)*inv_d.x;
        float t3 = (lb.y - p.y)*inv_d.y;
        float t4 = (ub.y - p.y)*inv_d.y;
        float t5 = (lb.z - p.z)*inv_d.z;
        float t6 = (ub.z - p.z)*inv_d.z;

        float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
        float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

        return vec2(tmin, tmax);
    }

    float distanceToLine( const glm::vec3& p, const glm::vec3& a, const glm::vec3& b ) 
	{
        glm::vec3 ab = b - a;
        float t = glm::clamp(glm::dot(p - a, ab) / glm::dot(ab, ab), 0.0f, 1.0f );
        return glm::length((ab*t + a) - p);
    }    

	
/*
	glm::mat2 inverse(glm::mat2 A) {
		return glm::mat2(A[1][1], -A[0][1], -A[1][0], A[0][0]) / glm::determinant(A);
	}

	float inverseSF(glm::vec3 p, float n) {
		float phi = glm::min(atan2(p.y, p.x), PI), cosTheta = p.z;

		float k = max(2, floor(
			log(n * PI * sqrt(5) * (1 - cosTheta*cosTheta))
				/ log(PHI*PHI)));

		float Fk = pow(PHI, k)/sqrt(5);
		float F0 = round(Fk), F1 = round(Fk * PHI);

		float2x2 B = float2x2(
			2*PI*madfrac(F0+1, PHI-1) - 2*PI*(PHI-1),
			2*PI*madfrac(F1+1, PHI-1) - 2*PI*(PHI-1),
			-2*F0/n,
			-2*F1/n);
		float2x2 invB = inverse(B);
		float2 c = floor(
			mul(invB, float2(phi, cosTheta - (1-1/n))));

		float d = INFINITY, j = 0;
		for (uint s = 0; s < 4; ++s) {
			float cosTheta =
				dot(B[1], float2(s%2, s/2) + c) + (1-1/n);
			cosTheta = clamp(cosTheta, -1, +1)*2 - cosTheta;

			float i = floor(n*0.5 - cosTheta*n*0.5);
			float phi = 2*PI*madfrac(i, PHI-1);
			cosTheta = 1 - (2*i + 1)*rcp(n);
			float sinTheta = sqrt(1 - cosTheta*cosTheta);

			float3 q = float3(
				cos(phi)*sinTheta,
				sin(phi)*sinTheta,
				cosTheta);

			float squaredDistance = dot(q-p, q-p);
			if (squaredDistance < d) {
				d = squaredDistance;
				j = i;
			}
		}

		return j;
	}
*/
	
	float mad( float A, float B, float C )
	{
		return (A*B+C);
	}
		

	float madfrac( const float &A, const float &B )
	{
		return mad((A), (B), -floorf((A)*(B)));
	}

	glm::vec3 sphericalFibonacci(const float &i, const float &n) 
	{
		const float TAU_HALF = 3.14159265358979323846f;
		const float PHI = (sqrtf(5.0f)*0.5f + 0.5f);
		
		#define madfrac(A,B) mad((A), (B), -floorf((A)*(B)))
		
		float phi = 2.0f*TAU_HALF*madfrac(i, PHI-1.0f);
		float cosTheta = 1.0f - (2.0f*i + 1.0f)*(1.0f/n);
		float sinTheta = sqrtf(glm::clamp((1.0f - cosTheta*cosTheta), 0.0f, 1.0f));
		return glm::vec3(
			cosf(phi)*sinTheta,
			sinf(phi)*sinTheta,
			cosTheta);
	}

	
	unsigned int hash(unsigned int x)
	{
		x = x * 1024 + x, x ^= x >> 6, x = x * 8 + x, x ^= x >> 11, x = x * 32768 + x;
		return x;
	}

	float rnd(unsigned int& rndSeed)
	{
		rndSeed = hash(rndSeed);
		unsigned int f = ((rndSeed & 0x007fffff) | 0x3f800000);
		return reinterpret_cast<float&>(f) - 1.0f;
	}

	glm::vec3 getRndColor(unsigned int k)
	{
		glm::uint rndSeed = k == 0 ? 0xffffffff : k;
		float R = rnd(rndSeed);
		float G = rnd(rndSeed);
		float B = rnd(rndSeed);
		return glm::vec3(R, G, B);
	}


	

	
}
