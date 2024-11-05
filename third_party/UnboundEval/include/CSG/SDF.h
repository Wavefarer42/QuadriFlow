////
// Unbound Alpha - sdf.cu.h
// Copyright (c) 2015 Unbound Technologies, Inc. All rights reserved.
// 
#pragma once

#include <glm/glm.hpp>
#include "MathAliases.h"

#if 0
static float smin(const float a, const float b, float r)
{
	using namespace glm;
	const float e(FMAX(r - glm::abs(a - b), 0.0f));
	return FMIN(a, b) - e * e * (0.25f / r);
}

static float smax(const float a, const float b, float r)
{
	using namespace glm;
    const float e(FMAX(r - glm::abs(a - b), 0.0f));
    return FMAX(a, b) + e * e * (0.25f / r);
}

#else

static float smin(float a, float b, float r)
{
    using namespace glm;
    float ux = FMAX(r - a, 0.0f);
    float uy = FMAX(r - b, 0.0f);
    float len = sqrt(ux*ux + uy*uy);
    float res =  FMAX(r, FMIN(a, b)) - len;
    return res;
}
static float smax(float a, float b, float r)
{
    using namespace glm;
    float ux = FMAX(r + a, 0.0f);
    float uy = FMAX(r + b, 0.0f);
    float len = sqrt(ux*ux + uy*uy);
    return FMIN(-r, FMAX(a, b)) + len;
}
#endif

/////////////////////
// Signed Distance Functions
/////////////////////
static float sdBox(__thread__ glm::vec3 const& p, glm::vec3 b)
{
	using namespace glm;

	const glm::vec3 d(glm::abs(p) - b);
	const float a(FMIN(FMAX(d.x, FMAX(d.y, d.z)), 0.0f));
	const glm::vec3 d3(
		FMAX(d.x, 0.0f),
		FMAX(d.y, 0.0f),
		FMAX(d.z, 0.0f)
	);
	return FLENGTH(d3) + a;
}

static float sdSphere(__thread__ glm::vec3 const& p, float r)
{
	using namespace glm;

	const float l(FLENGTH(p));
	const float sd(l - r);
	return sd;
}

static float sdSquishengon(__thread__ glm::vec3 const& pt, __thread__ glm::vec3 const& dim, float thickness, float cornerRadiusX, float sides)
{
    using namespace glm;
    // brute force :-)
    const int n = int(mix(3.0f, 10.0f, sides));
    float a = 6.283185307179586f/n;
    float sdm = 1e36f;
    glm::vec3 q = glm::vec3(pt.y, pt.x, pt.z);
    float s = 1.0f;
    float md = min(min(dim.x, dim.y), dim.z);
    float cr = md * cornerRadiusX;
    //float th = 1e36f;
    //float phi = 1.0f/(2.0f*tan(a*0.75f));
    
    for (int i=0; i<n; i++)
    {
        float th0 = (float(i))*a;
        float th1 = (fmod(float(i + 1), float(n))) * a;
        glm::vec2 p0 = glm::vec2( cos(th0)*(dim.y-cr), sin(th0)*(dim.x-cr) );
        glm::vec2 p1 = glm::vec2( cos(th1)*(dim.y-cr), sin(th1)*(dim.x-cr) );
        glm::vec2 e = p0 - p1;
        glm::vec2 w = glm::vec2(q.x - p1.x, q.y - p1.y);
        glm::vec2 b = w - e * clamp(dot(w, e) / dot(e, e), 0.0f, 1.0f);
        sdm = FMIN(sdm, dot(b, b));
        //th = FMIN(th, length(e)*phi);
        float c = float(q.y>=p0.y) + float(q.y<p1.y)*2.0f + float(e.x*w.y>e.y*w.x)*4.0f;
        s *= (c==0.0f || c==7.0f) ? -1.0f : 1.0f;
    }
    
    sdm = s*FSQRT(sdm);
    
    //th *= thickness;
    //sdm = abs(sdm + th) - th;
    
    glm::vec2 w = glm::vec2( sdm, abs(q.z) - (dim.z-cr) );
    glm::vec2 dd = glm::vec2( FMAX(w.x, 0.0f), FMAX(w.y, 0.0f) );
    sdm = FMIN(FMAX(w.x, w.y), 0.0f) + FLENGTH(dd);
    sdm = sdm - cr;
    return sdm;
}

static inline glm::vec2 bezier(const __thread__ glm::vec3& p0, const __thread__ glm::vec3& p1, const __thread__ glm::vec3& p2)
{
    using namespace glm;
    const glm::vec3 b = (-2.0f)*p1 + p2;
    const glm::vec3 c = p1 * 2.0f;
    const glm::vec3 d = - p0;

    const float kk = 1.0f/(dot(b, b));
    const float kx = kk * dot(p1, b);
    const float ky = kk * (2.0f*dot(p1, p1)+dot(d, b)) * 0.33333333333f;
    const float kz = kk * dot(d, p1);

    glm::vec2 res;

    float p = ky - kx*kx;
    float p3 = p*p*p;
    float q = kx*(2.0f * kx * kx - 3.0f * ky) + kz;
    float h = q * q + 4.0f * p3;

    if (h >= 0.0f)
    {
        h = FSQRT(h);
        const glm::vec2 x = (glm::vec2(h, -h) - q) / 2.0f;
        const glm::vec2 uv = sign(x)*glm::pow(abs(x), glm::vec2(0.33333333333f));
        float t = uv.x + uv.y - kx;
        t = clamp( t, 0.0f, 1.0f );

        // 1 root
        glm::vec3 qos = d + (c + b*t)*t;
        res = glm::vec2(FLENGTH(qos), t);
    }
    else
    {
        const float z = FSQRT(-p);
        const float v = acos(q * 1.0f/(p*z*2.0f)) * 0.33333333333f;
        const float m = cos(v);
        const float n = sin(v)*1.732050808f;
        glm::vec3 t = glm::vec3(m + m, -n - m, n - m) * z - kx;
        t = clamp( t, 0.0f, 1.0f );

        // 3 roots
        glm::vec3 qos = d + (c + b*t.x)*t.x;
        float dis = dot(qos, qos);
        
        res = glm::vec2(dis, t.x);

        qos = d + (c + b * t.y)*t.y;
        dis = dot(qos, qos);
        if (dis < res.x) res = glm::vec2(dis, t.y);

        qos = d + (c + b*t.z)*t.z;
        dis = dot(qos, qos);
        if (dis<res.x) res = glm::vec2(dis, t.z);

        res.x = FSQRT(res.x);
    }
    
    return res;
}

static float sdBezier(__thread__ glm::vec3 const& pt, __thread__ glm::vec3 const& dim, float __thread__ const& thickness, float __thread__ const& cornerRadiusX, float __thread__ const& cornerRadiusY)
{
    using namespace glm;
    
    float minxyz = FMIN(FMIN(dim.x, dim.y) * 0.5f, dim.z);
    float fatL = FMIN(minxyz, mix(minxyz*2.0f, 0.0f, cornerRadiusY))*thickness;
    float fatR = FMIN(minxyz, mix(minxyz*2.0f, 0.0f, 1.0f-cornerRadiusY))*thickness;
    float fatC = mix(fatL, fatR, 0.5f);
    glm::vec3 p0 = glm::vec3(
        -(dim.x - fatL),
        -(dim.y - fatL),
        0.0f
    );
    
    glm::vec3 p2 = glm::vec3(
        (dim.x - fatR),
        -(dim.y - fatR),
        0.0f
    ) - p0;

    glm::vec3 p1 = glm::vec3(
        mix(-dim.x + fatC, dim.x - fatC , cornerRadiusX),
        FMAX(p0.y, dim.y*3.0f-FMAX(fatL, fatR)*4.0f),
        0.0f
    ) - p0;
    
    glm::vec3 p = pt;
    float elongationLimit = FMAX(0.0f, dim.z - max(fatL, fatR));
    p.z -= clamp(p.z, -elongationLimit, elongationLimit);

    glm::vec2 distParam = bezier(
        p-p0,
        p1,
        p2
    );

    return distParam.x - mix(fatL, fatR, distParam.y);
}

static inline glm::vec2 solve_quadratic(const __thread__ float& a, const __thread__ float& b, const __thread__ float& c)
{
    using namespace glm;
    // ax^2 + bx + c = 0, a non-zero
    const float q = b*b - 4.0f*a*c;
    if (q < 0.0f)
    {
        return glm::vec2(1e20f);
    }
    
    const float aDoubleDiv = 1.0f/(2.0f * a);
    float r0 = -b * aDoubleDiv;
    float r1 = sqrt(q) * aDoubleDiv;
    return glm::vec2(r0 - r1, r0 + r1);
}

// Media Molecule's max norm sphere: https://www.shadertoy.com/view/Mt2XWG
inline float sdMaxNormSphere(glm::vec3 __thread__ const& p, glm::vec3 __thread__ const& r)
{
    using namespace glm;
    
    // move ellipse so that target point is at origin, centre in positive space
    // f(v) = (v.x - c.x)^2/r.x^2 + (v.y - c.y)^2/r.y^2 + (v.z - c.z)^2/r.z^2 - 1
    const glm::vec3 c = abs(p);
    const glm::vec3 c2 = c*c;
    const glm::vec3 r2 = r*r;
    float d = 1e20f;

    const glm::vec3 r2Div = 1.0f/r2;

    // gather terms of quadratic
    const glm::vec3 qa = r2Div;
    const glm::vec3 qb = -2.0f*c * r2Div;
    const glm::vec3 qc = c2 * r2Div;
    const float qcs = (qc.x + qc.y + qc.z) - 1.0f;

    // check corners:
    // solve f(v)=0 for v.x=v.y=v.z=t
    const glm::vec2 t0 = abs(solve_quadratic((qa.x+qa.y+qa.z), (qb.x+qb.y+qb.z), qcs));
    d = FMIN(FMIN(d, t0.x), t0.y);

    // interior of convex shape always hits corners first, so early out
    if (qcs <= 0.0f)
    {
        return -d;
    }

    // check edges:
    // df/dx=0 => v.x=c.x, solve f(v)=0 for v.x=c.x, v.y=v.z=t
    // then do the same for y and z cases
    glm::vec2 t = abs(solve_quadratic(qa.y + qa.z, qb.y + qb.z, qc.y + qc.z - 1.0f));
    d = FMIN(d, FMAX(FMIN(t.x, t.y), c.x));

    t = abs(solve_quadratic(qa.x + qa.z, qb.x + qb.z, qc.x + qc.z - 1.0f));
    d = FMIN(d, FMAX(FMIN(t.x, t.y), c.y));

    t = abs(solve_quadratic(qa.x + qa.y, qb.x + qb.y, qc.x + qc.y - 1.0f));
    d = FMIN(d, FMAX(FMIN(t.x, t.y), c.z));

    // check faces:
    // df/dx=df/dy=0 => v.xy=c.xy, so f(v)=0 => |v.z - c.z|=r.z
    d = FMIN(d, FMAX(FMAX(c.x, c.y), abs(c.z - r.z)));
    d = FMIN(d, FMAX(FMAX(c.x, abs(c.y - r.y)), c.z));
    d = FMIN(d, FMAX(FMAX(abs(c.x - r.x), c.y), c.z));

    // done
    return d;
}

static float sdEllipsoid(__thread__ glm::vec3 const& p, __thread__ glm::vec3 const& dim, float __thread__ const& thickness, float __thread__ const& cornerRadiusY, float __thread__ const& slice)
{
    using namespace glm;
    //const float round = cornerRadiusY - thickness;
            
    float dPerfectSphere = length(p) - dim.x;// + round;
    
    // mediamolecule's max norm sphere, needed for squashed spheroids
    const float dMaxNormSphere = sdMaxNormSphere(p, dim);// + round;
    
    const float maxScale = max(dim.x, max(dim.y, dim.z));
    const float minScale = min(dim.x, min(dim.y, dim.z));
    const float squash = smoothstep(0.9f, 1.0f, minScale/maxScale);
    const float d = mix(dMaxNormSphere, dPerfectSphere, squash);
    return d;
    
    // slicing
    //const float sliceY = (-dim.y * ((slice * 2.0f) - 1.0f));
    //const float slicedY = -p.y - sliceY;
    //return length(max(vec2(d, slicedY), 0.0f)) + min(max(d, slicedY), 0.0f) - round;
}


// s: width, height, depth, thickness
// r: xy corner radius, z corner radius
static float sdSuperprim(glm::vec3 __thread__ const& p, glm::vec3 __thread__ const& dim, float __thread__ const& thickness, float __thread__ const& cornerRadiusX, float __thread__ const& cornerRadiusY)
{
	using namespace glm;
	const float mdxy(FMAX(dim.x, dim.y));
	const float rx(cornerRadiusX*mdxy);
	const float ry(cornerRadiusY*dim.z);

	const glm::vec3 d(glm::abs(p.x) - dim.x,
                      glm::abs(p.y) - dim.y,
                      glm::abs(p.z) - dim.z);
	const glm::vec2 dr(
		FMAX(d.x + rx, 0.0f),
		FMAX(d.y + rx, 0.0f)
	);
	const float sw(thickness*mdxy);
	const float q(glm::abs(FLENGTH(dr) + FMIN(-rx, FMAX(d.x, d.y)) + sw) - sw);
	const glm::vec2 qq(FMAX(q + ry, 0.0f), FMAX(d.z + ry, 0.0f));
	return FLENGTH(qq) + FMIN(-ry, FMAX(q, d.z));
}

static float sdUberprim(glm::vec3 __thread__ const& p, glm::vec3 __thread__ const& dim, float __thread__ const& conef, float __thread__ const& cornerRadiusX, float __thread__ const& cornerRadiusY)
{
	using namespace glm;
    float mdxy = FMAX(dim.x, dim.y);
    float rx = cornerRadiusX * mdxy;
    float ry = cornerRadiusY * dim.z;

    const glm::vec3 dr (
        dim.x - rx,
        dim.y - rx, 
        dim.z - ry);
    
    const float cone = -(1.0f - conef)*FMIN(dim.x, dim.y);
    glm::vec2 ba = glm::vec2(cone, -2.0*dr.z);
    float sz2 = ba.y;
    ba /= dot(ba, ba);

    const glm::vec3 d(
        abs(p.x) - dr.x,
        abs(p.y) - dr.y,
        abs(p.z) - dr.z
    );
    const float q = FLENGTH(max(glm::vec2(d), 0.0f)) + FMIN(0.0f, FMAX(d.x, d.y)) - rx + ry;

    glm::vec2 pa = glm::vec2(q, p.z - dr.z);
    glm::vec2 diag = pa - glm::vec2(cone, sz2) * clamp(dot(pa, ba), 0.0f, 1.0f);
    glm::vec2 h0 = glm::vec2(FMAX(q - cone, 0.0f), p.z + dr.z);
    glm::vec2 h1 = glm::vec2(FMAX(q, 0.0f), p.z - dr.z);

    return FSQRT(FMIN(dot(diag, diag), FMIN(dot(h0, h0), dot(h1, h1)))) * sign(FMAX(dot(pa, glm::vec2(-ba.y, ba.x)), d.z)) - ry;
}
