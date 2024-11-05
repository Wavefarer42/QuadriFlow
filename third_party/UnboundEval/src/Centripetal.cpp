//
// Unbound - Centripetal.cpp
// Copyright (c) 2017 Unbound Technologies, Inc. All rights reserved.
// 

#include "Centripetal.h"


Math::CubicPoly Math::makeNonUniformCatmullRom(float x0, float x1, float x2, float x3, float dt0, float dt1, float dt2)
{
    // compute tangents when parameterized in [t1,t2]
    float t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
    float t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;

    // rescale tangents for parametrization in [0,1]
    t1 = t1 * dt1;
    t2 = t2 * dt1;

    return Math::CubicPoly(x1, x2, t1, t2);
}

Math::Centripedal3::Centripedal3(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
{
    float dt0 = glm::pow(glm::distance2(p0, p1), 0.25f);
    float dt1 = glm::pow(glm::distance2(p1, p2), 0.25f);
    float dt2 = glm::pow(glm::distance2(p2, p3), 0.25f);

    const float eps = 0.1f;
    // safety check for repeated points
    if (dt1 < eps) { dt1 = 1.0f; }
    if (dt0 < eps) { dt0 = dt1; }
    if (dt2 < eps) { dt2 = dt1; }

    x = makeNonUniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2);
    y = makeNonUniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2);
    z = makeNonUniformCatmullRom(p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2);
}

glm::vec3 Math::Centripedal3::eval(float t) const
{
    return glm::vec3(x.eval(t), y.eval(t), z.eval(t));
}

float Math::Centripedal3::getLength() const
{
    int num_pieces = 50;
    glm::vec3 last_pos = eval(0.0f);
    float d = 0.0;
    for (int i = 0; i < num_pieces; i++)
    {
        glm::vec3 pos = eval(float(i) / float(num_pieces));
        d = d + distance(last_pos, pos);
        last_pos = pos;
    }
    return d;
}


