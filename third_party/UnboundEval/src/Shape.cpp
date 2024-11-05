
//
// Unbound Framework
// Copyright (c) 2014-current Unbound Technologies, Inc. All rights reserved.
//


#include "Shape.h"
#include <iostream>
#define XXH_STATIC_LINKING_ONLY
#include <xxhash.h>

using namespace UB;


void Shape::switchSphere()
{
    dim = glm::vec3(glm::compMax(dim));
    cornerRadius = glm::vec2(1.0f);
    thickness = 1.0f;
}

void Shape::switchBox()
{
    cornerRadius = glm::vec2(0.0f, 0.0f);
    thickness = 1.0f;
}

void Shape::switchCapsule()
{
    float r = glm::compMax(dim);
    dim = glm::vec3(r*0.5f, r*0.5f, r);
    cornerRadius = glm::vec2(1.0f,0.5f);
    thickness = 1.0f;
}

void Shape::switchPellet()
{
    float r = glm::compMax(dim);
    dim = glm::vec3(r, r, r*0.25f);
    cornerRadius = glm::vec2(1.0f, 0.25f);
    thickness = 1.0f;
}

void Shape::switchTorus()
{
    float r = glm::compMax(dim);
    dim = glm::vec3(r, r, r*0.25f);
    cornerRadius = glm::vec2(1.0f, 1.0f);
    thickness = 0.25f;
}

void Shape::switchRoundedBox()
{
    cornerRadius = glm::vec2(0.25f, 0.25f);
    thickness = 1.0f;
}

void Shape::switchCylinder()
{
    cornerRadius = glm::vec2(1.0, 0.0f);
    thickness = 1.0f;
}

Shape UB::Shape::Sphere(float radius)
{
    Shape e;
    e.dim = glm::vec3(radius);
    e.cornerRadius = glm::vec2(1.0f);
    e.thickness = 1.0f;
    return e;
}

Shape UB::Shape::Box(float dx, float dy, float dz)
{
    Shape e;
    e.dim = glm::max(glm::vec3(dx, dy, dz), 0.0f);
    e.cornerRadius = glm::vec2(0.0f);
    e.thickness = 1.0f;
    return e;
}

Shape UB::Shape::Box(const glm::vec3 & dim)
{
    Shape e;
    e.dim = glm::max(dim, 0.0f);
    e.cornerRadius = glm::vec2(0.0f);
    e.thickness = 1.0f;
    return e;
}

Shape UB::Shape::Cylinder(float radius, float length)
{
    Shape e;
    e.dim = glm::max(glm::vec3(radius, radius, length), 0.0f);
    e.cornerRadius = glm::vec2(1.0f, 0.0f);
    e.thickness = 1.0f;
    return e;
}

std::ostream& UB::operator<<(std::ostream& o, const UB::Shape& e)
{
    o << "\tdim: (" << e.dim.x << ", " << e.dim.y << ", " << e.dim.z << ")" << std::endl;
    o << "\tcornerRadius: (" << e.cornerRadius.x << ", " << e.cornerRadius.y << ")" << std::endl;
    o << "\tthickness: " << e.thickness << std::endl;
    o << "\tcolor: " << (int)e.color.x << ", " << (int)e.color.y << ", " << (int)e.color.z << std::endl;
    return o;
}

std::ostream& operator<<(std::ostream& o, const UB::Shape* e)
{
    return UB::operator<<(o, *e);
}
