
//
// Unbound Framework
// Copyright (c) 2014-current Unbound Technologies, Inc. All rights reserved.
//

#pragma once
#include <vector>
#include <memory>
#include <sstream>
#include <iomanip>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/packing.hpp>
#include <glm/gtx/component_wise.hpp>

#include "CSG/IA.h"
#include "CSG/SDF.h"

namespace UB
{
    struct Shape
    {
        enum class Type
        {
            SUPERPRIM   = 0,
            ELLIPSOID   = 1,
            BEZIER      = 2,
            SQUISHENGON = 3,
            UBERPRIM    = 4
        };
        
        Type type = Type::SUPERPRIM;
        glm::vec3 dim = glm::vec3(1);
        float thickness = 1.0f;
        glm::vec2 cornerRadius = glm::vec2(0.0f);
        glm::vec3 color = glm::vec3(1.0f);
        
        Shape() = default;
        Shape(const Shape& other) = default;

        bool operator==(const Shape & other) const
        {
            return std::memcmp(this, &other, sizeof(Shape)) == 0;
        }

        bool operator!=(const Shape & other) const
        {
            return std::memcmp(this, &other, sizeof(Shape)) != 0;
        }

        void switchSphere();
        void switchBox();
        void switchCapsule();
        void switchPellet();
        void switchTorus();
        void switchRoundedBox();
        void switchCylinder();

        static Shape Sphere(float radius);
        static Shape Box(float dx, float dy, float dz);
        static Shape Box(const glm::vec3& dim);
        static Shape Cylinder(float radius, float length);
    };

    std::ostream& operator<<(std::ostream& o, const UB::Shape& e);
    std::ostream& operator<<(std::ostream& o, const UB::Shape* e);
    
}
