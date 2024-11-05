//
// Unbound - Transform.h
// Copyright (c) 2016 Unbound Technologies, Inc. All rights reserved.
// 

#pragma once

#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/hash.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/io.hpp>

#include "Util/CommonFunc.h"
#include "CSG/IA.h"

#include <functional>

namespace Math
{
    struct Transform
    {
        glm::quat	 orientation;
        glm::vec3	 position;
        glm::vec3	 scale;

        Transform(const Transform& other) :
            position(other.position),
            orientation(other.orientation),
            scale(other.scale)
        {}

        Transform() :
            position(0.0f, 0.0f, 0.0f),
            scale(1.0f),
            orientation()
        {}

        explicit Transform(float scale) :
            position(0.0f, 0.0f, 0.0f),
            scale(scale),
            orientation()
        {}

        explicit Transform(const glm::vec3& position) :
            position(position),
            scale(1.0f),
            orientation()
        {}

        explicit Transform(const glm::mat4& pose) :
            position(glm::vec3(glm::column(pose, 3))),
            scale(1.0f),
            orientation(glm::quat_cast(pose))
        {}

        Transform(const glm::vec3& position, const glm::quat& orientation) :
            position(position),
            scale(1.0f),
            orientation(orientation)
        {}

        Transform(const glm::vec3& position, const glm::quat& orientation, const float scale) :
            position(position),
            scale(scale),
            orientation(orientation)
        {}

        Transform(const glm::vec3& position, const glm::quat& orientation, const glm::vec3& scale) :
            position(position),
            scale(scale),
            orientation(orientation)
        {}

        Transform(const glm::vec3& position, const glm::vec3& direction, const glm::vec3& up)
        {
            this->orientation = Util::lookAtOrientation(glm::normalize(direction), glm::normalize(up));
            this->position = position;
            this->scale = glm::vec3(1.0f);
        }

        glm::mat4 getMatrix()
        {
            return glm::translate(position) * glm::mat4_cast(orientation) * glm::scale(glm::vec3(scale));
        }

        void applyParent(const Transform& parent)
        {
            orientation = parent.orientation * orientation;
            position = parent.orientation * (position * parent.scale) + parent.position;
            scale = scale * parent.scale;
        }

        void setRelativeTo(const Transform& parent)
        {
            const glm::vec3 scaleInv = 1.0f / parent.scale;
            const glm::quat invOrientation(glm::inverse(parent.orientation));

            scale = scale * scaleInv;
            position = (invOrientation * (position - parent.position)) * scaleInv;
            orientation = invOrientation * orientation;
        }

        Transform getRelativeTo(const Transform& other)
        {
            const glm::vec3 scaleInv = 1.0f / other.scale;
            const glm::quat invOrientation(glm::inverse(other.orientation));
            Transform ret;
            ret.scale = scale * scaleInv;
            ret.position = (invOrientation * (position - other.position)) * scaleInv;
            ret.orientation = invOrientation * orientation;
            return ret;
        }

        void mirrorXY()
        {
            orientation.x = -orientation.y;
            orientation.y = -orientation.y;
            position.z = -position.z;
        }

        void mirrorYZ()
        {
            orientation.y = -orientation.y;
            orientation.z = -orientation.z;
            position.x = -position.x;
        }

        void mirrorXZ()
        {
            orientation.x = -orientation.x;
            orientation.z = -orientation.z;
            position.y = -position.y;
        }

        glm::vec3 positionToLocal(const glm::vec3& pos) const
        {
            const glm::vec3 scaleInv = 1.0f / scale;
            const glm::quat invOrientation(glm::inverse(orientation));
            return (invOrientation * (pos - position)) * scaleInv;
        }

        glm::quat orientationToLocal(const glm::quat& rot) const
        {
            const glm::quat invOrientation(glm::inverse(orientation));
            return invOrientation * rot;
        }

        glm::vec3 directionToLocal(const glm::vec3& dir) const
        {
            const glm::quat invOrientation(glm::inverse(orientation));
            const glm::vec3 scaleInv = 1.0f / scale;
            return glm::normalize((invOrientation * dir) * scaleInv);
        }

        void operator=(const glm::mat4& m)
        {
            position = glm::vec3(glm::column(m, 3));
            orientation = glm::quat_cast(m);
        }

        void yawBy(float a)
        {
            glm::vec3 axis = orientation * glm::vec3(0.0f, 1.0f, 0.0f);
            orientation = glm::angleAxis(a, axis) * orientation;
        }

        void pitchBy(float a)
        {
            glm::vec3 axis = orientation * glm::vec3(1.0f, 0.0f, 0.0f);
            orientation = glm::angleAxis(a, axis) * orientation;
        }

        void rollBy(float a)
        {
            glm::vec3 axis = orientation * glm::vec3(0.0f, 0.0f, -1.0f);
            orientation = glm::angleAxis(a, axis) * orientation;
        }

        glm::vec3 forward() const
        {
            return orientation * glm::vec3(0.0f, 0.0f, -1.0f);
        }

        glm::vec3 backward() const
        {
            return orientation * glm::vec3(0.0f, 0.0f, 1.0f);
        }

        glm::vec3 up() const
        {
            return orientation * glm::vec3(0.0f, 1.0f, 0.0f);
        }

        glm::vec3 down() const
        {
            return orientation * glm::vec3(0.0f, -1.0f, 0.0f);
        }

        glm::vec3 right() const
        {
            return orientation * glm::vec3(1.0f, 0.0f, 0.0f);
        }

        glm::vec3 left() const
        {
            return orientation * glm::vec3(-1.0f, 0.0f, 0.0f);
        }

        size_t hash() const noexcept
        {
            return Util::hash(position, scale, orientation);
        }

        friend std::ostream& operator<<(std::ostream& stream, const Math::Transform& xform)
        {
            const auto angles = glm::eulerAngles(xform.orientation);

            stream << xform.position.x << ", " << xform.position.y << ", " << xform.position.z << " | ";
            stream << xform.scale << " | ";
            stream << angles.x << ", " << angles.y << ", " << angles.z;

            return stream;
        }
    };

    static Transform inverse(const Transform& t)
    {
        return Transform(-t.position, glm::inverse(t.orientation), 1.0f / t.scale);
    }

    static Transform operator*(const Transform& a, const Transform& b)
    {
        return Transform(
            a.orientation * (b.position * a.scale) + a.position,
            a.orientation * b.orientation,
            a.scale * b.scale
        );
    }

    static glm::vec3 operator*(const Transform& t, const glm::vec3& p)
    {
        return t.orientation * (p * t.scale) + t.position;
    }

    static bool operator==(const Transform& t0, const Transform& t1)
    {
        return glm::all(glm::equal(t0.position, t1.position)) and 
               glm::all(glm::equal(t0.orientation, t1.orientation)) and 
               glm::all(glm::equal(t0.scale, t1.scale));
    }

    static bool transformRoughlyEqual(const Transform& t0, const Transform& t1)
    {
        return glm::all(glm::lessThanEqual(glm::abs(t0.position - t1.position), glm::vec3(FLT_EPSILON))) and
               (glm::abs(t0.orientation.x - t1.orientation.x) < FLT_EPSILON) and
               (glm::abs(t0.orientation.y - t1.orientation.y) < FLT_EPSILON) and
               (glm::abs(t0.orientation.z - t1.orientation.z) < FLT_EPSILON) and
               (glm::abs(t0.orientation.w - t1.orientation.w) < FLT_EPSILON) and
               glm::all(glm::lessThanEqual(glm::abs(t0.scale - t1.scale), glm::vec3(FLT_EPSILON)));
    }

    static Transform mix(Transform const& x, Transform const& y, float a)
    {
        Transform ret;
        ret.position = glm::mix(x.position, y.position, a);
        ret.orientation = glm::slerp(x.orientation, y.orientation, a);
        ret.scale = glm::mix(x.scale, y.scale, a);
        return ret;
    }
    
    static glm::range3 operator*(const Transform& t, const glm::range3& aabb)
    {
        glm::vec3 lb = glm::vec3(aabb.x.lb, aabb.y.lb, aabb.z.lb) * t.scale;
        glm::vec3 ub = glm::vec3(aabb.x.ub, aabb.y.ub, aabb.z.ub) * t.scale;
        const glm::vec3 right = t.right();
        const glm::vec3 up = t.up();
        const glm::vec3 forward = -t.forward();
        const glm::vec3 xa = right * lb.x;
        const glm::vec3 xb = right * ub.x;
        const glm::vec3 ya = up * lb.y;
        const glm::vec3 yb = up * ub.y;
        const glm::vec3 za = forward * lb.z;
        const glm::vec3 zb = forward * ub.z;
        lb = glm::min(xa, xb) + glm::min(ya, yb) + glm::min(za, zb) + t.position;
        ub = glm::max(xa, xb) + glm::max(ya, yb) + glm::max(za, zb) + t.position;
        return glm::range3(lb, ub);
    }

}

namespace std
{

    template<>
    struct hash<Math::Transform>
    {
        size_t operator()(const Math::Transform& xform) const noexcept
        {
            return xform.hash();
        }
    };

} // namespace std
