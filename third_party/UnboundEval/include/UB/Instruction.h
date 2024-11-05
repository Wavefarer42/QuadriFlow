
#pragma once

#include "Shape.h"
#include "UB.h"
#include <vector>

namespace UB
{
    static float makeRotationBounds(const glm::quat& q)
    {
        using namespace glm;
        const vec3 qxyz(q.x, q.y, q.z);

        const vec3 w(q.w * qxyz);
        const vec3 x(q.x * qxyz);
        const vec3 y(q.y * qxyz);
        const float zz(q.z*q.z);

        const float a(fabs(0.5f - (y.y + zz)) + fabs(x.y - w.z) + fabs(x.z + w.y));
        const float b(fabs(x.y + w.z) + fabs(0.5f - (x.x + zz)) + fabs(y.z - w.x));
        const float c(fabs(x.z - w.y) + fabs(y.z + w.x) + fabs(0.5f - (x.x + y.y)));
        return 2.0f * max(a, max(b, c));
    }

    enum class InstructionType
    {
        NOOP = 0,
        EVAL = 1,
        PUSH = 2,
        POP  = 3
    };

    struct InstructionCommon
    {
        InstructionType type;
        int editIndex;
    };
    
    struct InstructionEval
    {
        InstructionType type;
        int editIndex;
        Math::Transform transform;
        Shape shape;
        EditOp op;
        
        glm::u16vec4 pack_orientation() const
        {
            glm::quat oq(glm::normalize(glm::inverse(transform.orientation)));
            return glm::u16vec4(
                glm::packHalf1x16(oq.x),
                glm::packHalf1x16(oq.y),
                glm::packHalf1x16(oq.z),
                glm::packHalf1x16(oq.w)
            );
        }

        glm::u16vec4 pack_center_blend() const
        {
            return glm::u16vec4(
                glm::packHalf1x16(transform.position.x),
                glm::packHalf1x16(transform.position.y),
                glm::packHalf1x16(transform.position.z),
                glm::packHalf1x16(op.blend)
            );
        }

        glm::u16vec4 pack_dim_bounds() const
        {
            return glm::u16vec4(
                glm::packHalf1x16(shape.dim.x),
                glm::packHalf1x16(shape.dim.y),
                glm::packHalf1x16(shape.dim.z),
                glm::packHalf1x16(makeRotationBounds(glm::normalize(transform.orientation)))
            );
        }

        glm::u8vec4 pack_cornerXY_thickness_op() const
        {
            return glm::u8vec4(
                glm::packUnorm1x8(shape.cornerRadius.x),
                glm::packUnorm1x8(shape.cornerRadius.y),
                glm::packUnorm1x8(shape.thickness),
                uint8_t((uint8_t(shape.type) << 4) | op.type)
            );
        }

        glm::u8vec4 pack_material() const
        {
            return glm::u8vec4(
                glm::clamp(shape.color.x, 0.0f, 1.0f)*255.0f,
                glm::clamp(shape.color.y, 0.0f, 1.0f)*255.0f,
                glm::clamp(shape.color.z, 0.0f, 1.0f)*255.0f,
                0
            );
        }

    };
    
    struct InstructionPush
    {
        InstructionType type;
        int editIndex;
    };
    
    struct InstructionPop
    {
        InstructionType type;
        int editIndex;
    };
    
    union Instruction
    {
        InstructionType type;
        InstructionCommon common;
        InstructionEval eval;
        InstructionPush push;
        InstructionPop pop;
        
        Instruction() : type(InstructionType::NOOP) {} // defaults to EditShape

        ~Instruction()
        {
        }
        
        void operator=(const Instruction& other)
        {
            switch (other.type)
            {
                case InstructionType::EVAL: eval=other.eval; break;
                case InstructionType::PUSH: push=other.push; break;
                case InstructionType::POP: pop=other.pop; break;
                case InstructionType::NOOP: break;
                default: break;
            }
        }

        Instruction(const Instruction& other)
        {
            *this = other;
        }
    };
    
    using InstructionList = std::vector<Instruction>;
    using InstructionListRef = std::shared_ptr<InstructionList>;
    
    bool compileEditList(const UB::EditList& source, InstructionList& destination);
    uint64_t hashEditList(const UB::EditList& source);
    
    glm::range evalRange(const glm::vec3& p, const float nodeHalfWidth, const InstructionEval &e);
    glm::range evalAppendRange(const glm::range& prevSd, const glm::range& sd, const int& op, const float& blend);
    
    float evalDistance(const glm::vec3& p, const InstructionEval& e);
    float evalDistanceLocal(const glm::vec3& p, const InstructionEval& e);
    float evalDistanceAndGradient(const glm::vec3& p, glm::vec3 &grad, const InstructionList& evals, float epsi = 0.001f);
    glm::range appendRange(const glm::range& prevSd, const glm::range& sd, const int& op, const float& blend);
    
    float appendDistance(const float prevSd, const float sd, const int op, const float blend);
    float evalDistance(const glm::vec3& p, const InstructionList& evals);
    glm::vec3 evalGradient(const glm::vec3& p, const InstructionList& evals);
    glm::vec3 evalAlbedo(const glm::vec3& p, const InstructionList& evals);
    
    int getClosestEditIndex(const glm::vec3& p, const InstructionList& evals);
    bool findClosestIntersection(const glm::vec3& pos, const glm::vec3 &dir, glm::vec3& pp, const InstructionList& evals);
    bool findClosestIntersection(const glm::vec3& pos, const glm::vec3 &dir, glm::vec3& pp, float tmin, float tmax, const InstructionList& evals);
    
    glm::range3 getAABB(const InstructionEval& inst);
    glm::range3 getAABB(const InstructionList& evals);
    glm::range3 getTightAABB(const InstructionList& evals);
    uint64_t getUUID(const InstructionList& evals);

     
}
