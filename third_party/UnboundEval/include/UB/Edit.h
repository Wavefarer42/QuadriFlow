
#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <memory>
#include <variant>

#include "Shape.h"
#include "Math/Centripetal.h"
#include "Math/Transform.h"

namespace UB
{
    struct EditList;
    using EditListRef = std::shared_ptr<EditList>;
        
    struct EditOp
    {
        enum Type
        {
            Sub   = 0,
            Add   = 1,
            Color = 2,
            Crop  = 3
        } type = Add;
        
        float blend = 0.0f;
        bool mirrorYZ = false;
        
        EditOp() = default;
        EditOp(const EditOp&) = default;
        EditOp(Type type) :type(type) {}
        EditOp(Type type, float blend) :type(type), blend(blend) {}
    };


    struct EditShape
    {
        EditOp op;
        Math::Transform transform;
        Shape shape;
    };

    struct EditGroup
    {
        EditOp op;
        Math::Transform transform;
        EditListRef edits;
    };

    struct EditRepeat
    {
        EditOp op;
        Math::Transform transform;
        Shape shape;
        int reps = 1;
        Math::Transform offsetTransform;
    };
    
    struct EditSpline
    {
        EditOp op;
        Math::Transform transform;
        bool loop = false;
        int reps = 16;
        struct ControlPoint
        {
            Math::Transform transform;
            Shape shape;
            ControlPoint() {};
            ControlPoint(const ControlPoint& other) : shape(other.shape), transform(other.transform) {}
            ControlPoint(const Math::Transform& transform) : transform(transform) {}
            ControlPoint(const glm::vec3& position) : transform(position) {}
            ControlPoint(const Math::Transform& transform, const UB::Shape& shape) : transform(transform), shape(shape) {}
            void operator=(const ControlPoint& other)
            {
                shape = other.shape;
                transform = other.transform;
            }
            ControlPoint& operator=(ControlPoint&& other)
            {
                shape = std::move(other.shape);
                transform = std::move(other.transform);
                return *this;
            }

        };
        std::vector<ControlPoint> controlpoints;
    };

    struct EditFilter
    {
        EditOp op;
    };

    struct EditGenerator
    {
        EditOp op;
    };

    using Edit = std::variant<std::monostate, EditShape, EditGroup, EditRepeat, EditSpline, EditFilter, EditGenerator>;
    
    struct EditList
    {
        static const uint64_t HASHSEED = 0x1374200094857563llu;
        std::vector<Edit> edits;
    
        EditList() = default;
        EditList(const EditList&) = default;

        inline void append(const Edit& e)           { edits.push_back(e); }
        inline void append(const EditShape& e)      { edits.push_back(e); }
        inline void append(const EditGroup& e)      { edits.push_back(e); }
        inline void append(const EditRepeat& e)     { edits.push_back(e); }
        inline void append(const EditSpline& e)     { edits.push_back(e); }
        inline void append(const EditFilter& e)     { edits.push_back(e); }
        inline void append(const EditGenerator& e)  { edits.push_back(e); }

        inline void insert(const Edit& e, int pos)          { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditShape& e, int pos)     { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditGroup& e, int pos)     { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditRepeat& e, int pos)    { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditSpline& e, int pos)    { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditFilter& e, int pos)    { edits.insert(edits.begin()+pos, e); }
        inline void insert(const EditGenerator& e, int pos) { edits.insert(edits.begin()+pos, e); }
        
        inline void replace(const Edit& e, int pos)             { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditShape& e, int pos)        { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditGroup& e, int pos)        { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditRepeat& e, int pos)       { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditSpline& e, int pos)       { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditFilter& e, int pos)       { if (pos < edits.size()) edits[pos] = e; }
        inline void replace(const EditGenerator& e, int pos)    { if (pos < edits.size()) edits[pos] = e; }

        void remove(int index);
        void remove(int indexStart, int indexEnd);
        void clear();
        
        inline Edit& get(int index) { return edits[index]; }
        inline const Edit& get(int index) const { return edits[index]; }
        inline const Edit& getConst(int index) const { return edits[index]; }
        inline int size() const { return (int)edits.size(); }
        void operator=(const EditList& other) { edits = other.edits; }
        uint64_t computeHash() const;

        friend std::ostream& operator<<(std::ostream& stream, const EditList& list);
    };
}
