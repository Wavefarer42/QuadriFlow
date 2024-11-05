#define XXH_STATIC_LINKING_ONLY
#include <xxhash.h>

#include "Edit.h"
#include "Util/FloatControl.h"

namespace UB
{

void EditList::remove(int index)
{
    if (index < edits.size())
    {
        edits.erase(edits.begin()+index);
    }
}

void EditList::remove(int indexStart, int indexEnd)
{
    for (int i=indexStart; i<indexEnd; i++)
    {
        remove(i);
    }
}

void EditList::clear()
{
    edits.clear();
}

uint64_t EditList::computeHash() const
{
    Util::SetAndRestoreFloatControlDownward fpCtrl;
    XXH64_state_t *const state = XXH64_createState();
    XXH64_hash_t const seed = HASHSEED;
    
    for (int i=0; i<size(); i++)
    {
        const UB::Edit& ee = getConst(i);
        
        if (std::holds_alternative<EditShape>(ee))
        {
            const EditShape& e = std::get<EditShape>(ee);
            XXH64_update(state, &e, sizeof(EditShape));
        }
        else if (std::holds_alternative<EditGroup>(ee))
        {
            const EditGroup& e = std::get<EditGroup>(ee);
            XXH64_update(state, &e, sizeof(EditGroup));
        }
        else if (std::holds_alternative<EditRepeat>(ee))
        {
            const EditRepeat& e = std::get<EditRepeat>(ee);
            XXH64_update(state, &e, sizeof(EditRepeat));
        }
        else if (std::holds_alternative<EditSpline>(ee))
        {
            const EditSpline& e = std::get<EditSpline>(ee);
            XXH64_update(state, &e.transform, sizeof(e.transform));
            XXH64_update(state, &e.op, sizeof(e.op));
            XXH64_update(state, &e.loop, sizeof(e.loop));
            XXH64_update(state, &e.reps, sizeof(e.reps));
            for (auto& cp : e.controlpoints)
            {
                XXH64_update(state, &cp, sizeof(EditSpline::ControlPoint));
            }
            
        }
        else if (std::holds_alternative<EditFilter>(ee))
        {
            const EditFilter& e = std::get<EditFilter>(ee);
            XXH64_update(state, &e, sizeof(EditFilter));
        }
        else if (std::holds_alternative<EditGenerator>(ee))
        {
            const EditGenerator& e = std::get<EditGenerator>(ee);
            XXH64_update(state, &e, sizeof(EditGenerator));
        }
    }
    
    uint64_t listHash = XXH64_digest(state);
    XXH64_freeState(state);
    return listHash;
}

std::ostream& operator<<(std::ostream& stream, const std::monostate& _)
{
    assert(false && "std::monostate should never actually occur!");
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditOp& edit)
{
    switch (edit.type) 
    {
        case EditOp::Type::Sub: stream << "sub"; break;
        case EditOp::Type::Add: stream << "add"; break;
        case EditOp::Type::Color: stream << "color"; break;
        case EditOp::Type::Crop: stream << "crop"; break;
        default: stream << "<Type>"; break;
    }

    stream << ", " << edit.blend << ", mirror: " << (edit.mirrorYZ ? "true" : "false");

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditShape& edit)
{
    stream << edit.op;
    stream << edit.transform;
    stream << edit.shape;

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditGroup& _)
{
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditRepeat& _)
{
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditSpline& _)
{
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditFilter& _)
{
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditGenerator& _)
{
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const EditList& list)
{
    for (int i=0; i<list.size(); i++)
    {
        const Edit& edit = list.getConst(i);
        std::visit([&stream](auto&& edit) { stream << edit; }, edit);
    }

    return stream;
}

} // namespace UB
