#include "Instruction.h"
#define XXH_STATIC_LINKING_ONLY
#include "xxhash.h"

//
// compiles InstructionList from EditList
// return false if there is guaranteed no surface (only color ops or only negative ops)
//
bool UB::compileEditList(const UB::EditList& source, UB::InstructionList& destination)
{
    using namespace glm;
    using namespace Math;
    
    destination.clear();
    destination.reserve(source.size());
    bool hasAdd = false;
    
    for (int i=0; i<source.size(); i++)
    {
        // save check for our max edit list length
        if (i > 65535) break;

        const UB::Edit& ee = source.getConst(i);
        Instruction inst;
        
        if (std::holds_alternative<EditShape>(ee))
        {
            //
            // Adds a single shape
            //
            const EditShape& e = std::get<EditShape>(ee);
            inst.type = InstructionType::EVAL;
            inst.eval.transform = e.transform;
            inst.eval.op = e.op;
            inst.eval.shape = e.shape;
            inst.eval.editIndex = i;
            if (e.op.type == EditOp::Add) hasAdd = true;
            destination.push_back(inst);
            if (e.op.mirrorYZ == true)
            {
                if (inst.eval.shape.type == UB::Shape::Type::BEZIER)
                {
                    inst.eval.shape.cornerRadius = 1.0f - inst.eval.shape.cornerRadius;
                }
                inst.eval.transform.mirrorYZ();
                destination.push_back(inst);
            }

        }
        else if (std::holds_alternative<EditGroup>(ee))
        {
            // TODO
            const EditGroup& e = std::get<EditGroup>(ee);
            if (e.op.type == EditOp::Add) hasAdd = true;
        }
        else if (std::holds_alternative<EditRepeat>(ee))
        {
            //
            // Adds a repeated chain of shapes via transform multiplication
            //
            const EditRepeat& e = std::get<EditRepeat>(ee);
            if (e.op.type == EditOp::Add) hasAdd = true;
            inst.type = InstructionType::EVAL;
            inst.eval.op = e.op;
            inst.eval.editIndex = i;
            inst.eval.transform = e.transform;
            inst.eval.shape = e.shape;
            inst.eval.shape.dim *= e.transform.scale;
            destination.push_back(inst);
            if (e.op.mirrorYZ == true)
            {
                if (inst.eval.shape.type == UB::Shape::Type::BEZIER)
                {
                    inst.eval.shape.cornerRadius = 1.0f - inst.eval.shape.cornerRadius;
                }

                inst.eval.transform.mirrorYZ();
                destination.push_back(inst);
            }
            Transform ct = e.transform;
            
            for (int j=0; j<e.reps; j++)
            {
                ct.position += ct.orientation * e.offsetTransform.position * e.transform.scale;
                ct.orientation = ct.orientation * e.offsetTransform.orientation;
                
                inst.eval.transform.position = ct.position;
                inst.eval.transform.orientation = ct.orientation;
                inst.eval.shape.dim *= e.offsetTransform.scale;
                destination.push_back(inst);
                if (e.op.mirrorYZ == true)
                {
                    if (inst.eval.shape.type == UB::Shape::Type::BEZIER)
                    {
                        inst.eval.shape.cornerRadius = 1.0f - inst.eval.shape.cornerRadius;
                    }

                    inst.eval.transform.mirrorYZ();
                    destination.push_back(inst);
                }

            }
        }
        else if (std::holds_alternative<EditSpline>(ee))
        {
            //
            // Put edits along a spline
            //
            const EditSpline& e = std::get<EditSpline>(ee);
            
            if (e.controlpoints.size() < 2)
                continue;
            
            inst.type = InstructionType::EVAL;
            inst.eval.op = e.op;
            if (e.op.type == EditOp::Add) hasAdd = true;
            inst.eval.editIndex = i;
            inst.eval.transform = e.transform;
            //inst.eval.shape = e.spline.shape;
            int numSegments = e.controlpoints.size() - 1 + int(e.loop);
            int lastIndex = e.controlpoints.size()-1;
            for (int j=0; j<numSegments; j++)
            {
                int idx0 = e.loop ? ((j-1 < 0) ? lastIndex : j-1) : glm::max(0, j-1);
                int idx1 = j;
                int idx2 = e.loop ? ((j+1) % int(e.controlpoints.size())) : j+1;
                int idx3 = e.loop ? ((j+2) % int(e.controlpoints.size())) : glm::min(j+2, lastIndex);
                Math::Centripedal3 curve(e.controlpoints[idx0].transform.position,
                                         e.controlpoints[idx1].transform.position,
                                         e.controlpoints[idx2].transform.position,
                                         e.controlpoints[idx3].transform.position);
                for (int k=0; k<e.reps; k++)
                {
                    const float t = clamp(float(k)/float(e.reps), 0.0f, 1.0f);
                    glm::vec3 pos = curve.eval(t);
                    inst.eval.shape = e.controlpoints[idx1].shape;
                    inst.eval.transform.position = e.transform.orientation * (pos * e.transform.scale) + e.transform.position;
                    inst.eval.transform.orientation = e.transform.orientation * normalize(slerp(e.controlpoints[idx1].transform.orientation, e.controlpoints[idx2].transform.orientation, t));
                    inst.eval.shape.color = mix(e.controlpoints[idx1].shape.color, e.controlpoints[idx2].shape.color, t);
                    inst.eval.shape.dim = mix(e.controlpoints[idx1].shape.dim*e.controlpoints[idx1].transform.scale, e.controlpoints[idx2].shape.dim*e.controlpoints[idx2].transform.scale, t) * e.transform.scale;
                    inst.eval.shape.cornerRadius = mix(e.controlpoints[idx1].shape.cornerRadius, e.controlpoints[idx2].shape.cornerRadius, t);
                    inst.eval.shape.thickness = mix(e.controlpoints[idx1].shape.thickness, e.controlpoints[idx2].shape.thickness, t);
                    destination.push_back(inst);
                    
                    if (e.op.mirrorYZ == true)
                    {
                        if (inst.eval.shape.type == UB::Shape::Type::BEZIER)
                        {
                            inst.eval.shape.cornerRadius = 1.0f - inst.eval.shape.cornerRadius;
                        }

                        inst.eval.transform.mirrorYZ();
                        destination.push_back(inst);
                    }
                }
                
            }
        }
        else if (std::holds_alternative<EditFilter>(ee))
        {
            const EditFilter& e = std::get<EditFilter>(ee);
            if (e.op.type == EditOp::Add) hasAdd = true;
            // TODO
        }
        else if (std::holds_alternative<EditGenerator>(ee))
        {
            const EditGenerator& e = std::get<EditGenerator>(ee);
            if (e.op.type == EditOp::Add) hasAdd = true;
            // TODO
        }
    }
    
    // everything else cannot result in a surface
    return (hasAdd && destination.size() > 0);
}

float UB::evalDistance(const glm::vec3 &p, const InstructionEval &e)
{
    using namespace glm;
    const vec3 pt(glm::conjugate(e.transform.orientation)*(p - e.transform.position));
    
    if (e.shape.type == Shape::Type::SUPERPRIM)
        return sdSuperprim(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::ELLIPSOID)
        return sdEllipsoid(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.y, e.shape.cornerRadius.x);
    else if (e.shape.type == Shape::Type::BEZIER)
        return sdBezier(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::SQUISHENGON)
        return sdSquishengon(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::UBERPRIM)
        return sdUberprim(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);

    else return FLT_MAX;
}

float UB::evalDistanceLocal(const glm::vec3 &p, const InstructionEval &e)
{
    using namespace glm;
    
    if (e.shape.type == Shape::Type::SUPERPRIM)
        return sdSuperprim(p, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::ELLIPSOID)
        return sdEllipsoid(p, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.y, e.shape.cornerRadius.x);
    else if (e.shape.type == Shape::Type::BEZIER)
        return sdBezier(p, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::SQUISHENGON)
        return sdSquishengon(p, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    else if (e.shape.type == Shape::Type::UBERPRIM)
        return sdUberprim(p, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);

    else return FLT_MAX;
}

float UB::appendDistance(const float prevSd, const float sd, const int op, const float blend)
{
    using namespace glm;
    switch (op)
    {
    case EditOp::Sub:
        return blend > 0.0f ? smax(prevSd, -sd, blend) : max(prevSd, -sd);
    case EditOp::Add:
        return blend > 0.0f ? smin(prevSd, sd, blend) : min(prevSd, sd);
    case EditOp::Crop:
        return blend > 0.0f ? smax(prevSd, sd, blend) : max(prevSd, sd);
    default:
        return FLT_MAX;
    }
}

glm::range UB::evalRange(const glm::vec3& p, const float nodeHalfWidth, const InstructionEval &e)
{
    using namespace glm;
    const float pushOffEdgeEpsilon = 0.01f;
    const vec3 pt(glm::conjugate(e.transform.orientation) * (p - e.transform.position));
    
    if (e.shape.type == Shape::Type::SUPERPRIM)
    {
        const float b = pushOffEdgeEpsilon + nodeHalfWidth * UB::makeRotationBounds(e.transform.orientation);
        const range3 r(pt - b, pt + b);
        return sdSuperprim(r, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
    }
    else if (e.shape.type == Shape::Type::ELLIPSOID)
    {
        float sd  = sdEllipsoid(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.y, e.shape.cornerRadius.x);
        const float b = pushOffEdgeEpsilon + nodeHalfWidth * 1.732050807568877f;
        return range(sd - b, sd + b);
    }
    else if (e.shape.type == Shape::Type::BEZIER)
    {
        float sd  = sdBezier(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
        const float b = pushOffEdgeEpsilon + nodeHalfWidth * 1.732050807568877f;
        return range(sd - b, sd + b);
    }
    else if (e.shape.type == Shape::Type::SQUISHENGON)
    {
        float sd  = sdSquishengon(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
        const float b = pushOffEdgeEpsilon + nodeHalfWidth * 1.732050807568877f;
        return range(sd - b, sd + b);
    }
    else if (e.shape.type == Shape::Type::UBERPRIM)
    {
        float sd = sdUberprim(pt, e.shape.dim, e.shape.thickness, e.shape.cornerRadius.x, e.shape.cornerRadius.y);
        const float b = pushOffEdgeEpsilon + nodeHalfWidth * 1.732050807568877f;
        return range(sd - b, sd + b);
    }

    else return glm::range(FLT_MAX);

}

glm::range UB::appendRange(const glm::range& prevSd, const glm::range& sd, const int& op, const float& blend)
{
    using namespace glm;
    switch (op)
    {
    case EditOp::Sub:
        return blend > 0.0f ? smax(prevSd, -sd, blend) : max(prevSd, -sd);
    case EditOp::Add:
        return blend > 0.0f ? smin(prevSd, sd, blend) : min(prevSd, sd);
    case EditOp::Crop:
        return blend > 0.0f ? smax(prevSd, sd, blend) : max(prevSd, sd);
    default:
        return glm::range(FLT_MAX);
    }
}

float UB::evalDistance(const glm::vec3& p, const InstructionList& evals)
{
    // brute force distance eval over edit list
    float sd = FLT_MAX;
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        
        // TODO: implement PUSH/POP
        
        if (e.type == InstructionType::EVAL && e.eval.op.type != 2)
        {
            const float sde = evalDistance(p, e.eval);
            sd = appendDistance(sd, sde, e.eval.op.type, e.eval.op.blend);
        }
    }
    return sd;
}

glm::vec3 UB::evalGradient(const glm::vec3& p, const InstructionList& evals)
{
    const float epsi = 0.01f;
    const float sd0 = evalDistance(p + glm::vec3(epsi, 0.0f, 0.0f), evals);
    const float sd1 = evalDistance(p + glm::vec3(0.0f, epsi, 0.0f), evals);
    const float sd2 = evalDistance(p + glm::vec3(0.0f, 0.0f, epsi), evals);
    const float sd3 = evalDistance(p - glm::vec3(epsi, 0.0f, 0.0f), evals);
    const float sd4 = evalDistance(p - glm::vec3(0.0f, epsi, 0.0f), evals);
    const float sd5 = evalDistance(p - glm::vec3(0.0f, 0.0f, epsi), evals);
    return glm::vec3(sd0 - sd3, sd1 - sd4, sd2 - sd5) / (2.0f * epsi);
}

glm::vec3 UB::evalAlbedo(const glm::vec3& p, const InstructionList& evals)
{
    float sd = FLT_MAX;
    glm::vec3 color(0.0f);
    
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        
        // TODO: implement PUSH/POP
        
        if (e.type == InstructionType::EVAL)
        {
            const float sde = evalDistance(p, e.eval);
            
            if (e.eval.op.blend == 0.0f)
            {
                if (fabsf(sde) <= fabsf(sd))
                {
                    color = e.eval.shape.color;
                }
            }
            else
            {
                color = mix(glm::vec3(e.eval.shape.color), color, glm::clamp(fabsf(sde) / e.eval.op.blend, 0.0f, 1.0f));
            }
            UB::appendDistance(sd, sde, e.eval.op.type, e.eval.op.blend);
        }
    }
    return color;

}

int UB::getClosestEditIndex(const glm::vec3& p, const InstructionList& evals)
{
    const UB::Instruction* closest = nullptr;
    float sd = FLT_MAX;
    
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        
        // TODO: implement PUSH/POP
        
        if (e.type == InstructionType::EVAL)
        {
            const float sde = evalDistance(p, e.eval);

            if (e.eval.op.type == UB::EditOp::Color)
            {
                const float T = (glm::abs(sd) > glm::abs(sde) && e.eval.op.blend < FLT_EPSILON) ? 0.0f : glm::max(glm::min((sde) / glm::max(FLT_EPSILON, e.eval.op.blend), 1.0f), 0.0f);
                if (T < 0.5f)
                {
                    closest = &e;
                }
            }
            else
            {
                const float eps = 0.005f;
                const float prev = sd;
                sd = UB::appendDistance(sd, sde, e.eval.op.type, e.eval.op.blend);
                const float T = 1.0f - glm::abs(glm::max(glm::min((glm::abs(prev) - glm::abs(sd)) / (glm::abs(prev) + glm::abs(sde) - 2.0f * glm::abs(sd)), 1.0f), -1.0f));
                if (T < 0.5f)
                {
                    closest = &e;
                }
            }
        }
    }
    if (closest != nullptr)
    {
        return closest->common.editIndex;
    }
    else
    {
        return -1;
    }
}

bool UB::findClosestIntersection(const glm::vec3& pos, const glm::vec3 &dir, glm::vec3& pp, const InstructionList& evals)
{
    // brute force sphere-tracing
    pp = pos;
    float sdPrev = FLT_MAX;
    for (int i = 0; i < 32; i++)
    {
        float sd = UB::evalDistance(pp, evals);
        if (fabsf(sd) < 0.1f)
        {
            return true;
        }
        pp += dir * sd * 0.5f;
    }
    return false;
}

float UB::evalDistanceAndGradient(const glm::vec3& p, glm::vec3 &grad, const InstructionList& evals, float epsi)
{
    float sd  = FLT_MAX;
    float sd0 = FLT_MAX;
    float sd1 = FLT_MAX;
    float sd2 = FLT_MAX;
    float sd3 = FLT_MAX;
    float sd4 = FLT_MAX;
    float sd5 = FLT_MAX;
    
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        
        // TODO: implement PUSH/POP
        
        if (e.type == InstructionType::EVAL && e.eval.op.type != 2)
        {
            const float sde  = evalDistance(p, e.eval);
            const float sde0 = evalDistance(p + glm::vec3(epsi, 0.0f, 0.0f), e.eval);
            const float sde1 = evalDistance(p + glm::vec3(0.0f, epsi, 0.0f), e.eval);
            const float sde2 = evalDistance(p + glm::vec3(0.0f, 0.0f, epsi), e.eval);
            const float sde3 = evalDistance(p - glm::vec3(epsi, 0.0f, 0.0f), e.eval);
            const float sde4 = evalDistance(p - glm::vec3(0.0f, epsi, 0.0f), e.eval);
            const float sde5 = evalDistance(p - glm::vec3(0.0f, 0.0f, epsi), e.eval);
            
            sd  = appendDistance(sd,  sde,  e.eval.op.type, e.eval.op.blend);
            sd0 = appendDistance(sd0, sde0, e.eval.op.type, e.eval.op.blend);
            sd1 = appendDistance(sd1, sde1, e.eval.op.type, e.eval.op.blend);
            sd2 = appendDistance(sd2, sde2, e.eval.op.type, e.eval.op.blend);
            sd3 = appendDistance(sd3, sde3, e.eval.op.type, e.eval.op.blend);
            sd4 = appendDistance(sd4, sde4, e.eval.op.type, e.eval.op.blend);
            sd5 = appendDistance(sd5, sde5, e.eval.op.type, e.eval.op.blend);
        }
    }
    
    grad = glm::vec3(sd0 - sd3,
                     sd1 - sd4,
                     sd2 - sd5) / (2.0f*epsi);
    return sd;
}


bool UB::findClosestIntersection(const glm::vec3& pos, const glm::vec3 &dir, glm::vec3& pp, float tmin, float tmax, const InstructionList& evals)
{
    pp = pos + dir * tmin;
    const float epsilon = 0.1f;
    float travel = tmin;
    float last_sd;
    float sd = FLT_MAX;
    glm::vec3 grad;
    
    for (int i = 0; i < 32; i++)
    {
        if (travel >= tmax) { return false; }
        last_sd = sd;
        sd = evalDistanceAndGradient(pp, grad, evals, epsilon);
        sd /= FLENGTH(grad);
        
        if (fabsf(sd) < epsilon)
        {
            return true;
        }
        
        if ((glm::sign(last_sd) == glm::sign(sd)))
        {
            travel += sd;
            pp = pp + dir * sd;
        }
        else
        {
            glm::vec3 pp0 = pp + dir * sd;
            travel += FLENGTH(pp - pp0)*0.5*glm::sign(sd);
            pp = (pp + pp0)*0.5f;
        }
    }
    return false;
}

glm::range3 UB::getAABB(const InstructionEval& inst)
{
    glm::range3 aabb(glm::vec3(1e36f), glm::vec3(-1e36f));

    for (int x = -1; x <= 1; x += 2)
    {
        for (int y = -1; y <= 1; y += 2)
        {
            for (int z = -1; z <= 1; z += 2)
            {
                if (inst.type == InstructionType::EVAL)
                {
                    const glm::vec3 p(inst.transform.orientation*((inst.op.blend+inst.shape.dim)*glm::vec3(x, y, z)) + inst.transform.position);
                    aabb.x.lb = glm::min(aabb.x.lb, p.x);
                    aabb.y.lb = glm::min(aabb.y.lb, p.y);
                    aabb.z.lb = glm::min(aabb.z.lb, p.z);
                    aabb.x.ub = glm::max(aabb.x.ub, p.x);
                    aabb.y.ub = glm::max(aabb.y.ub, p.y);
                    aabb.z.ub = glm::max(aabb.z.ub, p.z);
                }
            }
        }
    }
    return aabb;
}


glm::range3 UB::getAABB(const InstructionList& evals)
{
    glm::range3 totalAABB(1e32f, -1e32f);
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        if (e.type == InstructionType::EVAL)
        {
            const UB::InstructionEval& ev = e.eval;
            if (ev.op.type == EditOp::Color)
            {
                // exclude color ops
                continue;
            }

            const glm::range3 r = getAABB(ev);
            if (ev.op.type == EditOp::Crop)
            {
                totalAABB.x.lb = r.x.lb;
                totalAABB.y.lb = r.y.lb;
                totalAABB.z.lb = r.z.lb;
                totalAABB.x.ub = r.x.ub;
                totalAABB.y.ub = r.y.ub;
                totalAABB.z.ub = r.z.ub;
            }
            else
            {
                totalAABB.x.lb = glm::min(totalAABB.x.lb, r.x.lb);
                totalAABB.y.lb = glm::min(totalAABB.y.lb, r.y.lb);
                totalAABB.z.lb = glm::min(totalAABB.z.lb, r.z.lb);
                totalAABB.x.ub = glm::max(totalAABB.x.ub, r.x.ub);
                totalAABB.y.ub = glm::max(totalAABB.y.ub, r.y.ub);
                totalAABB.z.ub = glm::max(totalAABB.z.ub, r.z.ub);
            }
        }
    }
    return totalAABB;
}

glm::range3 UB::getTightAABB(const InstructionList& evals)
{
    glm::range3 totalAABB(1e32f, -1e32f);
    for (int i = 0; i < evals.size(); i++)
    {
        const UB::Instruction& e = evals[i];
        if (e.type == InstructionType::EVAL)
        {
            const UB::InstructionEval& ev = e.eval;
            if (ev.op.type == EditOp::Sub or ev.op.type == EditOp::Color)
            {
                // exclude color ops
                continue;
            }

            const glm::range3 r = getAABB(ev);
            if (ev.op.type == EditOp::Crop)
            {
                totalAABB.x.lb = r.x.lb;
                totalAABB.y.lb = r.y.lb;
                totalAABB.z.lb = r.z.lb;
                totalAABB.x.ub = r.x.ub;
                totalAABB.y.ub = r.y.ub;
                totalAABB.z.ub = r.z.ub;
            }
            else
            {
                totalAABB.x.lb = glm::min(totalAABB.x.lb, r.x.lb);
                totalAABB.y.lb = glm::min(totalAABB.y.lb, r.y.lb);
                totalAABB.z.lb = glm::min(totalAABB.z.lb, r.z.lb);
                totalAABB.x.ub = glm::max(totalAABB.x.ub, r.x.ub);
                totalAABB.y.ub = glm::max(totalAABB.y.ub, r.y.ub);
                totalAABB.z.ub = glm::max(totalAABB.z.ub, r.z.ub);
            }
        }
    }
    return totalAABB;
}


uint64_t UB::getUUID(const InstructionList& evals)
{
    XXH64_state_t *const state = XXH64_createState();
    XXH64_hash_t const seed = 12345678;
    if (XXH64_reset(state, seed) == XXH_ERROR)
    {
        printf("XXHash error");
    }
    XXH64_update(state, evals.data(), evals.size() * sizeof(UB::Instruction));
    uint64_t listHash = XXH64_digest(state);
    XXH64_freeState(state);
    return listHash;
}


