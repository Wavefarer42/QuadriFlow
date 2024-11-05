#pragma once

#include "Core/Id.h"
#include "Math/Transform.h"
#include "UB/UB.h"
#include <variant>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <uuid.h>
#include <memory>

namespace UB
{
    using AnimProperty = std::variant<float, int, bool, glm::vec2, glm::vec3, glm::vec4, glm::quat>;
    using AnimPropertyPtr = std::variant<float*, int*, bool*, glm::vec2*, glm::vec3*, glm::vec4*, glm::quat*>;
    
    enum class AnimInterpolant
    {
        STEP = 0,
        LINEAR = 1,
        CUBIC = 2
    };

    template <typename P>
    struct AnimFrame
    {
        AnimInterpolant interpolant;
        uint32_t timeMs;
        P value;
    };

    struct AnimTrackBase
    {
        virtual ~AnimTrackBase() {}
        virtual double duration() = 0;
        virtual AnimProperty getValue(float time) = 0;
        virtual bool hasFrame(uint32_t timeMs) = 0;
        virtual void sort() = 0;
        virtual AnimTrackBase* cloneEmpty() = 0;
        virtual AnimTrackBase* clone() = 0;
        virtual void removeFrame(uint32_t timeMs) = 0;
        virtual AnimInterpolant getFrameInterpolant(uint32_t timeMs) = 0;
        virtual void setFrameInterpolant(uint32_t timeMs, AnimInterpolant interpolant) = 0;
        virtual bool hasFrames() = 0;

        template <typename P>
        bool is() const;

        virtual bool copyFrameFrom(AnimTrackBase* other, uint32_t timeMs) = 0;

        template <typename P>
        bool setFrameValue(uint32_t timeMs, const P& value);

        template <typename P>
        P getFrameValue(uint32_t timeMs);
        
        template <typename P>
        P getValue(float time)
        {
            AnimProperty value = getValue(time);
            if (std::holds_alternative<P>(value))
            {
                return std::get<P>(value);
            }
            else
            {
                return P();
            }
        }
    };

    using AnimTrackPtr = std::shared_ptr<AnimTrackBase>;

    struct AnimBinding
    {
        Id entityId;
        AnimTrackPtr track;
        AnimPropertyPtr property;
    };

    template <typename P>
    struct AnimTrack : public AnimTrackBase
    {
        std::vector<AnimFrame<P>> frames;
        double duration()
        {
            sort();

            if (frames.size() == 0)
            {
                return 0.0;
            }

            return frames.back().timeMs * 0.001;
        }

        void sort()
        {
            std::sort(frames.begin(), frames.end(), [](const AnimFrame<P>& a, const AnimFrame<P>& b)
            {
                return a.timeMs < b.timeMs;
            });
        }

        AnimTrackBase* cloneEmpty()
        {
            return (AnimTrackBase*)(new AnimTrack<P>());
        }

        AnimTrackBase* clone()
        {
            auto* newTrack = new AnimTrack<P>();
            std::copy(frames.begin(), frames.end(), back_inserter(newTrack->frames));
            return (AnimTrackBase*)(newTrack);
        }

        bool copyFrameFrom(AnimTrackBase* other, uint32_t timeMs)
        {
            if (other->is<P>())
            {
                const AnimTrack<P>* src = dynamic_cast<const AnimTrack<P>*>(other);
                auto srcIt = std::find_if(src->frames.begin(), src->frames.end(), [&timeMs](const AnimFrame<P>& x) {
                    return x.timeMs == timeMs;
                });
                if (srcIt != src->frames.end())
                {
                    auto dstIt = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                        return x.timeMs == timeMs;
                    });
                    if (dstIt != frames.end())
                    {
                        *dstIt = *srcIt;
                    }
                    else
                    {
                        frames.push_back(*srcIt);
                    }
                    return true;
                }
            }
            return false;
        }

        void setFrameValue(uint32_t timeMs, const AnimProperty& value)
        {
            assert(std::holds_alternative<P>(value));
            if (std::holds_alternative<P>(value))
            {
                auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                    return x.timeMs == timeMs;
                });
                if (it != frames.end())
                {
                }
                else
                {
                    it->timeMs = timeMs;
                    it->value = std::get<P>(value);
                }
            }
        }

        AnimProperty getFrameValue(uint32_t timeMs)
        {
            if (frames.size() == 0)
                return P();
            
            auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (it != frames.end())
            {
                return it->value;
            }
            else return P();
        }

        AnimProperty getValue(float timeSec)
        {
            if (frames.size() == 0)
                return P();
            
            const uint32_t timeMs = timeSec*1000.0f;
            
            if (timeMs <= frames[0].timeMs)
                return frames[0].value;

            if (timeMs >= frames.back().timeMs)
                return frames.back().value;

            const AnimFrame<P>* frame0 = nullptr;
            const AnimFrame<P>* frame1 = nullptr;
            const AnimFrame<P>* frame2 = nullptr;
            const AnimFrame<P>* frame3 = nullptr;
            assert(frames.size() > 1);

            for (size_t i=0; i<frames.size() - 1; i++)
            {
                if (timeMs >= frames[i].timeMs && timeMs <= frames[i+1].timeMs)
                {
                    if ((int)(i) - 1 >= 0)
                    {
                        frame0 = &(frames[i - 1]);
                    }
                    else
                    {
                        frame0 = &(frames[i]);
                    }
                    
                    frame1 = &(frames[i]);
                    frame2 = &(frames[i + 1]);

                    if (i + 2 < frames.size())
                    {
                        frame3 = &(frames[i + 2]);
                    }
                    else
                    {
                        frame3 = &(frames[i + 1]);
                    }
                    
                    break;
                }
            }
            assert(frame0 != nullptr);
            assert(frame1 != nullptr);
            assert(frame2 != nullptr);
            assert(frame3 != nullptr);

             if (frame1->interpolant == AnimInterpolant::STEP)
            {
                if (timeMs >= frame2->timeMs)
                {
                    return frame2->value;
                }
                else
                {
                    return frame1->value;  
                }
            }

            const float duration = frame2->timeMs - frame1->timeMs;
            const float a = (timeMs - frame1->timeMs) / duration;
            float t = 0.0f;
            
            if (frame1->interpolant == AnimInterpolant::LINEAR)
            {
                t = a;
            }
            else if (frame1->interpolant == AnimInterpolant::CUBIC)
            {
                // actually a quadratic ease in-out function, cubic might be too strong on accel decel
                t = glm::clamp(0.5f * (sin((glm::clamp(a, 0.0f, 1.0f) - 0.5f) * glm::pi<float>()) + 1.0f), 0.0f, 1.0f);
            }
            
            if constexpr (std::is_same_v<P, glm::quat>)
            {
                return cubicInterpQuat(frame0->value, frame1->value, frame2->value, frame3->value, t);
            }
            else
            {
                return cubicInterp(frame0->value, frame1->value, frame2->value, frame3->value, t);
            }
        }

        bool hasFrame(uint32_t timeMs)
        {
            auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });

            return it != frames.end();
        }

        void removeFrame(uint32_t timeMs)
        {
            auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (it != frames.end())
            {
                frames.erase(it);
            }
        }

        void setFrameInterpolant(uint32_t timeMs, AnimInterpolant interpolant)
        {
            auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (it != frames.end())
            {
                it->interpolant = interpolant;
            }
        }

        AnimInterpolant getFrameInterpolant(uint32_t timeMs)
        {
            auto it = std::find_if(frames.begin(), frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (it != frames.end())
            {
                return it->interpolant;
            }
            return AnimInterpolant::STEP; // need an invalid enum entry?
        }

        bool hasFrames()
        {
            if (frames.size() > 0)
            {
                return true;
            }

            return false;
        }

        private:

        P cubicInterp(const P& p0, const P& p1, const P& p2, const P& p3, float t)
        {
            P a = 2.0f * p1;
            P b = p2 - p0;
            P c = 2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3;
            P d = -p0 + 3.0f * p1 - 3.0f * p2 + p3;

            return 0.5f * (a + (b * t) + (c * t * t) + (d * t * t * t));
        }

        glm::vec3 quat_log(const glm::quat &q, float eps=1e-8f)
        {
            float length = sqrtf(q.x*q.x + q.y*q.y + q.z*q.z);
            
            if (length < eps)
            {
                return glm::vec3(q.x, q.y, q.z);
            }
            else
            {
                float halfangle = acosf(glm::clamp(q.w, -1.0f, 1.0f));
                return halfangle * (glm::vec3(q.x, q.y, q.z) / length);
            }
        }

        glm::vec3 quat_to_scaled_angle_axis(const glm::quat &q, float eps=1e-8f)
        {
            return 2.0f * quat_log(q, eps);
        }

        glm::quat quat_exp(const glm::vec3 &v, float eps=1e-8f)
        {
            float halfangle = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
            
            if (halfangle < eps)
            {
                return glm::normalize(glm::quat(1.0f, v.x, v.y, v.z));
            }
            else
            {
                float c = cosf(halfangle);
                float s = sinf(halfangle) / halfangle;
                return glm::quat(c, s * v.x, s * v.y, s * v.z);
            }
        }

        glm::quat quat_from_scaled_angle_axis(const glm::vec3 &v, float eps=1e-8f)
        {
            return quat_exp(v * 0.5f, eps);
        }

        glm::quat quat_abs(const glm::quat &x)
        {
            return x.w < 0.0f ? -x : x;
        }

        glm::quat quat_mul_inv(const glm::quat &q, const glm::quat &p)
        {
            return q * glm::quat(-p.w, p.x, p.y, p.z);
        }

        P cubicInterpQuat(const P &r0, const P &r1, const P &r2, const P &r3, float t){
            // simply using cubicInterp on quaternions as if they were vec4 won't get us shortest path interpolation.
            // hence this guy's solution:
            //from https://theorangeduck.com/page/cubic-interpolation-quaternions
            //math repo https://github.com/orangeduck/QuaternionAverage/blob/637be07622066469662225a4f6fdee92749a110c/quat.h

            glm::vec3 r1_sub_r0 = quat_to_scaled_angle_axis(quat_abs(quat_mul_inv(r1, r0)));
            glm::vec3 r2_sub_r1 = quat_to_scaled_angle_axis(quat_abs(quat_mul_inv(r2, r1)));
            glm::vec3 r3_sub_r2 = quat_to_scaled_angle_axis(quat_abs(quat_mul_inv(r3, r2)));

            glm::vec3 v1 = (r1_sub_r0 + r2_sub_r1) * 0.5f;
            glm::vec3 v2 = (r2_sub_r1 + r3_sub_r2) * 0.5f;

            // hermite
            float w1 = 3*t*t - 2*t*t*t;
            float w2 = t*t*t - 2*t*t + t;
            float w3 = t*t*t - t*t;
            
            float q1 = 6*t - 6*t*t;
            float q2 = 3*t*t - 4*t + 1;
            float q3 = 3*t*t - 2*t;
            
            glm::vec3 r1_sub_r0_2 = quat_to_scaled_angle_axis(quat_abs(quat_mul_inv(r2, r1)));   
            
            glm::quat rot = quat_from_scaled_angle_axis(w1*r1_sub_r0_2 + w2*v1 + w3*v2) * r1;

            return rot;
        }
    };

    template <typename P>
    bool AnimTrackBase::is() const
    {
        return dynamic_cast<const AnimTrack<P>*>(this);
    }

    template <typename P>
    bool AnimTrackBase::setFrameValue(uint32_t timeMs, const P& value)
    {
        if (is<P>())
        {
            AnimTrack<P>* track = dynamic_cast<AnimTrack<P>*>(this);
            auto frame = std::find_if(track->frames.begin(), track->frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (frame != track->frames.end())
            {
                frame->value = value;
                return true;
            }
            else
            {
                track->frames.push_back(AnimFrame<P>{AnimInterpolant::LINEAR, timeMs, value});
                return true;
            }
        }
        return false;
    }

    template <typename P>
    P AnimTrackBase::getFrameValue(uint32_t timeMs)
    {
        if (is<P>())
        {
            AnimTrack<P>* track = dynamic_cast<AnimTrack<P>*>(this);
            auto frame = std::find_if(track->frames.begin(), track->frames.end(), [&timeMs](const AnimFrame<P>& x) {
                return x.timeMs == timeMs;
            });
            if (frame != track->frames.end())
            {
                return frame->value;
            }
        }
        return P();
    }
    
    struct AnimData
    {
        uuids::uuid uuid;
        std::string name;
        std::unordered_map<std::string, AnimTrackPtr> tracks;

        AnimData(const AnimData& other) = default;
        AnimData()
        {
            uuid = UB::generateUUID();
        }
        
        template <typename P>
        AnimTrackPtr getCreateTrack(const std::string& targetPath)
        {
            auto it = tracks.find(targetPath);
            if (it == tracks.end())
            {
                AnimTrackPtr track = std::make_shared<AnimTrack<P>>();
                tracks.emplace(targetPath, track);
                return track;
            }
            else
            {
                if (it->second->is<P>())
                {
                    return it->second;
                }
                else
                {
                    return nullptr;
                }
            }
        }

        bool trackHasFrame(const std::string& targetPath, uint32_t timeMs) const
        {
            const auto it = tracks.find(targetPath);
            if (it == tracks.end())
                return false;

            return it->second->hasFrame(timeMs);
        }

        double calculateDuration() const
        {
            double duration = 0.0;
            for (const auto& t : tracks)
            {
                duration = glm::max(duration, t.second->duration());
            }
            return duration;
        }

        void copyFrameFrom(const AnimData& src, uint32_t timeMs)
        {
            for (const auto& srcTrack : src.tracks)
            {
                auto it = tracks.find(srcTrack.first);
                if (it == tracks.end())
                {
                    tracks[srcTrack.first].reset(srcTrack.second->cloneEmpty());
                }
                tracks[srcTrack.first]->copyFrameFrom(srcTrack.second.get(), timeMs);
            }
        }

        void copyFrom(const AnimData& src)
        {
            name = src.name;
            uuid = src.uuid;
            tracks.clear();

            for (const auto& srcTrack : src.tracks)
            {
                tracks[srcTrack.first].reset(srcTrack.second->clone());
            }
        }

        std::shared_ptr<AnimData> clone() const
        {
            auto ret = std::make_shared<AnimData>();
            ret->name = name;
            ret->uuid = uuid;
            for (const auto& srcTrack : tracks)
            {
                ret->tracks[srcTrack.first].reset(srcTrack.second->clone());
            }
            return ret;
        }

        void removeFrame(const std::string& targetPath, uint32_t timeMs)
        {
            const auto it = tracks.find(targetPath);
            if (it == tracks.end())
                return;

            it->second->removeFrame(timeMs);
        }

        AnimInterpolant getFrameInterpolant(const std::string& targetPath, uint32_t timeMs)
        {
            const auto it = tracks.find(targetPath);
            if (it == tracks.end())
                return AnimInterpolant::STEP;

            return it->second->getFrameInterpolant(timeMs);
        }

        void setFrameInterpolant(const std::string& targetPath, uint32_t timeMs, AnimInterpolant interpolant)
        {
            const auto it = tracks.find(targetPath);
            if (it != tracks.end())
                it->second->setFrameInterpolant(timeMs, interpolant);
        }
    };
    
    using AnimDataPtr = std::shared_ptr<AnimData>;

   
}
