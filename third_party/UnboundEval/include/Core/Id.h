#pragma once

#pragma once
#include <cstdint>
#include <functional>
#include <string>

struct Id
{
    union
    {
        struct
        {
            uint32_t index   : 24;
            uint32_t version : 7;
            uint32_t valid   : 1;
        };
        uint32_t bits;
    };
    static const uint32_t versionMask = 0xfe;
    static const uint32_t maxIndex = 0xFFFFFF;
    
	inline bool isValid() const { return (valid > 0); }
	Id() :version(0), index(0), valid(0) {};
    Id(const Id& other) : bits(other.bits) {};
    Id(const uint32_t bits) : bits(bits) {};
	bool operator==(const Id& b) const { return (bits == b.bits); }
	bool operator!=(const Id& b) const { return (bits != b.bits); }
    std::string toString() const { return std::to_string(valid) + ":" + std::to_string(version) + ":" + std::to_string(index); }
};

namespace std
{
    template<> struct hash<Id>
    {
        size_t operator()(const Id& x) const
        {
            return hash<uint32_t>()(x.bits);
        }
    };
}
