//
// Unbound - CommonFunc.h
// Copyright (c) 2016 Unbound Technologies, Inc. All rights reserved.
// 

#pragma once
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <string>
#include "CSG/IA.h"
#include <type_traits>
#include <variant>

namespace Util
{
	glm::ivec3 saveNextMultipleOfN(const glm::ivec3& a, const int n);
	int saveNextMultipleOfN(const int a, const int n);
	glm::ivec3 nextMultipleOfN(const glm::ivec3& a, const int n);
	glm::ivec3 nextMultipleOfN(const glm::ivec3& a, const glm::ivec3& n);
	glm::ivec2 nextMultipleOfN(const glm::ivec2& a, const int n);
	int nextMultipleOfN(const int a, const int n);
	unsigned int roundUpToNextPowerOfTwo(unsigned int x);
	size_t roundUpToNextPowerOfTwo(size_t x);
	glm::vec3 hsv2rgb(float h, float s, float v);
	glm::vec3 rgb2hsv(float r, float g, float b);
	glm::vec3 spherical(float theta, float phi);
	float angularDifference(const glm::quat& a, const glm::quat& b);
	glm::quat lookAtOrientation(glm::vec3 dirNorm, glm::vec3 desiredUp);
	bool startsWith(const std::string& str, const std::string& prefix);
	bool strContains(const std::string& str, const std::string& to_find);
	void ltrim(std::string &s);
	void rtrim(std::string &s);
	void trim(std::string &s);
	glm::vec2 intersectRayAABB(const glm::vec3& p, const glm::vec3& d, const glm::vec3& lb, const glm::vec3& ub);
    glm::vec2 intersectRayAABB(const glm::vec3& p, const glm::vec3& d, const glm::range3& aabb);
    glm::vec2 intersectRayAABBInvDist(const glm::vec3& p, const glm::vec3& inv_d, const glm::vec3& lb, const glm::vec3& ub);
	float distanceToLine(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b);
	float mad(float A, float B, float C);
	float madfrac(const float &A, const float &B);
	glm::vec3 sphericalFibonacci(const float &i, const float &n);
	unsigned int hash(unsigned int x);
	float rnd(unsigned int& rndSeed);
	glm::vec3 getRndColor(unsigned int k);
    
    template <typename T, template <typename...> class C, typename ... Ts>
    constexpr auto isTypeInList (C<Ts...> const &) -> std::disjunction<std::is_same<T, Ts>...>;

    template <typename T, typename V>
    static constexpr bool isTypeInList_v = decltype(isTypeInList<T>(std::declval<V>()))::value;

	constexpr size_t combineHash(size_t a, size_t b) { return (a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2))); }
	template<typename key_T> inline size_t hash(const key_T& key) { return std::hash<key_T>{}(key); }
	template<typename key_T, typename... args_T> inline size_t hash(const key_T& key, args_T... args) { return combineHash(hash(key), hash(args...)); }


	template<class T, class TypeList>
	struct isContainedIn;

	template<class T, class... Ts>
	struct isContainedIn<T, std::variant<Ts...>> : std::bool_constant<(... || std::is_same<T, Ts>{})> {};

} // namespace Util
