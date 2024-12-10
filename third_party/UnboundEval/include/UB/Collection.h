
#pragma once
#include "Id.h"
#include "Edit.h"
#include "Instruction.h"
#include "Animation.h"
#include <map>
#include <memory>
#include <optional>

namespace UB
{
    struct SplatParams
    {
        uint32_t type = 255;
        float alpha = 1.0f;
        float density = 1.0f;
        float scale = 0.1f;
        float faceCamera = 0.0f;
        float jitterPosition = 0.0f;
        float jitterOrientation = 0.0f;
        float jitterScale = 0.0f;
        float jitterBump = 0.0f;
        float jitterHue = 0.0f;
        float jitterBrightness = 0.0f;
        float jitterMetalness = 0.0f;
        float jitterRoughness = 0.0f;
    };
    
    struct Material
    {
        float metalness;
        float roughness;
        float emission = 0.0f;
        glm::vec4 colorize = glm::vec4(0.0f);
        float unlit = 0.0f;
        float looseness = 0.0f;
        std::optional<SplatParams> splat;
    };
    
    struct Model
    {
        EditList editList;
        Material material;
        uint32_t maxRefineLevel;
        float shapeDensity = 1000.0f;
        float minBounds = 0.0f;
    };
    
    using ModelRef = uuids::uuid;
    
    struct Panel
    {
        // TODO: add members
        glm::vec3 dim;
    };
    
    struct Camera
    {
        float nearPlane = 0.005f;
        float fov = glm::radians(60.0f);
        float exposure = 1.0f;
        bool enableDOF = false;
        float focalDistance = 1.0f;
        float focalWidth = FLT_MAX;
        
        // Bloom
        float bloomIntensity = 0.08f;
        float bloomRadius = 0.005f;
        float bloomThreshold = 0.0f;
        
        // Image adjustments
        float adjustContrast = 0.0f;
        float adjustVibrance = 0.0f;
        float adjustSaturation = 0.0f;
        float adjustBlacks = 0.0f;
        
        // Vignette
        float vignetteMidPoint = 0.5f;
        float vignetteRoundness = 0.5f;
        float vignetteFeather = 0.5f;
        glm::vec4 vignetteColor = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f); // alpha 0 or 1 to turn effect off/on
    };
    
    struct Light
    {
        enum Type
        {
            DIRECTIONAL = 0,
            POINT       = 1,
            SPOT        = 2
        };

        Type type = DIRECTIONAL;
        bool castShadow = false;
        glm::vec3 color = glm::vec3(1.0f);
        float intensity = 1.0f;
        glm::vec2 dim = glm::vec2(1.0f);
        bool render = true;
        float penumbra = 1.0f;
        float radius = 1.0f;
    };
    
    struct Environment
    {
        struct GradientParams
        {
            // three colors with B at the horizon
            glm::vec4 A = glm::vec4(218.0f / 255.0f, 217.0f / 255.0f, 255.0f / 255.0f, 4.0f);
            glm::vec4 B = glm::vec4(218.0f / 255.0f, 217.0f / 255.0f, 255.0f / 255.0f, 5.0f);
            glm::vec4 C = glm::vec4(218.0f / 255.0f, 217.0f / 255.0f, 255.0f / 255.0f, 8.0f);
            float Q = 1.0f;
        };
        
        enum class Type : uint32_t
        {
            DEFAULT = 0,
            SKY     = 1
        };
        
        struct Fog
        {
            glm::vec4 color = glm::vec4(0.5f, 0.5f, 0.5f, 1.0f);
            bool colorFromIbl = false;
            float density = 0.0f;
            float height = 0.0f;
            float heightFalloff = 1.0f;
            float start = 0.0f;
            float maxOpacity = 1.0f;
            float inscatteringStart = 0.0f;
            float inscatteringSize = -1.0f;
        };

        Type type = Type::DEFAULT;
        GradientParams gradient;
        Fog fog;
        float ssaoRadius = 1.0f;
        float ssaoStrength = 1.0f;
        float emissiveStrength = 255.0f;
        float intensity = 1.0f;
    };
    
    enum class PhysicsType
    {
        NONE = 0,
        STATIC,
        KINEMATIC,
        DYNAMIC
    };

    struct PhysicsMaterial
    {
        float restitution = 0.0f;
        float friction = 0.2f;
        float linearDamping = 0.05f;
        float angularDamping = 0.05f;
        float gravityFactor = 1.0f;
    };
    
    struct PhysicsBody
    {
        PhysicsType type = PhysicsType::NONE;
        bool accurateCollision = false;
        uint8_t layer = 0;
        PhysicsMaterial material;
    };
        
    struct Constraint
    {
        uuids::uuid nodeA;
        uuids::uuid nodeB;
        
        enum Type
        {
            INVALID = 0,
            FIXED,
            DISTANCE,
            HINGE,
            POINT,
            CONE,
            SLIDER
        };

        static const std::string NAME_INVALID;
        static const std::string NAME_FIXED;
        static const std::string NAME_DISTANCE;
        static const std::string NAME_HINGE;
        static const std::string NAME_POINT;
        static const std::string NAME_CONE;
        static const std::string NAME_SLIDER;
        static const std::string& getTypeName(Type t);
        static const Type getType(const std::string& typeName);

        struct Settings
        {
            bool enabled = true;
            float frequency = 0.0f;
            float damping = 0.0f;
            float limitsMin = -glm::pi<float>();
            float limitsMax = glm::pi<float>();
            float maxFriction = 0.0f;
            float halfAngle = 0.0f;
            std::array<glm::vec3, 2> offset = { glm::vec3(0) };
        };
        
        Type type;
        Settings settings;
    };
    
    struct Logic
    {
        enum Type
        {
            TYPE_SCRIPT         = 0,
            TYPE_ANIMATION      = 1,
            TYPE_NATIVE_TEST    = 2 // Example only
        };
        Type type;
        uuids::uuid dataRef;
        std::map<std::string, std::pair<uuids::uuid, std::string>> inputs;
    };
    
    using NodeData = std::variant<std::monostate, ModelRef, Panel, Camera, Light, Environment, Constraint, Logic>;
    using Property = std::variant<std::monostate, Id, uuids::uuid, double, bool, std::string, glm::vec2, glm::vec3, glm::vec4, glm::quat, Math::Transform>;
    using PropertyList = std::map<std::string, Property>;
    using Tags = std::vector<std::string>;
    
    struct Node
    {
        uuids::uuid id;
        std::string name;
        bool hidden;
        bool overlay;
        bool shadowCaster;
        bool shadowReceiver;
        Tags tags;
        Math::Transform transform;
        PhysicsBody physicsBody;
        std::vector<Node> children;
        NodeData data;
        PropertyList properties;
    };
    
    struct NodeFilter
    {
        enum Condition
        {
            INVALID,
            IS_TYPE,
            IS_NAME,
            IS_UUID,
            HAS_PROPERTY,
            HAS_TAG,
        };
        Condition condition;
        std::string value;
        NodeFilter() = default;
        NodeFilter(const NodeFilter& other) = default;
        explicit NodeFilter(const std::string& input);
        bool check(const Node& node) const;
        using List = std::vector<NodeFilter>;
        static bool check(const Node& node, const List& filterList);
    };
    
    struct Script
    {
        std::string name;
        std::string code;
    };
    
    struct Collection
    {
        std::map<uuids::uuid, Script> scripts;
        std::map<ModelRef, Model> models;
        std::vector<Node> nodes;
        // std::map<uuids::uuid, AnimDataPtr> animations;
    };
    
    using CollectionPtr = std::shared_ptr<Collection>;
    void visitRecursive(UB::CollectionPtr collection, std::function<void(UB::Node*, UB::Node*)> cb);
    void visitNodeRecursive(UB::Node* node, UB::Node* parent, std::function<void(UB::Node*, UB::Node*)> cb);

}
