//
//  Serialization.h
//  Framework
//
//  Created by Florian Hoenig on 5/29/22.
//

#pragma once

#include "UB.h"
#include "Edit.h"
#include "Collection.h"
#include <nlohmann/json.hpp>

namespace uuids {
    static void to_json(nlohmann::json &j, const uuids::uuid &v) { j = uuids::to_string(v); }

    static void from_json(const nlohmann::json &j, uuids::uuid &v) {
        std::string s;
        if (j.is_string()) {
            j.get_to(s);
            auto val = uuids::uuid::from_string(s);
            if (val.has_value()) {
                v = *val;
                return;
            }
        } else if (j.is_number()) {
            std::mt19937 generator(j.get<unsigned int>());
            uuids::uuid_random_generator gen{generator};
            v = gen();
        }
    }
}

namespace glm {
    static void to_json(nlohmann::json &j, const u8vec3 &v) { j = nlohmann::json{v.x, v.y, v.z}; }

    static void from_json(const nlohmann::json &j, u8vec3 &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
        j.at(2).get_to(v.z);
    }

    static void to_json(nlohmann::json &j, const ivec3 &v) { j = nlohmann::json{v.x, v.y, v.z}; }

    static void from_json(const nlohmann::json &j, ivec3 &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
        j.at(2).get_to(v.z);
    }

    static void to_json(nlohmann::json &j, const vec2 &v) { j = nlohmann::json{v.x, v.y}; }

    static void from_json(const nlohmann::json &j, vec2 &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
    }

    static void to_json(nlohmann::json &j, const vec3 &v) { j = nlohmann::json{v.x, v.y, v.z}; }

    static void from_json(const nlohmann::json &j, vec3 &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
        j.at(2).get_to(v.z);
    }

    static void to_json(nlohmann::json &j, const vec4 &v) { j = nlohmann::json{v.x, v.y, v.z, v.w}; }

    static void from_json(const nlohmann::json &j, vec4 &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
        j.at(2).get_to(v.z);
        j.at(3).get_to(v.w);
    }

    static void to_json(nlohmann::json &j, const quat &v) { j = nlohmann::json{v.x, v.y, v.z, v.w}; }

    static void from_json(const nlohmann::json &j, quat &v) {
        j.at(0).get_to(v.x);
        j.at(1).get_to(v.y);
        j.at(2).get_to(v.z);
        j.at(3).get_to(v.w);
    }
}

namespace Math {
    inline void to_json(nlohmann::json &j, const Transform &t) {
        j["position"] = t.position;
        j["orientation"] = t.orientation;
        j["scale"] = t.scale;
    }

    inline void from_json(const nlohmann::json &j, Transform &t) {
        j.at("position").get_to(t.position);
        j.at("orientation").get_to(t.orientation);
        if (j.at("scale").is_number()) {
            t.scale = glm::vec3(j.at("scale").get<float>());
        } else {
            j.at("scale").get_to(t.scale);
        }
    }
}

namespace UB {
    NLOHMANN_JSON_SERIALIZE_ENUM(Shape::Type, {
        { Shape::Type::SUPERPRIM, "SUP" },
        { Shape::Type::ELLIPSOID, "ELL" },
        { Shape::Type::BEZIER, "BEZ" },
        { Shape::Type::SQUISHENGON, "SQU" },
        { Shape::Type::UBERPRIM, "UBR" }
    });

    NLOHMANN_JSON_SERIALIZE_ENUM(EditOp::Type, {
        { EditOp::Type::Sub, "SUB" },
        { EditOp::Type::Add, "ADD" },
        { EditOp::Type::Color, "COL" },
        { EditOp::Type::Crop, "CRP" },
    });


    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Shape, type, dim, thickness, cornerRadius, color)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditOp, type, blend, mirrorYZ)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditShape, op, transform, shape)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditSpline::ControlPoint, transform, shape)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditGroup, op, transform)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditRepeat, op, transform, shape, reps, offsetTransform)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditSpline, op, transform, loop, reps, controlpoints)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditFilter, op)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(EditGenerator, op)

    static void to_json(nlohmann::json &j, const Edit &v) {
        switch (v.index()) {
            case 0:
                break;
            case 1:
                j = std::get<EditShape>(v);
                j["type"] = "Shape";
                break;
            case 2:
                j = std::get<EditGroup>(v);
                j["type"] = "Group";
                break;
            case 3:
                j = std::get<EditRepeat>(v);
                j["type"] = "Repeat";
                break;
            case 4:
                j = std::get<EditSpline>(v);
                j["type"] = "Spline";
                break;
            case 5:
                j = std::get<EditFilter>(v);
                j["type"] = "Filter";
                break;
            case 6:
                j = std::get<EditGenerator>(v);
                j["type"] = "Generator";
                break;
        }
    }

    static void from_json(const nlohmann::json &j, Edit &v) {
        if (j == nullptr) {
            return;
        }
        if (j["type"] == "Shape") {
            v = j.get<EditShape>();
        } else if (j["type"] == "Group") {
            v = j.get<EditGroup>();
        } else if (j["type"] == "Repeat") {
            v = j.get<EditRepeat>();
        } else if (j["type"] == "Spline") {
            v = j.get<EditSpline>();
        } else if (j["type"] == "Filter") {
            v = j.get<EditFilter>();
        } else if (j["type"] == "Generator") {
            v = j.get<EditGenerator>();
        }
    }

    static void to_json(nlohmann::json &j, const EditList &v) { j = v.edits; }

    static void from_json(const nlohmann::json &j, EditList &v) { v.edits = j; }

    static void to_json(nlohmann::json &nlohmann_json_j, const Camera &camera) {
        nlohmann_json_j["nearPlane"] = camera.nearPlane;
        nlohmann_json_j["fov"] = camera.fov;
        nlohmann_json_j["focalDistance"] = camera.focalDistance;
        nlohmann_json_j["focalWidth"] = camera.focalWidth;
        nlohmann_json_j["exposure"] = camera.exposure;
        nlohmann_json_j["enableDOF"] = camera.enableDOF;
        nlohmann_json_j["bloomIntensity"] = camera.bloomIntensity;
        nlohmann_json_j["bloomRadius"] = camera.bloomRadius;
        nlohmann_json_j["bloomThreshold"] = camera.bloomThreshold;
        nlohmann_json_j["vignetteMidPoint"] = camera.vignetteMidPoint;
        nlohmann_json_j["vignetteRoundness"] = camera.vignetteRoundness;
        nlohmann_json_j["vignetteFeather"] = camera.vignetteFeather;
        nlohmann_json_j["vignetteColor"] = camera.vignetteColor;
        nlohmann_json_j["adjustContrast"] = camera.adjustContrast;
        nlohmann_json_j["adjustVibrance"] = camera.adjustVibrance;
        nlohmann_json_j["adjustSaturation"] = camera.adjustSaturation;
        nlohmann_json_j["adjustBlacks"] = camera.adjustBlacks;
    }

    static void from_json(const nlohmann::json &src, Camera &camera) {
        src.at("nearPlane").get_to(camera.nearPlane);
        src.at("fov").get_to(camera.fov);
        src.at("exposure").get_to(camera.exposure);
        src.at("enableDOF").get_to(camera.enableDOF);
        src.at("bloomIntensity").get_to(camera.bloomIntensity);
        src.at("bloomRadius").get_to(camera.bloomRadius);
        src.at("vignetteMidPoint").get_to(camera.vignetteMidPoint);
        src.at("vignetteRoundness").get_to(camera.vignetteRoundness);
        src.at("vignetteFeather").get_to(camera.vignetteFeather);
        src.at("vignetteColor").get_to(camera.vignetteColor);

        if (src.contains("bloomThreshold")) {
            src.at("bloomThreshold").get_to(camera.bloomThreshold);
        }

        if (src.contains("adjustContrast")) {
            src.at("adjustContrast").get_to(camera.adjustContrast);
        }

        if (src.contains("adjustVibrance")) {
            src.at("adjustVibrance").get_to(camera.adjustVibrance);
        }

        if (src.contains("adjustSaturation")) {
            src.at("adjustSaturation").get_to(camera.adjustSaturation);
        }

        if (src.contains("adjustBlacks")) {
            src.at("adjustBlacks").get_to(camera.adjustBlacks);
        }

        if (src.contains("focalDistance")) {
            src.at("focalDistance").get_to(camera.focalDistance);
        }

        if (src.contains("focalWidth")) {
            src.at("focalWidth").get_to(camera.focalWidth);
        }
    }

    //NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SplatParams, type, density, scale, faceCamera, jitterPosition, jitterOrientation, jitterScale, jitterBump, jitterHue, jitterBrightness, jitterMetalness, jitterRoughness)
    static void to_json(nlohmann::json &j, const SplatParams &s) {
        j["type"] = s.type;
        j["density"] = s.density;
        j["scale"] = s.scale;
        j["faceCamera"] = s.faceCamera;
        j["jitterPosition"] = s.jitterPosition;
        j["jitterOrientation"] = s.jitterOrientation;
        j["jitterScale"] = s.jitterScale;
        j["jitterBump"] = s.jitterBump;
        j["jitterHue"] = s.jitterHue;
        j["jitterBrightness"] = s.jitterBrightness;
        j["jitterMetalness"] = s.jitterMetalness;
        j["jitterRoughness"] = s.jitterRoughness;
        j["alpha"] = s.alpha;
    }

    static void from_json(const nlohmann::json &j, SplatParams &s) {
        j.at("type").get_to(s.type);
        j.at("density").get_to(s.density);
        j.at("scale").get_to(s.scale);
        j.at("faceCamera").get_to(s.faceCamera);
        j.at("jitterPosition").get_to(s.jitterPosition);
        j.at("jitterOrientation").get_to(s.jitterOrientation);
        j.at("jitterScale").get_to(s.jitterScale);
        j.at("jitterBump").get_to(s.jitterBump);
        j.at("jitterHue").get_to(s.jitterHue);
        j.at("jitterBrightness").get_to(s.jitterBrightness);
        j.at("jitterMetalness").get_to(s.jitterMetalness);
        j.at("jitterRoughness").get_to(s.jitterRoughness);
        if (j.contains("alpha")) {
            j.at("alpha").get_to(s.alpha);
        } else {
            s.alpha = 1.0f;
        }
    }

    static void to_json(nlohmann::json &j, const Material &m) {
        j["metalness"] = m.metalness;
        j["roughness"] = m.roughness;
        j["emission"] = m.emission;
        j["colorize"] = m.colorize;
        j["unlit"] = m.unlit;
        if (m.splat.has_value()) {
            j["splat"] = m.splat.value();
        }
        j["looseness"] = m.looseness;
    }

    static void from_json(const nlohmann::json &j, Material &m) {
        j.at("metalness").get_to(m.metalness);
        j.at("roughness").get_to(m.roughness);

        if (j.at("emission").is_array()) {
            glm::vec4 old;
            j.at("emission").get_to(old);
            m.emission = old.w;
            m.colorize = glm::vec4(old.x, old.y, old.z, 0.0f);
        } else {
            j.at("emission").get_to(m.emission);
            j.at("colorize").get_to(m.colorize);
        }

        m.unlit = j.contains("unlit") ? j.at("unlit").get<float>() : 0.0f;

        if (j.contains("splat")) {
            m.splat.emplace();
            j.at("splat").get_to(*m.splat);
        }
        if (j.contains("looseness")) {
            j.at("looseness").get_to(m.looseness);
        }
    }

    //NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Model, editList, material, maxRefineLevel, minBounds);

    static void to_json(nlohmann::json &j, const Model &m) {
        j["editList"] = m.editList;
        j["material"] = m.material;
        j["maxRefineLevel"] = m.maxRefineLevel;
        j["minBounds"] = m.minBounds;
        j["shapeDensity"] = m.shapeDensity;
    }

    static void from_json(const nlohmann::json &j, Model &m) {
        j.at("editList").get_to(m.editList);
        j.at("material").get_to(m.material);
        if (j.contains("maxRefineLevel")) j.at("maxRefineLevel").get_to(m.maxRefineLevel);
        m.minBounds = 0.0f;
        if (j.contains("minBounds")) j.at("minBounds").get_to(m.minBounds);
        if (j.contains("shapeDensity")) j.at("shapeDensity").get_to(m.shapeDensity);
    }

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Panel, dim)
    //NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Light, type, color, intensity, castShadow, dim, render, penumbra)

    static void to_json(nlohmann::json &nlohmann_json_j, const Light &light) {
        nlohmann_json_j["type"] = light.type;
        nlohmann_json_j["color"] = light.color;
        nlohmann_json_j["intensity"] = light.intensity;
        nlohmann_json_j["castShadow"] = light.castShadow;
        nlohmann_json_j["dim"] = light.dim;
        nlohmann_json_j["render"] = light.render;
        nlohmann_json_j["penumbra"] = light.penumbra;
        nlohmann_json_j["radius"] = light.radius;
    }

    static void from_json(const nlohmann::json &src, Light &light) {
        src.at("type").get_to(light.type);
        src.at("color").get_to(light.color);
        src.at("intensity").get_to(light.intensity);
        src.at("castShadow").get_to(light.castShadow);
        src.at("dim").get_to(light.dim);
        src.at("render").get_to(light.render);
        if (src.contains("penumbra")) {
            src.at("penumbra").get_to(light.penumbra);
        }
        if (src.contains("radius")) {
            src.at("radius").get_to(light.radius);
        }
    }

    inline void to_json(nlohmann::json &j, const Environment::Fog &fog) {
        j["color"] = fog.color;
        j["colorFromIbl"] = fog.colorFromIbl;
        j["density"] = fog.density;
        j["height"] = fog.height;
        j["heightFalloff"] = fog.heightFalloff;
        j["start"] = fog.start;
        j["maxOpacity"] = fog.maxOpacity;
        j["inscatteringStart"] = fog.inscatteringStart;
        j["inscatteringSize"] = fog.inscatteringSize;
    }

    inline void from_json(const nlohmann::json &j, Environment::Fog &fog) {
        j.at("color").get_to(fog.color);
        j.at("colorFromIbl").get_to(fog.colorFromIbl);
        if (j.contains("density")) {
            j.at("density").get_to(fog.density);
        }
        j.at("height").get_to(fog.height);
        j.at("heightFalloff").get_to(fog.heightFalloff);
        j.at("start").get_to(fog.start);
        j.at("maxOpacity").get_to(fog.maxOpacity);
        j.at("inscatteringStart").get_to(fog.inscatteringStart);
        j.at("inscatteringSize").get_to(fog.inscatteringSize);
    }

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Environment::GradientParams, A, B, C, Q)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Environment, type, gradient, fog, ssaoRadius, ssaoStrength, emissiveStrength,
                                       intensity)

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PhysicsMaterial,
                                       restitution,
                                       friction,
                                       linearDamping,
                                       angularDamping,
                                       gravityFactor)


    inline void to_json(nlohmann::json &j, const PhysicsBody &t) {
        j["type"] = t.type;
        j["material"] = t.material;
    }

    inline void from_json(const nlohmann::json &j, PhysicsBody &t) {
        j.at("type").get_to(t.type);
        j.at("material").get_to(t.material);
        if (j.contains("layer")) {
            j.at("layer").get_to(t.layer);
        } else {
            t.layer = 0;
        }
    }

    inline void to_json(nlohmann::json &j, const Constraint &c) {
        j["nodeA"] = c.nodeA;
        j["nodeB"] = c.nodeB;
        j["type"] = Constraint::getTypeName(c.type);
        j["enabled"] = c.settings.enabled;
        switch (c.type) {
            case Constraint::Type::FIXED:
                break;
            case Constraint::Type::DISTANCE:
                j["limitsMin"] = c.settings.limitsMin;
                j["limitsMax"] = c.settings.limitsMax;
                j["frequency"] = c.settings.frequency;
                j["damping"] = c.settings.damping;
                j["offset0"] = c.settings.offset[0];
                j["offset1"] = c.settings.offset[1];
                break;
            case Constraint::Type::HINGE:
                j["limitsMin"] = c.settings.limitsMin;
                j["limitsMax"] = c.settings.limitsMax;
                j["maxFriction"] = c.settings.maxFriction;
                break;
            case Constraint::Type::POINT:
                break;
            case Constraint::Type::CONE:
                j["halfAngle"] = c.settings.halfAngle;
                break;
            case Constraint::Type::SLIDER:
                j["limitsMin"] = c.settings.limitsMin;
                j["limitsMax"] = c.settings.limitsMax;
                j["frequency"] = c.settings.frequency;
                j["damping"] = c.settings.damping;
                j["maxFriction"] = c.settings.maxFriction;
                break;
            case Constraint::INVALID:
                break;
        }
    }

    inline void from_json(const nlohmann::json &j, Constraint &c) {
        j.at("nodeA").get_to(c.nodeA);
        j.at("nodeB").get_to(c.nodeB);
        j.at("enabled").get_to(c.settings.enabled);
        std::string typestring;
        j.at("type").get_to(typestring);
        c.type = Constraint::getType(typestring);

        switch (c.type) {
            case Constraint::Type::FIXED:
                break;
            case Constraint::Type::DISTANCE:
                j.at("limitsMin").get_to(c.settings.limitsMin);
                j.at("limitsMax").get_to(c.settings.limitsMax);
                j.at("frequency").get_to(c.settings.frequency);
                j.at("damping").get_to(c.settings.damping);
                if (j.contains("offset0")) {
                    j.at("offset0").get_to(c.settings.offset[0]);
                }
                if (j.contains("offset1")) {
                    j.at("offset1").get_to(c.settings.offset[1]);
                }
                break;
            case Constraint::Type::HINGE:
                j.at("limitsMin").get_to(c.settings.limitsMin);
                j.at("limitsMax").get_to(c.settings.limitsMax);
                j.at("maxFriction").get_to(c.settings.maxFriction);
                break;
            case Constraint::Type::POINT:
                break;
            case Constraint::Type::CONE:
                j.at("halfAngle").get_to(c.settings.halfAngle);
                break;
            case Constraint::Type::SLIDER:
                j.at("limitsMin").get_to(c.settings.limitsMin);
                j.at("limitsMax").get_to(c.settings.limitsMax);
                j.at("frequency").get_to(c.settings.frequency);
                j.at("damping").get_to(c.settings.damping);
                j.at("maxFriction").get_to(c.settings.maxFriction);
                break;
            default:
                break;
        }
    }

    static inline void to_json(nlohmann::json &j, const Logic &t) {
        j["dataRef"] = t.dataRef;
        j["inputs"] = t.inputs;
        j["type"] = int(t.type);
    }

    static inline void from_json(const nlohmann::json &j, Logic &t) {
        if (j.contains("type")) {
            j.at("dataRef").get_to(t.dataRef);
            j.at("inputs").get_to(t.inputs);
            j.at("type").get_to(t.type);
        } else {
            j.at("scriptRef").get_to(t.dataRef);
            j.at("scriptInputs").get_to(t.inputs);
            t.type = UB::Logic::TYPE_SCRIPT;
        }
    };

    static void to_json(nlohmann::json &j, const NodeData &v) {
        switch (v.index()) {
            case 0:
                break;
            case 1:
                j["modelRef"] = std::get<ModelRef>(v);
                j["class"] = "ModelRef";
                break;
            case 2:
                j = std::get<Panel>(v);
                j["class"] = "Panel";
                break;
            case 3:
                j = std::get<Camera>(v);
                j["class"] = "Camera";
                break;
            case 4:
                j = std::get<Light>(v);
                j["class"] = "Light";
                break;
            case 5:
                j = std::get<Environment>(v);
                j["class"] = "Environment";
                break;
            case 6:
                j = std::get<Constraint>(v);
                j["class"] = "Constraint";
                break;
            case 7:
                j = std::get<Logic>(v);
                j["class"] = "Logic";
                break;

        }
    }

    static void from_json(const nlohmann::json &j, NodeData &v) {
        if (j == nullptr) {
            return;
        }

        if (j["class"] == "ModelRef") {
            v = j["modelRef"].get<ModelRef>();
        } else if (j["class"] == "Panel") {
            v = j.get<Panel>();
        } else if (j["class"] == "Camera") {
            v = j.get<Camera>();
        } else if (j["class"] == "Light") {
            v = j.get<Light>();
        } else if (j["class"] == "Environment") {
            v = j.get<Environment>();
        } else if (j["class"] == "Constraint") {
            v = j.get<Constraint>();
        } else if (j["class"] == "Logic") {
            v = j.get<Logic>();
        }

    }

    // start at 2 since we neither serialize monostate nor Id
    template<std::size_t I = 2>
    static void to_json(nlohmann::json &j, const UB::Property &v) {
        if constexpr (I < std::variant_size_v<UB::Property>) {
            if (v.index() == I) {
                j["type"] = I;
                j["value"] = std::get<I>(v);
                return;
            } else {
                to_json<I + 1>(j, v);
            }
        }
    }

    template<std::size_t I = 2>
    static void from_json(const nlohmann::json &j, UB::Property &v) {
        if constexpr (I < std::variant_size_v<UB::Property>) {
            if (j.contains("type") && j.contains("value")) {
                const size_t type = j.at("type").get<size_t>();
                if (type == I) {
                    v = j.at("value").get<std::variant_alternative_t<I, UB::Property>>();
                    return;
                } else {
                    from_json<I + 1>(j, v);
                }
            }
        }
    }

    static void to_json(nlohmann::json &j, const UB::PropertyList &v) {
        for (auto it: v) {
            nlohmann::json entry = nlohmann::json::object();
            to_json(entry, it.second);
            j[it.first] = entry;
        }
    }


    static void from_json(const nlohmann::json &j, UB::PropertyList &v) {
        v.clear();
        if (j.is_object()) {
            for (auto &el: j.items()) {
                if (el.value().is_object()) {
                    UB::Property p;
                    from_json(el.value(), p);
                    v[el.key()] = p;
                } else {
                }
            }
        }
    }

    static void to_json(nlohmann::json &j, const Node &n) {
        j["id"] = uuids::to_string(n.id);
        j["name"] = n.name;
        j["tags"] = n.tags;
        j["transform"] = n.transform;
        j["children"] = n.children;
        j["data"] = n.data;
        j["hidden"] = n.hidden;
        j["overlay"] = n.overlay;
        j["shadowCaster"] = n.shadowCaster;
        j["shadowReceiver"] = n.shadowReceiver;
        if (n.physicsBody.type != PhysicsType::NONE) {
            j["body"] = n.physicsBody;
        }

        auto properties = nlohmann::json::object();
        to_json(properties, n.properties);
        j["properties"] = properties;
    }

    static void from_json(const nlohmann::json &j, Node &n) {
        j.at("id").get_to(n.id);
        j.at("name").get_to(n.name);
        j.at("tags").get_to(n.tags);
        j.at("transform").get_to(n.transform);
        j.at("children").get_to(n.children);
        j.at("data").get_to(n.data);
        n.hidden = j.contains("hidden") ? j.at("hidden").get<bool>() : false;
        n.overlay = j.contains("overlay") ? j.at("overlay").get<bool>() : false;
        n.shadowCaster = j.contains("shadowCaster") ? j.at("shadowCaster").get<bool>() : true;
        n.shadowReceiver = j.contains("shadowReceiver") ? j.at("shadowReceiver").get<bool>() : true;

        if (j.contains("body")) {
            n.physicsBody = j.at("body").get<PhysicsBody>();
        }

        if (j.contains("properties")) {
            from_json(j.at("properties"), n.properties);
        }
    }

#include "Legacy.h"

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(UB::Script, name, code);

    //template <typename T>
    //static inline void to_json(nlohmann::json& j, const UB::AnimFrame<T>& frame)
    //{
    //	j["time"] = frame.timeMs;
    //	j["interp"] = frame.interpolant;
    //	j["value"] = frame.value;
    //
    //}

    template<typename V, size_t I = 0>
    static inline void doframes(nlohmann::json &j, const AnimTrackBase *track) {
        if constexpr (I < std::variant_size_v<V>) {
            using T = std::variant_alternative_t<I, V>;
            if (track->is<T>()) {
                auto *ts = dynamic_cast<const AnimTrack<T> *>(track);
                j["type"] = I;
                j["times"] = nlohmann::json::array();
                j["interps"] = nlohmann::json::array();
                j["values"] = nlohmann::json::array();
                for (auto &f: ts->frames) {
                    j["times"].push_back(f.timeMs);
                    j["interps"].push_back(f.interpolant);
                    j["values"].push_back(f.value);
                }
            } else {
                doframes<V, I + 1>(j, track);
            }
        }
    }

    static inline void to_json(nlohmann::json &j, const std::map<uuids::uuid, AnimDataPtr> &anims) {
        for (const auto &animData: anims) {
            const std::string uuidstr = uuids::to_string(animData.first);
            j[uuidstr] = nlohmann::json::object();
            nlohmann::json &a = j[uuidstr];
            a["name"] = animData.second->name;
            a["tracks"] = nlohmann::json::object();
            nlohmann::json &tracks = a["tracks"];
            for (const auto &t: animData.second->tracks) {
                tracks[t.first] = nlohmann::json::object();
                nlohmann::json &frames = tracks[t.first];
                doframes<UB::AnimProperty>(frames, t.second.get());
            }
        }
    }

    template<size_t I = 0>
    AnimTrackBase *createAnimTrackFromTypeIndex(const nlohmann::json::object_t &trackj) {
        if constexpr (I < std::variant_size_v<AnimProperty>) {
            using T = std::variant_alternative_t<I, AnimProperty>;
            const size_t ti = trackj.at("type").get<size_t>();

            if (ti == I) {
                auto *t = new AnimTrack<T>();
                const auto &times = trackj.at("times").get<nlohmann::json::array_t>();
                const auto &interps = trackj.at("interps").get<nlohmann::json::array_t>();
                const auto &values = trackj.at("values").get<nlohmann::json::array_t>();
                const size_t numFrames = times.size();
                for (size_t ts = 0; ts < numFrames; ts++) {
                    t->frames.push_back({
                                                .interpolant = AnimInterpolant(interps[ts].get<uint32_t>()),
                                                .timeMs = times[ts].get<uint32_t>(),
                                                .value = values[ts].get<T>()
                                        });
                }
                return (AnimTrackBase *) (t);
            } else {
                return createAnimTrackFromTypeIndex<I + 1>(trackj);
            }
        }
        return nullptr;
    }

    static inline void from_json(const nlohmann::json &j, std::map<uuids::uuid, AnimDataPtr> &anims) {
        // loops over all anims in json collection
        if (j.is_object()) {
            for (const auto &it: j.get<nlohmann::json::object_t>()) {
                const std::string &uuidstr = it.first;
                const nlohmann::json &ad = it.second;
                AnimDataPtr adp = std::make_shared<AnimData>();
                adp->uuid = uuids::uuid::from_string(uuidstr).value();
                if (ad.contains("name")) {
                    ad.at("name").get_to(adp->name);
                }

                if (ad.contains("tracks")) {
                    const nlohmann::json::object_t &tracksj = ad.at("tracks");
                    for (const auto &ti: tracksj) {
                        adp->tracks[ti.first].reset(createAnimTrackFromTypeIndex(ti.second));
                    }
                }
                anims[adp->uuid] = adp;
            }
        }
    };

    static void to_json(nlohmann::json &j, const Collection &c) {
        j["schemaVersion"] = 1;
        j["engineVersion"] = "1";
        j["models"] = c.models;
        j["nodes"] = c.nodes;
        j["scripts"] = c.scripts;
        // j["animations"] = c.animations;
    }

    static void from_json(const nlohmann::json &j, Collection &c) {
        if (j.contains("schemaVersion") == false) {
        } else {
            j.at("models").get_to(c.models);
            j.at("nodes").get_to(c.nodes);
            if (j.contains("scripts")) {
                j.at("scripts").get_to(c.scripts);
            }
            // if (j.contains("animations")) {
                // j.at("animations").get_to(c.animations);
            // }
        }
    }
}
