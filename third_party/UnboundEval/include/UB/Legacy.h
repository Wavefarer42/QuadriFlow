//
//  Serialization.h
//  Framework
//
//  Created by Florian Hoenig on 5/29/22.
//

#pragma once

namespace Legacy
{
    struct Model
    {
        uint64_t editListId;
        float metalness;
        float roughness;
        glm::vec4 emission;
        std::optional<SplatParams> splat;
    };
    
    using NodeData = std::variant<std::monostate, Model, Panel, Camera, Light, Environment>;
    struct Node
    {
        uint64_t id;
        std::string name;
        Tags tags;
        Math::Transform transform;
        std::vector<Legacy::Node> children;
        Legacy::NodeData data;
    };
    
    static void to_json(nlohmann::json& j, const Legacy::Model& m)
    {
        j["editListId"] = m.editListId;
        j["metalness"] = m.metalness;
        j["roughness"] = m.roughness;
        j["emission"] = m.emission;
        if (m.splat.has_value())
        {
            j["splat"] = m.splat.value();
        }
    }
    
    static void from_json(const nlohmann::json& j, Legacy::Model& m)
    {
        j.at("editListId").get_to(m.editListId);
        j.at("metalness").get_to(m.metalness);
        j.at("roughness").get_to(m.roughness);
        j.at("emission").get_to(m.emission);
        if (j.contains("splat"))
        {
            m.splat.emplace();
            j.at("splat").get_to(*m.splat);
        }
    }
            
    static void to_json(nlohmann::json& j, const Legacy::NodeData& v)
    {
        switch (v.index())
        {
            case 0:
                break;
            case 1:
                j = std::get<Legacy::Model>(v);
                j["class"] = "Model";
                break;
            case 2:
                j = std::get<UB::Panel>(v);
                j["class"] = "Panel";
                break;
            case 3:
                j = std::get<UB::Camera>(v);
                j["class"] = "Camera";
                break;
            case 4:
                j = std::get<UB::Light>(v);
                j["class"] = "Light";
                break;
            case 5:
                j = std::get<UB::Environment>(v);
                j["class"] = "Environment";
                break;
        }
    }
    static void from_json(const nlohmann::json& j, Legacy::NodeData& v)
    {
        if (j == nullptr)
        {
            return;
        }
        
        if (j["class"] == "Model")
        {
            v = j.get<Legacy::Model>();
        }
        else if (j["class"] == "Panel")
        {
            v = j.get<UB::Panel>();
        }
        else if (j["class"] == "Camera")
        {
            v = j.get<UB::Camera>();
        }
        else if (j["class"] == "Light")
        {
            v = j.get<UB::Light>();
        }
        else if (j["class"] == "Environment")
        {
            v = j.get<UB::Environment>();
        }
    }
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Legacy::Node, id, name, tags, transform, children, data)
}

