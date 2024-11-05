
#include "Collection.h"
#include <regex>

namespace
{
    bool iequals(const std::string& a, const std::string& b)
    {
        return std::equal(a.begin(), a.end(),
                          b.begin(), b.end(),
                          [](char a, char b) {
            return tolower(a) == tolower(b);
        });
    }
}

namespace UB
{
    const std::string Constraint::NAME_INVALID = "Invalid";
    const std::string Constraint::NAME_FIXED = "Fixed";
    const std::string Constraint::NAME_DISTANCE = "Distance";
    const std::string Constraint::NAME_HINGE = "Hinge";
    const std::string Constraint::NAME_POINT = "Point";
    const std::string Constraint::NAME_CONE = "Cone";
    const std::string Constraint::NAME_SLIDER = "Slider";

    NodeFilter::NodeFilter(const std::string& input) : value(), condition(INVALID)
    {
        std::smatch sm;
        const std::regex r(R"((\w+)\:(.+))");
        if (std::regex_match(input, sm, r) && sm.size() == 3)
        {
            if (iequals(sm[1], "type")) condition = IS_TYPE;
            else if (iequals(sm[1], "name")) condition = IS_NAME;
            else if (iequals(sm[1], "uuid")) condition = IS_UUID;
            else if (iequals(sm[1], "property")) condition = HAS_PROPERTY;
            else if (iequals(sm[1], "tag")) condition = HAS_TAG;
            value = sm[2];
        }
    }
    
    bool NodeFilter::check(const Node& node) const
    {
        switch(condition)
        {
            case INVALID:
                return true;
            case IS_TYPE:
            {
                if (std::holds_alternative<UB::ModelRef>(node.data)         && iequals(value, "model"))          return true;
                else if (std::holds_alternative<UB::Panel>(node.data)       && iequals(value, "panel"))          return true;
                else if (std::holds_alternative<UB::Camera>(node.data)      && iequals(value, "camera"))         return true;
                else if (std::holds_alternative<UB::Environment>(node.data) && iequals(value, "environment"))    return true;
                else if (std::holds_alternative<UB::Light>(node.data)       && iequals(value, "light"))          return true;
                // else if (std::holds_alternative<UB::PhysicsBody>(node.data) && iequals(value, "body"))           return true;
                else if (std::holds_alternative<std::monostate>(node.data)  && iequals(value, "entity"))         return true;
                else return false;
            }
            case IS_NAME:
            {
                //Log::debug("%s == %s %d", node.name.c_str(), value.c_str(), node.name == value);
                return node.name == value;
            }
            case IS_UUID:
            {
                return uuids::to_string(node.id) == value;
            }
            case HAS_PROPERTY:
            {
                return node.properties.find(value) != node.properties.end();
            }
            case HAS_TAG:
            {
                for (auto& tag : node.tags)
                {
                    if (tag == value) return true;
                }
                return false;
            }
        }
    }
    
    bool NodeFilter::check(const Node& node, const List& filterList)
    {
        if (filterList.size() == 0) return true;
        
        bool pass = false;
        for (const auto& f : filterList)
        {
            pass = pass || f.check(node);
        }
        
        return pass;
    }


    void visitNodeRecursive(UB::Node* node, UB::Node* parent, std::function<void(UB::Node*, UB::Node*)> cb)
    {
        cb(node, parent);
        for (auto& child : node->children)
        {
            visitNodeRecursive(&child, node, cb);
        }
    }

    void visitRecursive(UB::CollectionPtr collection, std::function<void(UB::Node*, UB::Node*)> cb)
    {
        for (UB::Node& n : collection->nodes)
        {
            visitNodeRecursive(&n, nullptr, cb);
        }
    }

    const std::string& Constraint::getTypeName(Type type)
    {
        if (type == UB::Constraint::FIXED)
        {
            return  UB::Constraint::NAME_FIXED;
        }
        else if (type == UB::Constraint::DISTANCE)
        {
            return  UB::Constraint::NAME_DISTANCE;
        }
        else if (type == UB::Constraint::HINGE)
        {
            return  UB::Constraint::NAME_HINGE;
        }
        else if (type == UB::Constraint::POINT)
        {
            return  UB::Constraint::NAME_POINT;
        }
        else if (type == UB::Constraint::CONE)
        {
            return  UB::Constraint::NAME_CONE;
        }
        else if (type == UB::Constraint::SLIDER)
        {
            return  UB::Constraint::NAME_SLIDER;
        }
        else
        {
            return  UB::Constraint::NAME_INVALID;
        }
    }

    const Constraint::Type Constraint::getType(const std::string& typeName)
    {
        if (typeName == UB::Constraint::NAME_FIXED)
        {
            return  UB::Constraint::FIXED;
        }
        else if (typeName == UB::Constraint::NAME_DISTANCE)
        {
            return  UB::Constraint::DISTANCE;
        }
        else if (typeName == UB::Constraint::NAME_HINGE)
        {
            return  UB::Constraint::HINGE;
        }
        else if (typeName == UB::Constraint::NAME_POINT)
        {
            return  UB::Constraint::POINT;
        }
        else if (typeName == UB::Constraint::NAME_CONE)
        {
            return  UB::Constraint::CONE;
        }
        else if (typeName == UB::Constraint::NAME_SLIDER)
        {
            return  UB::Constraint::SLIDER;
        }
        else
        {
            return  UB::Constraint::INVALID;
        }
    }

}
