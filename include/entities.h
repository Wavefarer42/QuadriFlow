#pragma once

#include <vector>
#include <list>
#include <map>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "Collection.h"
#include "Serialization.h"
#include "uuid.h"

namespace entities {
    using namespace Eigen;

    typedef OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::DefaultTraits> Mesh;

    struct PointHash {
        size_t operator()(const Mesh::Point &p) const {
            // A simple hash function for 3D points
            return std::hash<float>()(p[0]) ^ std::hash<float>()(p[1]) ^ std::hash<float>()(p[2]);
        }
    };

    struct PointEqual {
        bool operator()(const Mesh::Point &p1, const Mesh::Point &p2) const {
            // Component-wise equality comparison
            return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]);
        }
    };

    typedef std::function<VectorXf(const MatrixXf &)> SDFn;

    class UnboundModel {
        const UB::Collection _collection;
        std::vector<UB::InstructionList> _models;
        std::map<std::string, SDFn> _sdfns;
        std::vector<std::string> _keys;

    public:
        explicit UnboundModel(UB::Collection collection)
            : _collection(collection) {
            _models.resize(collection.models.size());
            _keys.resize(collection.models.size());

            int idx = 0;
            for (auto const &[id_model, model]: collection.models) {
                UB::InstructionList evalList;
                UB::compileEditList(model.editList, evalList);

                const auto sampler = [evalList](MatrixXf domain) {
                    VectorXf distances(domain.rows());

                    tbb::parallel_for(
                        tbb::blocked_range<int>(0, domain.rows()),
                        [&](const tbb::blocked_range<int> &range) {
                            for (int i = range.begin(); i < range.end(); ++i) {
                                const auto point = glm::vec3(domain(i, 0), domain(i, 1), domain(i, 2));
                                distances[i] = UB::evalDistance(point, evalList);
                            }
                        }
                    );
                    return distances;
                };

                _models[idx] = evalList;
                _sdfns[uuids::to_string(id_model)] = sampler;
                _keys[idx] = uuids::to_string(id_model);
                idx++;
            }
        }

        size_t size() const {
            return _sdfns.size();
        }

        [[nodiscard]] std::vector<std::string> keys() const {
            std::vector<std::string> keys_new;
            for (auto const key: _keys) {
                keys_new.push_back(key);
            }
            return keys_new;
        }

        SDFn &operator[](const std::string &key) {
            return _sdfns[key];
        }

        SDFn &operator[](const uuids::uuid &key) {
            return _sdfns[uuids::to_string(key)];
        }

        operator std::vector<SDFn>() const {
            std::vector<SDFn> sfdns;
            for (auto const &[key, val]: _sdfns) {
                sfdns.push_back(val);
            }
            return sfdns;
        }

        SDFn &operator[](const int &key) {
            assert(key < _keys.size());
            return _sdfns[_keys[key]];
        }

        AlignedBox3f bounding_box(const int &idx, const float margin = 1) {
            if (idx >= _models.size()) {
                throw std::runtime_error("Index out of bounds");
            }

            const auto bb = UB::getTightAABB(_models[idx]);
            const auto extends = glm::compMax(glm::vec3(
                glm::max(glm::abs(bb.x.lb), glm::abs(bb.x.ub)),
                glm::max(glm::abs(bb.y.lb), glm::abs(bb.y.ub)),
                glm::max(glm::abs(bb.z.lb), glm::abs(bb.z.ub))
            ));
            return AlignedBox3f(
                Vector3f(-extends, -extends, -extends),
                Vector3f(extends, extends, extends)
            );
        }
    };
}
