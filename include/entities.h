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

    typedef std::function<VectorXf(MatrixXf)> SDFn;

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

            return AlignedBox3f(
                Vector3f(bb.x.lb, bb.y.lb, bb.z.lb).array().ceil() - margin,
                Vector3f(bb.x.ub, bb.y.ub, bb.z.ub).array().ceil() + margin
            );
        }
    };
}
