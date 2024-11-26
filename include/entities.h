#pragma once

#include <utility>
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

    struct CLIArgs {
        std::string path_in;
        std::string path_out;
        int face_count = -1;
        int preserve_edges = 0;
        int preserve_boundaries = 0;
        int use_adaptive_meshing = 0;
        int seed = 0;
        int resolution = 100;
        bool is_valid = true;
    };

    // Adjacency matrix
    struct Link {
        Link() {
        }

        Link(int _id, double _w = 1)
            : id(_id), weight(_w) {
        }

        inline bool operator<(const Link &link) const { return id < link.id; }

        int id;
        double weight;
    };

    typedef std::vector<std::vector<Link> > AdjacentMatrix;

    // Tree
    class DisjointTree {
    public:
        DisjointTree() {
        }

        DisjointTree(int n) {
            parent.resize(n);
            rank.resize(n, 1);
            for (int i = 0; i < n; ++i) parent[i] = i;
        }

        int Parent(int x) {
            if (x == parent[x]) return x;
            int y = Parent(parent[x]);
            parent[x] = y;
            return y;
        }

        int Index(int x) { return indices[x]; }

        int IndexToParent(int x) { return indices_to_parent[x]; };

        void MergeFromTo(int x, int y) {
            int px = Parent(x);
            int py = Parent(y);
            if (px == py) return;
            rank[py] += rank[px];
            parent[px] = py;
        }

        void Merge(int x, int y) {
            int px = Parent(x);
            int py = Parent(y);
            if (px == py) return;
            if (rank[px] < rank[py]) {
                rank[py] += rank[px];
                parent[px] = py;
            } else {
                rank[px] += rank[py];
                parent[py] = px;
            }
        }

        // renumber the root so that it is consecutive.
        void BuildCompactParent() {
            std::vector<int> compact_parent;
            compact_parent.resize(parent.size());
            compact_num = 0;
            for (int i = 0; i < parent.size(); ++i) {
                if (parent[i] == i) {
                    compact_parent[i] = compact_num++;
                    indices_to_parent.push_back(i);
                }
            }
            indices.resize(parent.size());
            for (int i = 0; i < parent.size(); ++i) {
                indices[i] = compact_parent[Parent(i)];
            }
        }

        int CompactNum() { return compact_num; }

        int compact_num;
        std::vector<int> parent;
        std::vector<int> indices, indices_to_parent;
        std::vector<int> rank;
    };

    class DisjointOrientTree {
    public:
        DisjointOrientTree() {
        }

        DisjointOrientTree(int n) {
            parent.resize(n);
            rank.resize(n, 1);
            for (int i = 0; i < n; ++i) parent[i] = std::make_pair(i, 0);
        }

        int Parent(int j) {
            if (j == parent[j].first) return j;
            int k = Parent(parent[j].first);
            parent[j].second = (parent[j].second + parent[parent[j].first].second) % 4;
            parent[j].first = k;
            return k;
        }

        int Orient(int j) {
            if (j == parent[j].first) return parent[j].second;
            return (parent[j].second + Orient(parent[j].first)) % 4;
        }

        int Index(int x) { return indices[x]; }

        void MergeFromTo(int v0, int v1, int orient0, int orient1) {
            int p0 = Parent(v0);
            int p1 = Parent(v1);
            if (p0 == p1) return;
            int orientp0 = Orient(v0);
            int orientp1 = Orient(v1);

            if (p0 == p1) {
                return;
            }
            rank[p1] += rank[p0];
            parent[p0].first = p1;
            parent[p0].second = (orient0 - orient1 + orientp1 - orientp0 + 8) % 4;
        }

        void Merge(int v0, int v1, int orient0, int orient1) {
            int p0 = Parent(v0);
            int p1 = Parent(v1);
            if (p0 == p1) {
                return;
            }
            int orientp0 = Orient(v0);
            int orientp1 = Orient(v1);

            if (p0 == p1) {
                return;
            }
            if (rank[p1] < rank[p0]) {
                rank[p0] += rank[p1];
                parent[p1].first = p0;
                parent[p1].second = (orient1 - orient0 + orientp0 - orientp1 + 8) % 4;
            } else {
                rank[p1] += rank[p0];
                parent[p0].first = p1;
                parent[p0].second = (orient0 - orient1 + orientp1 - orientp0 + 8) % 4;
            }
        }

        void BuildCompactParent() {
            std::vector<int> compact_parent;
            compact_parent.resize(parent.size());
            compact_num = 0;
            for (int i = 0; i < parent.size(); ++i) {
                if (parent[i].first == i) {
                    compact_parent[i] = compact_num++;
                }
            }
            indices.resize(parent.size());
            for (int i = 0; i < parent.size(); ++i) {
                indices[i] = compact_parent[Parent(i)];
            }
        }

        int CompactNum() { return compact_num; }

        int compact_num;
        std::vector<std::pair<int, int> > parent;
        std::vector<int> indices;
        std::vector<int> rank;
    };

    // Dedge
    struct DEdge {
        DEdge()
            : x(0), y(0) {
        }

        DEdge(int _x, int _y) {
            if (_x > _y)
                x = _y, y = _x;
            else
                x = _x, y = _y;
        }

        bool operator<(const DEdge &e) const {
            return (x < e.x) || (x == e.x && y < e.y);
        }

        bool operator==(const DEdge &e) const {
            return x == e.x && y == e.y;
        }

        bool operator!=(const DEdge &e) const {
            return x != e.x || y != e.y;
        }

        int x, y;
    };

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

        AlignedBox3f bounding_box(const int &idx) {
            if (idx >= _models.size()) {
                throw std::runtime_error("Index out of bounds");
            }

            const auto bb = UB::getTightAABB(_models[idx]);

            return AlignedBox3f(
                Vector3f(bb.x.lb, bb.y.lb, bb.z.lb).array().ceil(),
                Vector3f(bb.x.ub, bb.y.ub, bb.z.ub).array().ceil()
            );
        }
    };
}
