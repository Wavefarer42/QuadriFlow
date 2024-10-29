#pragma once

#include <vector>
#include <list>
#include <map>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace qflow {
    using namespace Eigen;

    // Adjacency matrix
    struct Link {
        Link() {}

        Link(int _id, double _w = 1)
                : id(_id), weight(_w) {}

        inline bool operator<(const Link &link) const { return id < link.id; }

        int id;
        double weight;
    };

    struct TaggedLink {
        int id;
        unsigned char flag;

        TaggedLink() {}

        TaggedLink(int id) : id(id), flag(0) {}

        bool used() const { return flag & 1; }

        void markUsed() { flag |= 1; }

        TaggedLink &operator=(const Link &l) {
            flag = 0;
            id = l.id;
            return *this;
        }
    };

    typedef std::vector<std::vector<Link> > AdjacentMatrix;


    // Compare Key
    struct Key2i {
        Key2i(int x, int y)
                : key(std::make_pair(x, y)) {}

        bool operator==(const Key2i &other) const {
            return key == other.key;
        }

        bool operator<(const Key2i &other) const {
            return key < other.key;
        }

        std::pair<int, int> key;
    };

    struct Key3i {
        Key3i(int x, int y, int z)
                : key(std::make_pair(x, std::make_pair(y, z))) {}

        bool operator==(const Key3i &other) const {
            return key == other.key;
        }

        bool operator<(const Key3i &other) const {
            return key < other.key;
        }

        std::pair<int, std::pair<int, int> > key;
    };

    struct Key3f {
        Key3f(double x, double y, double z, double threshold)
                : key(std::make_pair(x / threshold, std::make_pair(y / threshold, z / threshold))) {}

        bool operator==(const Key3f &other) const {
            return key == other.key;
        }

        bool operator<(const Key3f &other) const {
            return key < other.key;
        }

        std::pair<int, std::pair<int, int> > key;
    };

    struct KeySorted2i {
        KeySorted2i(int x, int y)
                : key(std::make_pair(x, y)) {
            if (x > y)
                std::swap(key.first, key.second);
        }

        bool operator==(const KeySorted2i &other) const {
            return key == other.key;
        }

        bool operator<(const KeySorted2i &other) const {
            return key < other.key;
        }

        std::pair<int, int> key;
    };

    struct KeySorted3i {
        KeySorted3i(int x, int y, int z)
                : key(std::make_pair(x, std::make_pair(y, z))) {
            if (key.first > key.second.first)
                std::swap(key.first, key.second.first);
            if (key.first > key.second.second)
                std::swap(key.first, key.second.second);
            if (key.second.first > key.second.second)
                std::swap(key.second.first, key.second.second);
        }

        bool operator==(const Key3i &other) const {
            return key == other.key;
        }

        bool operator<(const Key3i &other) const {
            return key < other.key;
        }

        std::pair<int, std::pair<int, int> > key;
    };

    // Tree
    class DisjointTree {
    public:
        DisjointTree() {}

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
        DisjointOrientTree() {}

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
        std::vector<std::pair<int, int>> parent;
        std::vector<int> indices;
        std::vector<int> rank;
    };

    // Dedge
    struct DEdge {
        DEdge()
                : x(0), y(0) {}

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

}
