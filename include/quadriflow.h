#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

#include "entities.h"

#define RCPOVERFLOW 2.93873587705571876e-39f

namespace quadriflow {
    using namespace boost;
    using namespace Eigen;

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

    class Hierarchy {
    public:
        Hierarchy();

        void Initialize(double scale, int with_scale = 0);

        void DownsampleGraph(
            const AdjacentMatrix adj,
            const MatrixXd &V,
            const MatrixXd &N,
            const VectorXd &A,
            MatrixXd &V_p,
            MatrixXd &N_p,
            VectorXd &A_p,
            MatrixXi &to_upper,
            VectorXi &to_lower,
            AdjacentMatrix &adj_p
        );

        void generate_graph_coloring_deterministic(
            const AdjacentMatrix &adj,
            int size,
            std::vector<std::vector<int> > &phases
        );

        void FixFlip();

        void PropagateEdge();

        void DownsampleEdgeGraph(
            std::vector<Vector3i> &FQ,
            std::vector<Vector3i> &F2E,
            std::vector<Vector2i> &edge_diff,
            std::vector<int> &allow_changes,
            int level
        );

        void UpdateGraphValue(
            std::vector<Vector3i> &FQ,
            std::vector<Vector3i> &F2E,
            std::vector<Vector2i> &edge_diff
        );

        enum {
            MAX_DEPTH = 25
        };

        void clearConstraints();

        void propagateConstraints();

        double m_scale = 1;
        int rng_seed;

        MatrixXi m_faces; // m_faces(i, j) i \in [0, 3) ith index in face j
        VectorXi m_E2E; // inverse edge
        std::vector<AdjacentMatrix> m_adjacency;
        std::vector<MatrixXd> m_vertices;
        std::vector<MatrixXd> m_normals;
        std::vector<VectorXd> m_vertex_area;
        std::vector<std::vector<std::vector<int> > > m_phases;

        // parameters
        std::vector<MatrixXd> m_orientation;
        std::vector<MatrixXd> m_positions;
        std::vector<VectorXi> mToLower;
        std::vector<MatrixXi> mToUpper; // mToUpper[h](i, j) \in m_vertices; i \in [0, 2); j \in m_vertices
        std::vector<MatrixXd> m_scales;
        std::vector<MatrixXd> m_areas;

        // constraints
        std::vector<MatrixXd> m_orientation_constraint;
        std::vector<MatrixXd> m_position_constraints;
        std::vector<VectorXd> m_orientation_constraint_weight;
        std::vector<VectorXd> m_position_constraint_weights;

        int with_scale;

        // upper: fine to coarse
        std::vector<std::vector<int> > mToUpperFaces; // face correspondance
        std::vector<std::vector<int> > mSing;
        std::vector<std::vector<int> > mToUpperEdges; // edge correspondance
        std::vector<std::vector<int> > mToUpperOrients; // rotation of edges from fine to coarse
        std::vector<std::vector<Vector3i> > mFQ; // m_face_edge_orientation
        std::vector<std::vector<Vector3i> > mF2E; // m_face_edge_ids
        std::vector<std::vector<Vector2i> > mE2F; // undirect edges to face ID
        std::vector<std::vector<int> > mAllowChanges;
        std::vector<std::vector<Vector2i> > mEdgeDiff; // face_edgeDiff
    };

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

    class Parametrizer {
    public:
        entities::Mesh m_mesh;
        /**
         * Maps the face to it's singularity valence, i.e., 1 for valence=3 or 3 for valence=5
         */
        std::map<int, int> m_singularities;

        /**
         * Maps the singularity face to the vertex index.
         */
        std::map<int, Vector2i> m_singularity_position;

        /**
         * Stores the re-oriented singularity vertex such that all face in the same direction.
         *
         * Given m_singularity_rank(i, j) with the i'th vertex (i in [0, 3)) in the j'th face,
         * rotate the face until all three vertices are orientated in the same direction.
         */
        MatrixXi m_singularity_rank;

        /**
         * Stores the singularity index for each face.
         * m_singularity_index(i x 2 + dim, j) with the i'th vertex (i in [0, 6)) in the j'th face and the (t_ij-t_ji)'s dimenstion in the paper.
         * TODO : Look up in paper
         */
        MatrixXi m_singularity_index;

        // input mesh
        MatrixXd m_vertices;
        MatrixXd m_normals_vertices;
        MatrixXd m_faces_normals;
        MatrixXi m_faces;
        MatrixXd m_faces_slope;
        MatrixXd m_faces_orientation;

        double m_normalize_scale;
        Vector3d m_normalize_offset;

        // data structures
        VectorXd m_rho;
        VectorXi m_V2E;
        VectorXi m_E2E;
        VectorXi m_boundary;
        VectorXi m_non_manifold; // m_non_manifold vertices, in boolean
        AdjacentMatrix m_adjacency_matrix;
        Hierarchy m_hierarchy;

        // Mesh Status;
        double m_surface_area;
        double m_scale;
        double m_average_edge_length;
        double m_max_edge_length;
        VectorXd m_vertex_area;

        DisjointTree m_disjoint_tree;

        int compact_num_v;
        std::vector<std::vector<int> > m_vertices_set;
        std::vector<Vector3d> m_positions_compact;
        std::vector<Vector3d> m_orientations_compact;
        std::vector<Vector3d> m_normals_compact;
        std::vector<Vector4i> m_faces_compact;
        std::set<std::pair<int, int> > m_quad_edges;
        std::vector<int> m_V2E_compact;
        std::vector<int> m_E2E_compact;
        VectorXi m_boundary_compact;
        VectorXi m_non_manifold_compact;

        std::vector<int> m_bad_vertices;
        std::vector<double> m_counter;

        /**
         * m_edges_preserve[deid]: Whether the edge is a edge or feature that should be preserved.
         * Edge indices in a linear flatten array. Assumes that faces are triangles.
         * Size is 3 * m_faces.cols()
         */
        std::vector<int> m_edges_preserve;

        /**
         * m_allow_changes[variable_id]: Whether variable can be changed based on sharp edges
         */
        std::vector<int> m_allow_changes;

        /**
         * m_edge_difference[edgeIds[i](j)]:  t_ij+t_ji under m_edge_values[edgeIds[i](j)].x's Q value
         */
        std::vector<Vector2i> m_edge_difference;
        std::vector<DEdge> m_edge_values; // see above

        /**
         *  m_face_edge_ids[i](j): "undirected edge ID" of the i'th face and the j'th edge
         */
        std::vector<Vector3i> m_face_edge_ids;

        /**
         * m_face_edge_orientation[i](j): Rotate from m_edge_difference space
         * a) initially, to F(0, i)'s Q space
         * b) later on, to a global Q space where some edges are fixed
         */
        std::vector<Vector3i> m_face_edge_orientation;

        /**
         * variable[i].first: indices of the two equations corresponding to variable i
         * variable[i].second: number of positive minus negative of variables' occurrences
         */
        std::vector<std::pair<Vector2i, int> > m_variables;

        struct QuadInfo {
            QuadInfo() : patchId(-1), coordinate(0x10000000, 0x10000000), singular(0), edge(0) {
            }

            int patchId;
            Vector2i coordinate;
            int singular;
            int edge;
        };

        std::vector<QuadInfo> m_quad_info;


        std::vector<MatrixXd> m_triangle_space;


        // Mesh Initialization


        /**
         * Analyzes mesh basic topological statistics such as, surface area, average edge length, max edge length.
         */
        void analyze_mesh();

        /**
         * Compute the vertex area between a vertex and its adjacent vertex bounded by the midpoints.
         *
         * TODO : Belongs to analyze mesh. Used for vertex scale and edge length.
         */
        void compute_vertex_area();

        /**
         * Compute the face and vertex normals.
         * TODO : Use the normals from the field.
         */
        void compute_normals();

        /**
         * Finds edges in the mesh and sets them for preservation if active.
         */
        void find_and_set_edges();

        /**
         * Finds the edges in the mesh based on the underlying model and sets them for preservation if active.
         */
        void find_and_set_edges(
            entities::SDFn sdfn,
            float max_angle
        );

        /**
         * Initializes the mesh data structures.
         * TODO: Service level method to initialize datastructures
         */
        void initialize(
            entities::SDFn sdfn,
            int target_face_count,
            bool with_scale
        );

        // Integer Grid Map

        /**
         * Maps edges and create a list of singularity rotations for each face.
         * TODO : Use half edge for easier navigation.
         */
        void build_edge_info();

        void build_integer_constraints();

        /**
         * Computes the constraints for the edges and features and optimizes the integer map. Afterwards it fixes several issues arising from the optimization.
         *
         * TODO : Service level method for the surfacenets flow. Split and simplify
         */
        void compute_index_map(Hierarchy &hierarchy, int with_scale = 0);

        // Post-Processing

        /**
         * Fix the valence of vertices by removing valence 0, 2 connections and reducing higher valence connections.
         * TODO : Cleanup implementation with half edge.
         */
        void fix_valence();

        void fix_flip_hierarchy();

        void find_fix_holes();

        void close_hole(std::vector<int> &loop_vertices);

        /**
         * Computes the quad energy of the vertex ring and returns the resulting quads in the list.
         * TODO: Change return type
         */
        double compute_quad_energy(
            std::vector<int> &loop_vertices,
            std::vector<Vector4i> &res_quads,
            int level
        );

        // Mesh Extraction
        /**
         * Extracts the quad mesh from the fields.
         */
        void extract_quad();

        /**
         * Builds the triangle
         */
        void build_triangle_manifold(
            DisjointTree &disajoint_tree,
            std::vector<int> &edge,
            std::vector<int> &face,
            std::vector<DEdge> &edge_values,
            std::vector<Vector3i> &F2E,
            std::vector<Vector2i> &E2F,
            std::vector<Vector2i> &EdgeDiff,
            std::vector<Vector3i> &FQ
        );
    };

    class Optimizer {
    public:
        Optimizer();

        static void optimize_orientations(Hierarchy &hierarchy);

        static void optimize_scale(Hierarchy &mRes, VectorXd &rho, int adaptive);

        static void optimize_positions(Hierarchy &mRes, int iterations = 6, bool with_scale = true);

        static void optimize_integer_constraints(Hierarchy &mRes, std::map<int, int> &singularities);

        static void optimize_positions_fixed(
            Hierarchy &mRes,
            std::vector<DEdge> &edge_values,
            std::vector<Vector2i> &edge_diff,
            std::set<int> &sharp_vertices,
            std::map<int, std::pair<Vector3d, Vector3d> > &sharp_constraints,
            bool with_scale = false
        );

        static void optimize_positions_sharp(
            Hierarchy &mRes,
            std::vector<DEdge> &edge_values,
            std::vector<Vector2i> &edge_diff,
            std::vector<int> &sharp_edges,
            std::set<int> &sharp_vertices,
            std::map<int, std::pair<Vector3d, Vector3d> > &sharp_constraints,
            bool with_scale = false
        );

        static void optimize_positions_dynamic(
            MatrixXi &F,
            MatrixXd &V,
            MatrixXd &N,
            MatrixXd &Q,
            std::vector<std::vector<int> > &Vset,
            std::vector<Vector3d> &O_compact,
            std::vector<Vector4i> &F_compact,
            std::vector<int> &V2E_compact,
            std::vector<int> &E2E_compact,
            double mScale,
            std::vector<Vector3d> &diffs, std::vector<int> &diff_count,
            std::map<std::pair<int, int>, int> &o2e,
            std::vector<int> &sharp_o,
            std::map<int, std::pair<Vector3d, Vector3d> > &compact_sharp_constraints,
            bool with_scale
        );
    };

    class MaxFlowHelper {
    public:
        MaxFlowHelper() {
        }

        virtual ~MaxFlowHelper() {
        };

        virtual void resize(int n, int m) = 0;

        virtual void addEdge(int x, int y, int c, int rc, int v, int cost = 1) = 0;

        virtual int compute() = 0;

        virtual void applyTo(std::vector<Vector2i> &edge_diff) = 0;
    };

    class BoykovMaxFlowHelper : public MaxFlowHelper {
    public:
        typedef int EdgeWeightType;
        typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
        typedef adjacency_list<vecS, vecS, directedS,
            property<vertex_name_t, std::string,
                property<vertex_index_t, long,
                    property<vertex_color_t, boost::default_color_type,
                        property<vertex_distance_t, long,
                            property<vertex_predecessor_t, Traits::edge_descriptor> > > > >,

            property<edge_capacity_t, EdgeWeightType,
                property<edge_residual_capacity_t, EdgeWeightType,
                    property<edge_reverse_t, Traits::edge_descriptor> > > >
        Graph;

    public:
        BoykovMaxFlowHelper() { rev = get(edge_reverse, g); }

        void resize(int n, int m) {
            vertex_descriptors.resize(n);
            for (int i = 0; i < n; ++i) vertex_descriptors[i] = add_vertex(g);
        }

        int compute() {
            EdgeWeightType flow =
                    boykov_kolmogorov_max_flow(g, vertex_descriptors.front(), vertex_descriptors.back());
            return flow;
        }

        void addDirectEdge(
            Traits::vertex_descriptor &v1,
            Traits::vertex_descriptor &v2,
            property_map<Graph, edge_reverse_t>::type &rev,
            const int capacity,
            const int inv_capacity,
            Graph &g,
            Traits::edge_descriptor &e1,
            Traits::edge_descriptor &e2
        ) {
            e1 = add_edge(v1, v2, g).first;
            e2 = add_edge(v2, v1, g).first;
            put(edge_capacity, g, e1, capacity);
            put(edge_capacity, g, e2, inv_capacity);

            rev[e1] = e2;
            rev[e2] = e1;
        }

        void addEdge(int x, int y, int c, int rc, int v, int cost = 1) {
            Traits::edge_descriptor e1, e2;
            addDirectEdge(vertex_descriptors[x], vertex_descriptors[y], rev, c, rc, g, e1, e2);
            if (v != -1) {
                edge_to_variables[e1] = std::make_pair(v, -1);
                edge_to_variables[e2] = std::make_pair(v, 1);
            }
        }

        void applyTo(std::vector<Vector2i> &edge_diff) {
            property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
            property_map<Graph, edge_residual_capacity_t>::type residual_capacity =
                    get(edge_residual_capacity, g);

            graph_traits<Graph>::vertex_iterator u_iter, u_end;
            graph_traits<Graph>::out_edge_iterator ei, e_end;
            for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
                for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
                    if (capacity[*ei] > 0) {
                        int flow = (capacity[*ei] - residual_capacity[*ei]);
                        if (flow > 0) {
                            auto it = edge_to_variables.find(*ei);
                            if (it != edge_to_variables.end()) {
                                edge_diff[it->second.first / 2][it->second.first % 2] +=
                                        it->second.second * flow;
                            }
                        }
                    }
        }

    private:
        Graph g;
        property_map<Graph, edge_reverse_t>::type rev;
        std::vector<Traits::vertex_descriptor> vertex_descriptors;
        std::map<Traits::edge_descriptor, std::pair<int, int> > edge_to_variables;
    };

    class ECMaxFlowHelper : public MaxFlowHelper {
    public:
        struct FlowInfo {
            int id;
            int capacity, flow;
            int v, d;
            FlowInfo *rev;
        };

        struct SearchInfo {
            SearchInfo(int _id, int _prev_id, FlowInfo *_info)
                : id(_id), prev_id(_prev_id), info(_info) {
            }

            int id;
            int prev_id;
            FlowInfo *info;
        };

        ECMaxFlowHelper() { num = 0; }

        int num;
        std::vector<FlowInfo *> variable_to_edge;

        void resize(int n, int m) {
            graph.resize(n);
            variable_to_edge.resize(m, 0);
            num = n;
        }

        void addEdge(int x, int y, int c, int rc, int v, int cost = 0) {
            FlowInfo flow;
            flow.id = y;
            flow.capacity = c;
            flow.flow = 0;
            flow.v = v;
            flow.d = -1;
            graph[x].push_back(flow);
            auto &f1 = graph[x].back();
            flow.id = x;
            flow.capacity = rc;
            flow.flow = 0;
            flow.v = v;
            flow.d = 1;
            graph[y].push_back(flow);
            auto &f2 = graph[y].back();
            f2.rev = &f1;
            f1.rev = &f2;
        }

        int compute() {
            int total_flow = 0;
            int count = 0;
            while (true) {
                count += 1;
                std::vector<int> vhash(num, 0);
                std::vector<SearchInfo> q;
                q.push_back(SearchInfo(0, -1, 0));
                vhash[0] = 1;
                int q_front = 0;
                bool found = false;
                while (q_front < q.size()) {
                    int vert = q[q_front].id;
                    for (auto &l: graph[vert]) {
                        if (vhash[l.id] || l.capacity <= l.flow) continue;
                        q.push_back(SearchInfo(l.id, q_front, &l));
                        vhash[l.id] = 1;
                        if (l.id == num - 1) {
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                    q_front += 1;
                }
                if (q_front == q.size()) break;
                int loc = q.size() - 1;
                while (q[loc].prev_id != -1) {
                    q[loc].info->flow += 1;
                    q[loc].info->rev->flow -= 1;
                    loc = q[loc].prev_id;
                    // int prev_v = q[loc].id;
                    // applyFlow(prev_v, current_v, 1);
                    // applyFlow(current_v, prev_v, -1);
                }
                total_flow += 1;
            }
            return total_flow;
        }

        void applyTo(std::vector<Vector2i> &edge_diff) {
            for (int i = 0; i < graph.size(); ++i) {
                for (auto &flow: graph[i]) {
                    if (flow.flow > 0 && flow.v != -1) {
                        if (flow.flow > 0) {
                            edge_diff[flow.v / 2][flow.v % 2] += flow.d * flow.flow;
                            if (abs(edge_diff[flow.v / 2][flow.v % 2]) > 2) {
                            }
                        }
                    }
                }
            }
        }

        void applyFlow(int v1, int v2, int flow) {
            for (auto &it: graph[v1]) {
                if (it.id == v2) {
                    it.flow += flow;
                    break;
                }
            }
        }

        std::vector<std::list<FlowInfo> > graph;
    };
}
