#pragma once

#include <atomic>
#include <list>
#include <map>
#include <set>
#include <condition_variable>
#include <unordered_set>

#include "entities.h"
#include "field-math.h"
#include "hierarchy.h"
#include "persistence.h"

namespace services {
    using namespace Eigen;

    class MeshService {
    public:
        explicit MeshService(persistence::MeshDao mesh_dao) : mesh_dao(mesh_dao) {}

        entities::QuadMesh load_trimesh_from_file(const std::string &filename) const;

        void save_quadmesh_to_file(const std::string &filename, const entities::QuadMesh &mesh);

        void set_boundary_constraints(Hierarchy &hierarchy) const;

        std::map<int, int> find_orientation_singularities(Hierarchy &hierarchy) const;

        std::tuple<std::map<int, Vector2i>, MatrixXi, MatrixXi> find_position_singularities(
                Hierarchy &m_hierarchy,
                bool with_scale
        ) const;

        std::tuple<MatrixXd, MatrixXd> estimate_slope(
                Hierarchy &hierarchy,
                std::vector<MatrixXd> &triangle_space,
                MatrixXd &normals_faces
        ) const;

    private:
        persistence::MeshDao mesh_dao;

    };

    class Parametrizer {
    public:

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
        VectorXi m_non_manifold;  // m_non_manifold vertices, in boolean
        entities::AdjacentMatrix m_adjacency_matrix;
        Hierarchy m_hierarchy;

        // Mesh Status;
        double m_surface_area;
        double m_scale;
        double m_average_edge_length;
        double m_max_edge_length;
        VectorXd m_vertex_area;

        entities::DisjointTree m_disjoint_tree;

        int compact_num_v;
        std::vector<std::vector<int>> m_vertices_set;
        std::vector<Vector3d> m_positions_compact;
        std::vector<Vector3d> m_orientations_compact;
        std::vector<Vector3d> m_normals_compact;
        std::vector<Vector4i> m_faces_compact;
        std::set<std::pair<int, int>> m_quad_edges;
        std::vector<int> m_V2E_compact;
        std::vector<int> m_E2E_compact;
        VectorXi m_boundary_compact;
        VectorXi m_non_manifold_compact;

        std::vector<int> m_bad_vertices;
        std::vector<double> m_counter;

        /**
         * m_sharp_edges[deid]: Whether the edge is a edge or feature that should be preserved
         */
        std::vector<int> m_sharp_edges;

        /**
         * m_allow_changes[variable_id]: Whether variable can be changed based on sharp edges
         */
        std::vector<int> m_allow_changes;

        /**
         * m_edge_difference[edgeIds[i](j)]:  t_ij+t_ji under m_edge_values[edgeIds[i](j)].x's Q value
         */
        std::vector<Vector2i> m_edge_difference;
        std::vector<entities::DEdge> m_edge_values;   // see above

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
        std::vector<std::pair<Vector2i, int>> m_variables;

        struct QuadInfo {
            QuadInfo() : patchId(-1), coordinate(0x10000000, 0x10000000), singular(0), edge(0) {}

            int patchId;
            Vector2i coordinate;
            int singular;
            int edge;
        };

        std::vector<QuadInfo> m_quad_info;


        std::vector<MatrixXd> m_triangle_space;


        // Mesh IO

        void save_to_obj(const char *obj_name);

        // Mesh Initialization

        /**
         * Rescale and center the mesh to bounding box.
         * TODO : Use the bounding box from the field.
         */
        void normalize_mesh();

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
         * Finds the edges, features and boundaries of the mesh.
         * TODO : Use the half edge iteration for boundaries. Use variance of normals for edges.
         */
        void find_edges_and_features_and_boundaries(bool should_preserve_boundaries, bool should_preserve_edges);

        /**
         * Initializes the mesh data structures.
         * TODO: Service level method to initialize datastructures
         */
        void initialize_parameterizer(
                bool should_preserve_boundaries,
                bool should_preserve_edges,
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
         * TODO : Service level method for the meshing flow. Split and simplify
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
                entities::DisjointTree &disajoint_tree,
                std::vector<int> &edge,
                std::vector<int> &face,
                std::vector<entities::DEdge> &edge_values,
                std::vector<Vector3i> &F2E,
                std::vector<Vector2i> &E2F,
                std::vector<Vector2i> &EdgeDiff,
                std::vector<Vector3i> &FQ
        );

        // scale
        void estimate_slope();

    };

}
