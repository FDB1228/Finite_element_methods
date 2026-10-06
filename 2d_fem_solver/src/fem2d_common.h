#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pipes2d {

    using vec2 = Eigen::Vector2d;
    using vec3 = Eigen::Vector3d;
    using mat33 = Eigen::Matrix3d;
    using mat36 = Eigen::Matrix<double, 3, 6>;
    using mat66 = Eigen::Matrix<double, 6, 6>;
    using triplet = Eigen::Triplet<double>;

    /// Узел двумерной сетки.
    struct Node {
        double x = 0.0;
        double y = 0.0;
    };

    /// Линейный треугольный конечный элемент.
    struct Triangle {
        std::array<int, 3> node{};
        std::array<int, 3> neigh{};
        int subdomain = 0;
    };

    /// Граничное ребро.
    struct BoundaryEdge {
        int n1 = -1;
        int n2 = -1;
        int owner = -1;
        int type = 0;
    };

    /// Конечно-элементная сетка.
    struct Mesh {
        std::vector<Node> nodes;
        std::vector<Triangle> elements;
        std::vector<BoundaryEdge> boundary;
        std::vector<BoundaryEdge> constraints;
    };

    /// Параметры материала.
    struct Material {
        double E = 2.0e11;
        double nu = 0.33;
        double alpha = 1.0e-5;
        double k = 30.0;
        double rho = 1.0e4;
        double cp = 300.0;
        bool plane_stress = false;
    };

    /// Контактная пара узлов.
    struct ContactPair {
        int slave = -1;
        int master = -1;
        vec2 normal = vec2::Zero();
        double initial_gap = 0.0;
    };

    /// Граничные условия.
    struct BoundaryCondition2D {
        std::map<int, double> temperature_by_type;
        std::set<int> insulated_types;
        std::map<int, vec2> traction_by_type;
        std::map<int, std::array<bool, 2>> fix_mask_by_type;
        std::map<int, vec2> fix_value_by_type;
    };

    /// Площадь треугольника.
    inline double triangle_area(const Mesh& mesh, const Triangle& tri) {
        const auto& a = mesh.nodes[tri.node[0]];
        const auto& b = mesh.nodes[tri.node[1]];
        const auto& c = mesh.nodes[tri.node[2]];
        const double det = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        return 0.5 * std::abs(det);
    }

    /// Градиенты линейных функций формы.
    inline Eigen::Matrix<double, 2, 3> grad_shape(const Mesh& mesh, const Triangle& tri) {
        const auto& a = mesh.nodes[tri.node[0]];
        const auto& b = mesh.nodes[tri.node[1]];
        const auto& c = mesh.nodes[tri.node[2]];
        const double twoA = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);

        if (std::abs(twoA) < 1e-20) {
            throw std::runtime_error("Обнаружен вырожденный треугольник");
        }

        Eigen::Matrix<double, 2, 3> g;
        g(0, 0) = (b.y - c.y) / twoA;
        g(1, 0) = (c.x - b.x) / twoA;
        g(0, 1) = (c.y - a.y) / twoA;
        g(1, 1) = (a.x - c.x) / twoA;
        g(0, 2) = (a.y - b.y) / twoA;
        g(1, 2) = (b.x - a.x) / twoA;
        return g;
    }

    /// B-матрица для плоской задачи упругости.
    inline mat36 make_B(const Mesh& mesh, const Triangle& tri) {
        const auto g = grad_shape(mesh, tri);
        mat36 B = mat36::Zero();
        for (int i = 0; i < 3; ++i) {
            B(0, 2 * i) = g(0, i);
            B(1, 2 * i + 1) = g(1, i);
            B(2, 2 * i) = g(1, i);
            B(2, 2 * i + 1) = g(0, i);
        }
        return B;
    }

    /// Матрица упругости для плоской деформации или плоского напряженного состояния.
    inline mat33 elasticity_matrix(const Material& m) {
        const double E = m.E;
        const double nu = m.nu;
        mat33 C = mat33::Zero();

        if (m.plane_stress) {
            const double c = E / (1.0 - nu * nu);
            C << c, c* nu, 0.0,
                c* nu, c, 0.0,
                0.0, 0.0, c* (1.0 - nu) * 0.5;
        }
        else {
            const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            C << c * (1.0 - nu), c* nu, 0.0,
                c* nu, c* (1.0 - nu), 0.0,
                0.0, 0.0, c* (1.0 - 2.0 * nu) * 0.5;
        }
        return C;
    }

    /// Вектор температурных деформаций.
    inline vec3 thermal_strain_vector(const Material& m, double dT) {
        const double a = m.alpha * dT;
        vec3 epsT;
        epsT << a, a, 0.0;
        return epsT;
    }

    /// Длина ребра.
    inline double edge_length(const Mesh& mesh, const BoundaryEdge& e) {
        const auto& a = mesh.nodes[e.n1];
        const auto& b = mesh.nodes[e.n2];
        return std::hypot(b.x - a.x, b.y - a.y);
    }

    /// Радиус точки относительно начала координат.
    inline double radius(const Node& p) {
        return std::hypot(p.x, p.y);
    }

    /// Уникальные узлы граничных ребер заданного типа.
    inline std::vector<int> unique_nodes_for_boundary_type(const Mesh& mesh, int type) {
        std::set<int> uniq;
        for (const auto& e : mesh.boundary) {
            if (e.type == type) {
                uniq.insert(e.n1);
                uniq.insert(e.n2);
            }
        }
        return std::vector<int>(uniq.begin(), uniq.end());
    }

    /// Все граничные ребра заданного типа.
    inline std::vector<BoundaryEdge> boundary_edges_of_type(const Mesh& mesh, int type) {
        std::vector<BoundaryEdge> out;
        for (const auto& e : mesh.boundary) {
            if (e.type == type) {
                out.push_back(e);
            }
        }
        return out;
    }

    /// Сортировка узлов по углу вокруг начала координат.
    inline std::vector<int> sort_nodes_by_angle(const Mesh& mesh, const std::vector<int>& ids) {
        std::vector<int> out = ids;
        std::sort(out.begin(), out.end(), [&](int a, int b) {
            const auto& A = mesh.nodes[a];
            const auto& B = mesh.nodes[b];
            return std::atan2(A.y, A.x) < std::atan2(B.y, B.x);
            });
        return out;
    }

    /// Вес узла на границе: половина суммы длин инцидентных ребер заданного типа.
    inline std::map<int, double> nodal_boundary_weights(const Mesh& mesh, int boundary_type) {
        const auto edges = boundary_edges_of_type(mesh, boundary_type);
        if (edges.empty()) {
            throw std::runtime_error("Не найдены граничные ребра для вычисления весов узлов");
        }

        std::map<int, double> w;
        for (const auto& e : edges) {
            const double L = edge_length(mesh, e);
            w[e.n1] += 0.5 * L;
            w[e.n2] += 0.5 * L;
        }
        return w;
    }

} // namespace pipes2d
