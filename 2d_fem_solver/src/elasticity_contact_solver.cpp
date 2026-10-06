#include "elasticity_contact_solver.h"

#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace pipes2d {

    namespace {

        struct ContactEdgeConstraint {
            int slave_n1 = -1;
            int slave_n2 = -1;
            int master_n1 = -1;
            int master_n2 = -1;
            double master_t = 0.5;
            vec2 normal = vec2::Zero();
            double initial_gap = 0.0;
            double slave_length = 0.0;
        };

        struct FixedDofGlobal {
            int dof = -1;
            double value = 0.0;
        };

        double clamp01(double x) {
            if (x < 0.0) return 0.0;
            if (x > 1.0) return 1.0;
            return x;
        }

        vec2 edge_midpoint(const Mesh& mesh, const BoundaryEdge& e) {
            const auto& a = mesh.nodes[e.n1];
            const auto& b = mesh.nodes[e.n2];
            return vec2(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
        }

        vec2 radial_normal_from_point(const vec2& p) {
            vec2 n = p;
            if (n.norm() == 0.0) {
                throw std::runtime_error("Контактная точка совпала с началом координат");
            }
            n.normalize();
            return n;
        }

        double point_segment_distance_sq(const vec2& P,
            const vec2& A,
            const vec2& B,
            double& t_out) {
            const vec2 AB = B - A;
            const double len2 = AB.squaredNorm();
            if (len2 <= 1e-30) {
                t_out = 0.0;
                return (P - A).squaredNorm();
            }
            const double t = clamp01((P - A).dot(AB) / len2);
            t_out = t;
            const vec2 Q = A + t * AB;
            return (P - Q).squaredNorm();
        }

        // Для каждого slave-ребра ищется ближайшее master-ребро по проекции середины slave-ребра.
        std::vector<ContactEdgeConstraint> build_contact_constraints_segment_to_segment(
            const Mesh& inner,
            int inner_type,
            const Mesh& outer,
            int outer_type) {

            const auto slave_edges = boundary_edges_of_type(inner, inner_type);
            const auto master_edges = boundary_edges_of_type(outer, outer_type);

            if (slave_edges.empty()) {
                throw std::runtime_error("Не найдены slave-ребра контактной границы");
            }
            if (master_edges.empty()) {
                throw std::runtime_error("Не найдены master-ребра контактной границы");
            }

            std::vector<ContactEdgeConstraint> constraints;
            constraints.reserve(slave_edges.size());

            for (const auto& se : slave_edges) {
                const vec2 P = edge_midpoint(inner, se);

                int best_idx = -1;
                double best_t = 0.0;
                double best_dist2 = std::numeric_limits<double>::max();

                for (int i = 0; i < static_cast<int>(master_edges.size()); ++i) {
                    const auto& me = master_edges[i];
                    const auto& a = outer.nodes[me.n1];
                    const auto& b = outer.nodes[me.n2];
                    const vec2 A(a.x, a.y);
                    const vec2 B(b.x, b.y);

                    double t = 0.0;
                    const double dist2 = point_segment_distance_sq(P, A, B, t);
                    if (dist2 < best_dist2) {
                        best_dist2 = dist2;
                        best_idx = i;
                        best_t = t;
                    }
                }

                if (best_idx < 0) {
                    throw std::runtime_error("Не удалось сопоставить slave-ребро master-ребру");
                }

                const auto& me = master_edges[best_idx];

                const auto& ma = outer.nodes[me.n1];
                const auto& mb = outer.nodes[me.n2];
                const vec2 A(ma.x, ma.y);
                const vec2 B(mb.x, mb.y);
                const vec2 Q = (1.0 - best_t) * A + best_t * B;

                ContactEdgeConstraint c;
                c.slave_n1 = se.n1;
                c.slave_n2 = se.n2;
                c.master_n1 = me.n1;
                c.master_n2 = me.n2;
                c.master_t = best_t;
                c.normal = radial_normal_from_point(P); // направление от внутренней трубы к внешней
                // Геометрический начальный зазор между серединой slave-ребра
                // и ее проекцией на master-ребро. Для согласованных окружностей
                // он должен быть близок к нулю, но не обязан быть точно нулевым
                // из-за дискретизации окружности ломаной.
                c.initial_gap = c.normal.dot(Q - P);
                c.slave_length = edge_length(inner, se);
                constraints.push_back(c);
            }

            return constraints;
        }

        std::map<int, double> collect_fixed_dofs_for_body(const BodyState& body, int global_offset) {
            std::map<int, double> fixed;

            // Старые закрепления по type границы.
            for (const auto& item : body.bc.fix_mask_by_type) {
                const int type = item.first;
                const auto& mask = item.second;
                vec2 value = vec2::Zero();
                const auto itv = body.bc.fix_value_by_type.find(type);
                if (itv != body.bc.fix_value_by_type.end()) {
                    value = itv->second;
                }

                const auto ids = unique_nodes_for_boundary_type(body.mesh, type);
                for (int node : ids) {
                    if (mask[0]) fixed[global_offset + 2 * node] = value[0];
                    if (mask[1]) fixed[global_offset + 2 * node + 1] = value[1];
                }
            }

            // Новые точечные закрепления.
            for (const auto& f : body.fixed_dofs) {
                if (f.node < 0 || f.node >= static_cast<int>(body.mesh.nodes.size())) {
                    throw std::runtime_error("Некорректный номер узла в fixed_dofs");
                }
                if (f.component != 0 && f.component != 1) {
                    throw std::runtime_error("Некорректная компонента в fixed_dofs");
                }
                fixed[global_offset + 2 * f.node + f.component] = f.value;
            }

            return fixed;
        }

        std::vector<FixedDofGlobal> collect_fixed_dofs(const BodyState& inner,
            const BodyState& outer,
            int n_in) {
            std::map<int, double> all = collect_fixed_dofs_for_body(inner, 0);
            const auto outer_fixed = collect_fixed_dofs_for_body(outer, n_in);
            for (const auto& item : outer_fixed) {
                all[item.first] = item.second;
            }

            std::vector<FixedDofGlobal> out;
            out.reserve(all.size());
            for (const auto& item : all) {
                out.push_back(FixedDofGlobal{ item.first, item.second });
            }
            return out;
        }

        // Решение СЛАУ с точным исключением закрепленных степеней свободы.
        // Это заменяет штраф 1e30 и резко улучшает обусловленность седловой системы.
        Eigen::VectorXd solve_sparse_system_with_fixed_dofs(
            const Eigen::SparseMatrix<double>& A,
            const Eigen::VectorXd& rhs,
            const std::vector<FixedDofGlobal>& fixed_dofs) {

            const int n = static_cast<int>(A.rows());
            if (A.cols() != n || rhs.size() != n) {
                throw std::runtime_error("Некорректный размер системы");
            }

            std::vector<char> is_fixed(n, 0);
            Eigen::VectorXd fixed_value = Eigen::VectorXd::Zero(n);
            for (const auto& f : fixed_dofs) {
                if (f.dof < 0 || f.dof >= n) {
                    throw std::runtime_error("Некорректный номер закрепленной степени свободы");
                }
                is_fixed[f.dof] = 1;
                fixed_value[f.dof] = f.value;
            }

            std::vector<int> full_to_reduced(n, -1);
            int n_reduced = 0;
            for (int i = 0; i < n; ++i) {
                if (!is_fixed[i]) {
                    full_to_reduced[i] = n_reduced++;
                }
            }

            Eigen::VectorXd rhs_reduced = Eigen::VectorXd::Zero(n_reduced);
            for (int i = 0; i < n; ++i) {
                if (!is_fixed[i]) {
                    rhs_reduced[full_to_reduced[i]] = rhs[i];
                }
            }

            std::vector<triplet> trips_reduced;
            trips_reduced.reserve(static_cast<std::size_t>(A.nonZeros()));

            for (int col = 0; col < A.outerSize(); ++col) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
                    const int i = it.row();
                    const int j = it.col();
                    const double v = it.value();

                    if (is_fixed[i]) {
                        continue;
                    }

                    const int ir = full_to_reduced[i];
                    if (is_fixed[j]) {
                        rhs_reduced[ir] -= v * fixed_value[j];
                    }
                    else {
                        trips_reduced.emplace_back(ir, full_to_reduced[j], v);
                    }
                }
            }

            Eigen::SparseMatrix<double> A_reduced(n_reduced, n_reduced);
            A_reduced.setFromTriplets(trips_reduced.begin(), trips_reduced.end());

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(A_reduced);
            solver.factorize(A_reduced);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Не удалось факторизовать редуцированную систему");
            }

            const Eigen::VectorXd x_reduced = solver.solve(rhs_reduced);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Не удалось решить редуцированную систему");
            }

            Eigen::VectorXd x = fixed_value;
            for (int i = 0; i < n; ++i) {
                if (!is_fixed[i]) {
                    x[i] = x_reduced[full_to_reduced[i]];
                }
            }
            return x;
        }

        vec2 displacement_at_node(const Eigen::VectorXd& u, int node) {
            return vec2(u[2 * node], u[2 * node + 1]);
        }

        double constraint_gap(const ContactEdgeConstraint& c,
            const Eigen::VectorXd& u_inner,
            const Eigen::VectorXd& u_outer) {

            const vec2 us1 = displacement_at_node(u_inner, c.slave_n1);
            const vec2 us2 = displacement_at_node(u_inner, c.slave_n2);
            const vec2 um1 = displacement_at_node(u_outer, c.master_n1);
            const vec2 um2 = displacement_at_node(u_outer, c.master_n2);

            const vec2 us_mid = 0.5 * (us1 + us2);
            const vec2 um_proj = (1.0 - c.master_t) * um1 + c.master_t * um2;

            // g > 0 — раскрытие, g < 0 — проникновение.
            return c.initial_gap + c.normal.dot(um_proj - us_mid);
        }

    } // namespace

    std::vector<ContactPair> CoupledThermoMechanicalContact2D::build_contact_pairs(
        const Mesh& inner,
        int inner_type,
        const Mesh& outer,
        int outer_type) const {

        const auto constraints = build_contact_constraints_segment_to_segment(
            inner, inner_type, outer, outer_type);

        std::vector<ContactPair> pairs;
        pairs.reserve(constraints.size());
        for (const auto& c : constraints) {
            pairs.push_back(ContactPair{ c.slave_n1, c.master_n1, c.normal, c.initial_gap });
        }
        return pairs;
    }

    std::vector<triplet> CoupledThermoMechanicalContact2D::assemble_body(
        const BodyState& body,
        Eigen::VectorXd& rhs) const {

        const int ndofs = static_cast<int>(body.mesh.nodes.size()) * 2;
        rhs = Eigen::VectorXd::Zero(ndofs);

        std::vector<triplet> trips;
        const mat33 C = elasticity_matrix(body.material);

        for (const auto& tri : body.mesh.elements) {
            const double A = triangle_area(body.mesh, tri);
            const mat36 B = make_B(body.mesh, tri);
            const mat66 Ke = A * (B.transpose() * C * B);

            double dT = 0.0;
            if (body.thermal.temperature.size() == static_cast<int>(body.mesh.nodes.size())) {
                dT = (body.thermal.temperature[tri.node[0]]
                    + body.thermal.temperature[tri.node[1]]
                    + body.thermal.temperature[tri.node[2]]) / 3.0
                    - body.reference_temperature;
            }

            const vec3 epsT = thermal_strain_vector(body.material, dT);
            const Eigen::Matrix<double, 6, 1> fe = A * (B.transpose() * C * epsT);

            std::array<int, 6> dof{};
            for (int i = 0; i < 3; ++i) {
                dof[2 * i] = 2 * tri.node[i];
                dof[2 * i + 1] = 2 * tri.node[i] + 1;
            }

            for (int i = 0; i < 6; ++i) {
                rhs[dof[i]] += fe[i];
                for (int j = 0; j < 6; ++j) {
                    trips.emplace_back(dof[i], dof[j], Ke(i, j));
                }
            }
        }

        for (const auto& edge : body.mesh.boundary) {
            auto it = body.bc.traction_by_type.find(edge.type);
            if (it == body.bc.traction_by_type.end()) {
                continue;
            }

            const vec2 t = it->second;
            const double L = edge_length(body.mesh, edge);
            Eigen::Matrix<double, 4, 1> fe;
            fe << t[0], t[1], t[0], t[1];
            fe *= 0.5 * L;

            const std::array<int, 4> dof = {
                2 * edge.n1, 2 * edge.n1 + 1,
                2 * edge.n2, 2 * edge.n2 + 1
            };

            for (int i = 0; i < 4; ++i) {
                rhs[dof[i]] += fe[i];
            }
        }

        return trips;
    }

    double CoupledThermoMechanicalContact2D::radial_displacement(
        const Mesh& mesh,
        int node,
        const Eigen::VectorXd& u) const {

        const auto& p = mesh.nodes[node];
        vec2 n(p.x, p.y);
        if (n.norm() == 0.0) {
            return 0.0;
        }
        n.normalize();
        return n.dot(vec2(u[2 * node], u[2 * node + 1]));
    }

    CoupledResult CoupledThermoMechanicalContact2D::solve(
        const BodyState& inner,
        const BodyState& outer,
        int inner_contact_type,
        int outer_contact_type,
        double contact_mode_parameter) const {

        CoupledResult out;

        const auto constraints = build_contact_constraints_segment_to_segment(
            inner.mesh, inner_contact_type,
            outer.mesh, outer_contact_type);

        out.pairs.reserve(constraints.size());
        out.contact_edges.reserve(constraints.size());
        for (const auto& c : constraints) {
            out.pairs.push_back(ContactPair{ c.slave_n1, c.master_n1, c.normal, c.initial_gap });
            out.contact_edges.push_back(ContactEdgeInfo{
                c.slave_n1, c.slave_n2,
                c.master_n1, c.master_n2,
                c.master_t,
                c.slave_length,
                c.normal
                });
        }

        const int n_in = static_cast<int>(inner.mesh.nodes.size()) * 2;
        const int n_out = static_cast<int>(outer.mesh.nodes.size()) * 2;
        const int n_all = n_in + n_out;

        Eigen::VectorXd rhs_in, rhs_out;
        auto trips = assemble_body(inner, rhs_in);
        const auto trips_out = assemble_body(outer, rhs_out);

        Eigen::VectorXd rhs_mech = Eigen::VectorXd::Zero(n_all);
        rhs_mech.segment(0, n_in) = rhs_in;
        rhs_mech.segment(n_in, n_out) = rhs_out;

        for (const auto& t : trips_out) {
            trips.emplace_back(t.row() + n_in, t.col() + n_in, t.value());
        }

        const std::vector<FixedDofGlobal> fixed_dofs = collect_fixed_dofs(inner, outer, n_in);

        // Механическая задача без контакта.
        if (contact_mode_parameter <= 0.0) {
            Eigen::SparseMatrix<double> K(n_all, n_all);
            K.setFromTriplets(trips.begin(), trips.end());

            const Eigen::VectorXd u = solve_sparse_system_with_fixed_dofs(K, rhs_mech, fixed_dofs);

            out.inner.displacement = u.segment(0, n_in);
            out.outer.displacement = u.segment(n_in, n_out);
            out.contact_pressure.assign(constraints.size(), 0.0);
            out.contact_gap.reserve(constraints.size());

            for (const auto& c : constraints) {
                out.contact_gap.push_back(constraint_gap(c, out.inner.displacement, out.outer.displacement));
            }
            return out;
        }

        // Полный контакт через рёберные множители Лагранжа.
        // Все контактные рёбра считаются активными. Это соответствует текущей
        // постановке задачи: контакт задан на всей известной границе Gamma_c.
        const int n_contact = static_cast<int>(constraints.size());
        const int n_total = n_all + n_contact;

        std::vector<triplet> trips_total = trips;
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_total);
        rhs.segment(0, n_all) = rhs_mech;

        for (int i = 0; i < n_contact; ++i) {
            const auto& c = constraints[i];
            const int row = n_all + i;

            const vec2 n = c.normal;
            const double L = c.slave_length;
            const double alpha = 1.0 - c.master_t;
            const double beta = c.master_t;

            const int is1x = 2 * c.slave_n1;
            const int is1y = 2 * c.slave_n1 + 1;
            const int is2x = 2 * c.slave_n2;
            const int is2y = 2 * c.slave_n2 + 1;

            const int im1x = n_in + 2 * c.master_n1;
            const int im1y = n_in + 2 * c.master_n1 + 1;
            const int im2x = n_in + 2 * c.master_n2;
            const int im2y = n_in + 2 * c.master_n2 + 1;

            // Зазор принят как
            //     g = g0 + n · (u_master - u_slave).
            // Полный контакт означает g = 0.
            //
            // В строке ограничения используем нормализованную форму g = 0,
            // без умножения на длину ребра. Это улучшает численное выполнение
            // контактного условия, потому что строки не становятся слишком малыми.
            //
            // В столбце множителя длина ребра сохраняется: множитель lambda
            // по-прежнему интерпретируется как контактное давление, а его вклад
            // в уравнения равновесия является интегралом давления по ребру.

            const double slave_shape = 0.5;
            const double master_shape_1 = alpha;
            const double master_shape_2 = beta;

            // Строка ограничения g = 0.
            // Slave midpoint: -0.5 * u_s1 - 0.5 * u_s2.
            trips_total.emplace_back(row, is1x, -slave_shape * n[0]);
            trips_total.emplace_back(row, is1y, -slave_shape * n[1]);
            trips_total.emplace_back(row, is2x, -slave_shape * n[0]);
            trips_total.emplace_back(row, is2y, -slave_shape * n[1]);

            // Master projection: alpha * u_m1 + beta * u_m2.
            trips_total.emplace_back(row, im1x, master_shape_1 * n[0]);
            trips_total.emplace_back(row, im1y, master_shape_1 * n[1]);
            trips_total.emplace_back(row, im2x, master_shape_2 * n[0]);
            trips_total.emplace_back(row, im2y, master_shape_2 * n[1]);

            // Столбец множителя: контактные силы от давления.
            const double cs = 0.5 * L;
            const double cm1 = L * alpha;
            const double cm2 = L * beta;

            // На внутреннюю трубу действует -lambda*n.
            trips_total.emplace_back(is1x, row, -cs * n[0]);
            trips_total.emplace_back(is1y, row, -cs * n[1]);
            trips_total.emplace_back(is2x, row, -cs * n[0]);
            trips_total.emplace_back(is2y, row, -cs * n[1]);

            // На внешнюю трубу действует +lambda*n.
            trips_total.emplace_back(im1x, row, cm1 * n[0]);
            trips_total.emplace_back(im1y, row, cm1 * n[1]);
            trips_total.emplace_back(im2x, row, cm2 * n[0]);
            trips_total.emplace_back(im2y, row, cm2 * n[1]);

            rhs[row] = -c.initial_gap;
        }

        Eigen::SparseMatrix<double> A(n_total, n_total);
        A.setFromTriplets(trips_total.begin(), trips_total.end());

        const Eigen::VectorXd sol = solve_sparse_system_with_fixed_dofs(A, rhs, fixed_dofs);
        const Eigen::VectorXd solution_all = sol.segment(0, n_all);

        out.inner.displacement = solution_all.segment(0, n_in);
        out.outer.displacement = solution_all.segment(n_in, n_out);

        out.contact_pressure.assign(n_contact, 0.0);
        out.contact_gap.assign(n_contact, 0.0);

        for (int i = 0; i < n_contact; ++i) {
            out.contact_pressure[i] = sol[n_all + i];
            out.contact_gap[i] = constraint_gap(constraints[i], out.inner.displacement, out.outer.displacement);
        }

        return out;
    }

} // namespace pipes2d
