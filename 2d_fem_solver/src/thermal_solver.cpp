#include "thermal_solver.h"

#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace pipes2d {

    namespace {

        struct DirichletData {
            std::vector<char> is_dirichlet;
            Eigen::VectorXd value;
        };

        /// Собрать информацию о граничных условиях Дирихле.
        DirichletData build_dirichlet_data(const Mesh& mesh, const BoundaryCondition2D& bc) {
            const int n = static_cast<int>(mesh.nodes.size());

            DirichletData out;
            out.is_dirichlet.assign(n, 0);
            out.value = Eigen::VectorXd::Zero(n);

            for (const auto& item : bc.temperature_by_type) {
                const int type = item.first;
                const double val = item.second;
                const auto ids = unique_nodes_for_boundary_type(mesh, type);
                for (int id : ids) {
                    out.is_dirichlet[id] = 1;
                    out.value[id] = val;
                }
            }
            return out;
        }

        /// Собрать матрицу теплопроводности и диагональную матрицу теплоемкости.
        void assemble_thermal_matrices(
            const Mesh& mesh,
            const Material& material,
            Eigen::SparseMatrix<double>& Kfull,
            Eigen::VectorXd& Mdiag) {

            const int n = static_cast<int>(mesh.nodes.size());
            std::vector<triplet> Ktrips;
            Ktrips.reserve(mesh.elements.size() * 9);
            Mdiag = Eigen::VectorXd::Zero(n);

            for (const auto& tri : mesh.elements) {
                const double A = triangle_area(mesh, tri);
                const auto g = grad_shape(mesh, tri);
                const Eigen::Matrix3d k = material.k * A * (g.transpose() * g);

                for (int i = 0; i < 3; ++i) {
                    const int gi = tri.node[i];
                    Mdiag[gi] += material.rho * material.cp * A / 3.0;
                    for (int j = 0; j < 3; ++j) {
                        Ktrips.emplace_back(gi, tri.node[j], k(i, j));
                    }
                }
            }

            Kfull.resize(n, n);
            Kfull.setFromTriplets(Ktrips.begin(), Ktrips.end());
        }

        /// Один неявный шаг Эйлера.
        Eigen::VectorXd implicit_euler_step(
            const Eigen::SparseMatrix<double>& Kfull,
            const Eigen::VectorXd& Mdiag,
            const DirichletData& dir,
            const Eigen::VectorXd& Tprev,
            double dt) {

            const int n = static_cast<int>(Tprev.size());
            std::vector<triplet> Atrips;
            Atrips.reserve(Kfull.nonZeros() + n);
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

            for (int k = 0; k < Kfull.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(Kfull, k); it; ++it) {
                    const int i = static_cast<int>(it.row());
                    const int j = static_cast<int>(it.col());
                    const double aij = it.value();

                    if (dir.is_dirichlet[i]) {
                        continue;
                    }
                    else if (dir.is_dirichlet[j]) {
                        rhs[i] -= aij * dir.value[j];
                    }
                    else {
                        Atrips.emplace_back(i, j, aij);
                    }
                }
            }

            for (int i = 0; i < n; ++i) {
                if (dir.is_dirichlet[i]) {
                    Atrips.emplace_back(i, i, 1.0);
                    rhs[i] = dir.value[i];
                }
                else {
                    Atrips.emplace_back(i, i, Mdiag[i] / dt);
                    rhs[i] += Mdiag[i] / dt * Tprev[i];
                }
            }

            Eigen::SparseMatrix<double> A(n, n);
            A.setFromTriplets(Atrips.begin(), Atrips.end());

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Не удалось факторизовать матрицу нестационарной тепловой задачи");
            }

            const Eigen::VectorXd Tnew = solver.solve(rhs);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Не удалось решить нестационарную тепловую задачу");
            }
            return Tnew;
        }

        double max_temperature_error(const Eigen::VectorXd& T, double target) {
            double max_error = 0.0;
            for (int i = 0; i < T.size(); ++i) {
                max_error = (std::max)(max_error, std::abs(T[i] - target));
            }
            return max_error;
        }

        void save_temperature_history(ThermalResult& out,
            const Eigen::VectorXd& T,
            double time,
            int history_stride,
            int step,
            bool force) {

            if (history_stride <= 0) {
                return;
            }

            if (force || step == 0 || step % history_stride == 0) {
                out.history.push_back(T);
                out.history_time.push_back(time);
            }
        }

    } // namespace

    ThermalSolver2D::ThermalSolver2D(const Mesh& mesh, const Material& material)
        : mesh_(mesh), material_(material) {}

    ThermalResult ThermalSolver2D::solve_steady(
        const BoundaryCondition2D& bc,
        const std::function<double(double, double)>& source) const {

        const int n = static_cast<int>(mesh_.nodes.size());
        const auto dir = build_dirichlet_data(mesh_, bc);

        std::vector<triplet> trips;
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
        trips.reserve(mesh_.elements.size() * 9 + n);

        for (const auto& tri : mesh_.elements) {
            const double A = triangle_area(mesh_, tri);
            const auto g = grad_shape(mesh_, tri);
            const Eigen::Matrix3d k = material_.k * A * (g.transpose() * g);

            Eigen::Vector3d f = Eigen::Vector3d::Zero();
            if (source) {
                const auto& a = mesh_.nodes[tri.node[0]];
                const auto& b = mesh_.nodes[tri.node[1]];
                const auto& c = mesh_.nodes[tri.node[2]];
                const double xc = (a.x + b.x + c.x) / 3.0;
                const double yc = (a.y + b.y + c.y) / 3.0;
                f.setConstant(source(xc, yc) * A / 3.0);
            }

            for (int i = 0; i < 3; ++i) {
                const int gi = tri.node[i];
                if (!dir.is_dirichlet[gi]) {
                    rhs[gi] += f[i];
                }

                for (int j = 0; j < 3; ++j) {
                    const int gj = tri.node[j];
                    const double aij = k(i, j);

                    if (dir.is_dirichlet[gi]) {
                        continue;
                    }
                    else if (dir.is_dirichlet[gj]) {
                        rhs[gi] -= aij * dir.value[gj];
                    }
                    else {
                        trips.emplace_back(gi, gj, aij);
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            if (dir.is_dirichlet[i]) {
                trips.emplace_back(i, i, 1.0);
                rhs[i] = dir.value[i];
            }
        }

        Eigen::SparseMatrix<double> K(n, n);
        K.setFromTriplets(trips.begin(), trips.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(K);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Не удалось факторизовать матрицу стационарной тепловой задачи");
        }

        ThermalResult out;
        out.temperature = solver.solve(rhs);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Не удалось решить стационарную тепловую задачу");
        }
        return out;
    }

    ThermalResult ThermalSolver2D::solve_transient_implicit(
        const BoundaryCondition2D& bc,
        double initial_temperature,
        double dt,
        int steps,
        const std::function<double(double, double, double)>& source) const {

        if (source) {
            throw std::runtime_error("В текущей версии solve_transient_implicit не поддерживает внутренний источник q(x,y,t)");
        }
        if (dt <= 0.0 || steps <= 0) {
            throw std::runtime_error("Некорректные параметры нестационарной тепловой задачи");
        }

        const int n = static_cast<int>(mesh_.nodes.size());
        const auto dir = build_dirichlet_data(mesh_, bc);

        Eigen::SparseMatrix<double> Kfull;
        Eigen::VectorXd Mdiag;
        assemble_thermal_matrices(mesh_, material_, Kfull, Mdiag);

        Eigen::VectorXd Tprev = Eigen::VectorXd::Constant(n, initial_temperature);
        for (int i = 0; i < n; ++i) {
            if (dir.is_dirichlet[i]) {
                Tprev[i] = dir.value[i];
            }
        }

        for (int step = 0; step < steps; ++step) {
            Tprev = implicit_euler_step(Kfull, Mdiag, dir, Tprev, dt);
        }

        ThermalResult out;
        out.temperature = Tprev;
        out.performed_steps = steps;
        out.final_time = dt * static_cast<double>(steps);
        return out;
    }

    ThermalResult ThermalSolver2D::solve_unsteady_cooling(
        int cold_boundary_type,
        double cold_temperature,
        double dt,
        int steps,
        double initial_temperature) const {

        BoundaryCondition2D bc;
        bc.temperature_by_type[cold_boundary_type] = cold_temperature;
        return solve_transient_implicit(bc, initial_temperature, dt, steps);
    }

    ThermalResult ThermalSolver2D::solve_unsteady_cooling_until_converged(
        int cold_boundary_type,
        double cold_temperature,
        double dt,
        int max_steps,
        double initial_temperature,
        double temperature_tolerance,
        int history_stride) const {

        if (dt <= 0.0) {
            throw std::runtime_error("Шаг по времени должен быть положительным");
        }
        if (max_steps <= 0) {
            throw std::runtime_error("Максимальное число шагов должно быть положительным");
        }
        if (temperature_tolerance < 0.0) {
            throw std::runtime_error("Температурный критерий останова не может быть отрицательным");
        }

        BoundaryCondition2D bc;
        bc.temperature_by_type[cold_boundary_type] = cold_temperature;

        const int n = static_cast<int>(mesh_.nodes.size());
        const auto dir = build_dirichlet_data(mesh_, bc);

        Eigen::SparseMatrix<double> Kfull;
        Eigen::VectorXd Mdiag;
        assemble_thermal_matrices(mesh_, material_, Kfull, Mdiag);

        ThermalResult out;
        Eigen::VectorXd Tprev = Eigen::VectorXd::Constant(n, initial_temperature);
        for (int i = 0; i < n; ++i) {
            if (dir.is_dirichlet[i]) {
                Tprev[i] = dir.value[i];
            }
        }

        save_temperature_history(out, Tprev, 0.0, history_stride, 0, true);

        for (int step = 1; step <= max_steps; ++step) {
            const Eigen::VectorXd Tnew = implicit_euler_step(Kfull, Mdiag, dir, Tprev, dt);
            const double max_error = max_temperature_error(Tnew, cold_temperature);
            const double time = dt * static_cast<double>(step);

            save_temperature_history(out, Tnew, time, history_stride, step, false);

            Tprev = Tnew;
            if (max_error <= temperature_tolerance) {
                out.temperature = Tprev;
                out.performed_steps = step;
                out.final_time = time;
                out.final_max_error = max_error;
                save_temperature_history(out, Tprev, time, history_stride, step, true);
                return out;
            }
        }

        throw std::runtime_error(
            "Нестационарная тепловая задача не достигла критерия max|T - Tc| за заданное число шагов");
    }

} // namespace pipes2d
