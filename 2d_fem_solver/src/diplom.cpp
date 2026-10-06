#include "elasticity_contact_solver.h"
#include "gridder2d_io.h"
#include "thermal_solver.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace pipes2d;

namespace {

    struct GapStats {
        double minimum = 0.0;
        double maximum = 0.0;
        double average = 0.0;
    };

    struct PressureStats {
        double minimum_signed = 0.0;
        double maximum_signed = 0.0;
        double average_signed = 0.0;
        double average_physical = 0.0;
    };

    struct ComparisonStats {
        double relative_to_delta_h_percent = 0.0;
        double relative_to_delta_an_percent = 0.0;
    };

    struct TemperatureStats {
        double minimum = 0.0;
        double maximum = 0.0;
        double average = 0.0;
    };

    struct NodalStressField {
        std::vector<double> sigma_xx;
        std::vector<double> sigma_yy;
        std::vector<double> tau_xy;
        std::vector<double> sigma_rr;
        std::vector<double> sigma_tt;
        std::vector<double> tau_rt;
        std::vector<double> von_mises;
    };

    // Масштаб изображения деформированной формы в Tecplot.
    // Для реальной деформации поставить 1.0.
    const double DEFORMATION_SCALE_FOR_TECPLOT = 50.0;

    // Сохранять каждый N-й шаг температурной истории для Tecplot-анимации.
    // Если поставить 0, анимация не будет записываться.
    const int THERMAL_HISTORY_STRIDE = 1;

    std::string join_path(const std::string& folder, const std::string& filename) {
        if (folder.empty()) {
            return filename;
        }
        const char last = folder[folder.size() - 1];
        if (last == '/' || last == '\\') {
            return folder + filename;
        }
        return folder + "/" + filename;
    }

    bool file_exists(const std::string& filename) {
        std::ifstream in(filename.c_str());
        return static_cast<bool>(in);
    }

    std::string select_grid_folder(int argc, char* argv[]) {
        const std::string default_grid_folder = "grid_0125";
        if (argc < 2) {
            return default_grid_folder;
        }

        const std::string arg = argv[1];
        if (arg.empty()) {
            return default_grid_folder;
        }

        const bool already_grid_name = (arg.compare(0, 5, "grid_") == 0);
        const bool looks_like_path =
            (arg.find('/') != std::string::npos) ||
            (arg.find('\\') != std::string::npos);

        if (already_grid_name || looks_like_path) {
            return arg;
        }
        return "grid_" + arg;
    }

    vec2 node_position(const Mesh& mesh, int node) {
        const Node& p = mesh.nodes[node];
        return vec2(p.x, p.y);
    }

    vec2 displacement_at_node(const Eigen::VectorXd& u, int node) {
        return vec2(u[2 * node], u[2 * node + 1]);
    }

    double radial_displacement_value(const Mesh& mesh, int node, const Eigen::VectorXd& u) {
        vec2 n = node_position(mesh, node);
        if (n.norm() == 0.0) {
            return 0.0;
        }
        n.normalize();
        return n.dot(displacement_at_node(u, node));
    }

    TemperatureStats temperature_stats(const Eigen::VectorXd& T) {
        TemperatureStats out;
        if (T.size() == 0) {
            return out;
        }

        out.minimum = T[0];
        out.maximum = T[0];
        out.average = 0.0;

        for (int i = 0; i < T.size(); ++i) {
            out.minimum = (std::min)(out.minimum, T[i]);
            out.maximum = (std::max)(out.maximum, T[i]);
            out.average += T[i];
        }
        out.average /= static_cast<double>(T.size());
        return out;
    }

    GapStats gap_stats(const std::vector<double>& gaps) {
        GapStats out;
        if (gaps.empty()) {
            return out;
        }

        out.minimum = std::numeric_limits<double>::max();
        out.maximum = -std::numeric_limits<double>::max();

        for (double g : gaps) {
            out.minimum = (std::min)(out.minimum, g);
            out.maximum = (std::max)(out.maximum, g);
            out.average += g;
        }
        out.average /= static_cast<double>(gaps.size());
        return out;
    }

    PressureStats pressure_stats_by_edges(const CoupledResult& result) {
        PressureStats out;
        if (result.contact_pressure.empty() || result.contact_edges.empty()) {
            return out;
        }

        out.minimum_signed = std::numeric_limits<double>::max();
        out.maximum_signed = -std::numeric_limits<double>::max();

        double sum_signed = 0.0;
        double total_length = 0.0;

        const std::size_t n = (std::min)(result.contact_pressure.size(), result.contact_edges.size());
        for (std::size_t i = 0; i < n; ++i) {
            const double L = result.contact_edges[i].length;
            if (L <= 0.0) {
                continue;
            }

            const double p = result.contact_pressure[i];
            out.minimum_signed = (std::min)(out.minimum_signed, p);
            out.maximum_signed = (std::max)(out.maximum_signed, p);
            sum_signed += L * p;
            total_length += L;
        }

        if (total_length > 0.0) {
            out.average_signed = sum_signed / total_length;
            out.average_physical = -out.average_signed;
        }

        if (out.minimum_signed == std::numeric_limits<double>::max()) {
            out.minimum_signed = 0.0;
            out.maximum_signed = 0.0;
        }
        return out;
    }

    ComparisonStats make_comparison(double p_num,
        double p_ref_delta_h,
        double p_ref_delta_an) {

        ComparisonStats out;
        if (std::abs(p_ref_delta_h) > 0.0) {
            out.relative_to_delta_h_percent =
                100.0 * std::abs(p_num - p_ref_delta_h) / std::abs(p_ref_delta_h);
        }
        if (std::abs(p_ref_delta_an) > 0.0) {
            out.relative_to_delta_an_percent =
                100.0 * std::abs(p_num - p_ref_delta_an) / std::abs(p_ref_delta_an);
        }
        return out;
    }

    double explicit_contact_pressure_formula(double delta,
        double a,
        double c,
        double b,
        const Material& inner_material,
        const Material& outer_material) {

        const double E1 = inner_material.E;
        const double E2 = outer_material.E;
        const double nu1 = inner_material.nu;
        const double nu2 = outer_material.nu;

        const double compliance_inner = c / E1 *
            ((1.0 - nu1) * c * c + (1.0 + nu1) * a * a) /
            (c * c - a * a);

        const double compliance_outer = c / E2 *
            ((1.0 - nu2) * c * c + (1.0 + nu2) * b * b) /
            (b * b - c * c);

        const double compliance = compliance_inner + compliance_outer;
        if (compliance <= 0.0) {
            throw std::runtime_error("Некорректная податливость в аналитической формуле давления");
        }
        return delta / compliance;
    }

    const char* plane_model_name(const Material& m) {
        return m.plane_stress ? "плоское напряженное состояние" : "плоское деформированное состояние";
    }

    bool is_on_phi0_ray(const Node& p) {
        const double r = std::hypot(p.x, p.y);
        return r > 1.0e-14 && p.x > 0.0 && std::abs(p.y) <= 1.0e-9 * r + 1.0e-14;
    }

    bool is_on_phi_half_pi_ray(const Node& p) {
        const double r = std::hypot(p.x, p.y);
        return r > 1.0e-14 && p.y > 0.0 && std::abs(p.x) <= 1.0e-9 * r + 1.0e-14;
    }

    void add_unique_fixed_dof(BodyState& body, int node, int component, double value) {
        for (const auto& f : body.fixed_dofs) {
            if (f.node == node && f.component == component) {
                return;
            }
        }
        body.fixed_dofs.push_back(FixedDof2D{ node, component, value });
    }

    /// Условия симметрии для четверти области:
    /// phi = 0: uy = 0; phi = pi/2: ux = 0.
    void apply_quarter_symmetry_constraints(BodyState& body) {
        int count_phi0 = 0;
        int count_phi90 = 0;

        for (int i = 0; i < static_cast<int>(body.mesh.nodes.size()); ++i) {
            const Node& p = body.mesh.nodes[i];

            if (is_on_phi0_ray(p)) {
                add_unique_fixed_dof(body, i, 1, 0.0);
                ++count_phi0;
            }
            if (is_on_phi_half_pi_ray(p)) {
                add_unique_fixed_dof(body, i, 0, 0.0);
                ++count_phi90;
            }
        }

        if (count_phi0 == 0) {
            throw std::runtime_error("Не найдены узлы на луче phi = 0 для условия uy = 0");
        }
        if (count_phi90 == 0) {
            throw std::runtime_error("Не найдены узлы на луче phi = pi/2 для условия ux = 0");
        }
    }

    double weighted_average_radial_displacement(
        const Mesh& mesh,
        int boundary_type,
        const Eigen::VectorXd& u) {

        const auto weights = nodal_boundary_weights(mesh, boundary_type);
        double sum_w = 0.0;
        double sum_u = 0.0;

        for (const auto& item : weights) {
            const int node = item.first;
            const double w = item.second;
            sum_w += w;
            sum_u += w * radial_displacement_value(mesh, node, u);
        }

        if (sum_w <= 0.0) {
            throw std::runtime_error("Нулевая длина контактной границы при усреднении перемещений");
        }
        return sum_u / sum_w;
    }

    /// Среднее по элементам напряжение переносится в узлы с весом площади элемента.
    NodalStressField compute_nodal_stress_field(
        const Mesh& mesh,
        const Material& material,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& T,
        double reference_temperature) {

        const int n_nodes = static_cast<int>(mesh.nodes.size());
        NodalStressField s;
        s.sigma_xx.assign(n_nodes, 0.0);
        s.sigma_yy.assign(n_nodes, 0.0);
        s.tau_xy.assign(n_nodes, 0.0);
        s.sigma_rr.assign(n_nodes, 0.0);
        s.sigma_tt.assign(n_nodes, 0.0);
        s.tau_rt.assign(n_nodes, 0.0);
        s.von_mises.assign(n_nodes, 0.0);

        std::vector<double> w(n_nodes, 0.0);
        const mat33 C = elasticity_matrix(material);

        for (const auto& tri : mesh.elements) {
            const double A = triangle_area(mesh, tri);
            const mat36 B = make_B(mesh, tri);

            Eigen::Matrix<double, 6, 1> ue;
            for (int i = 0; i < 3; ++i) {
                const int node = tri.node[i];
                ue[2 * i] = u[2 * node];
                ue[2 * i + 1] = u[2 * node + 1];
            }

            double Tavg = reference_temperature;
            if (T.size() == n_nodes) {
                Tavg = 0.0;
                for (int i = 0; i < 3; ++i) {
                    Tavg += T[tri.node[i]];
                }
                Tavg /= 3.0;
            }

            const vec3 strain = B * ue;
            const vec3 epsT = thermal_strain_vector(material, Tavg - reference_temperature);
            const vec3 sigma = C * (strain - epsT);

            for (int i = 0; i < 3; ++i) {
                const int node = tri.node[i];
                s.sigma_xx[node] += A * sigma[0];
                s.sigma_yy[node] += A * sigma[1];
                s.tau_xy[node] += A * sigma[2];
                w[node] += A;
            }
        }

        for (int i = 0; i < n_nodes; ++i) {
            if (w[i] > 0.0) {
                s.sigma_xx[i] /= w[i];
                s.sigma_yy[i] /= w[i];
                s.tau_xy[i] /= w[i];
            }

            const Node& p = mesh.nodes[i];
            vec2 er(p.x, p.y);
            if (er.norm() == 0.0) {
                er = vec2(1.0, 0.0);
            }
            else {
                er.normalize();
            }
            const vec2 et(-er.y(), er.x());

            Eigen::Matrix2d S;
            S << s.sigma_xx[i], s.tau_xy[i],
                s.tau_xy[i], s.sigma_yy[i];

            s.sigma_rr[i] = er.dot(S * er);
            s.sigma_tt[i] = et.dot(S * et);
            s.tau_rt[i] = er.dot(S * et);

            const double sx = s.sigma_xx[i];
            const double sy = s.sigma_yy[i];
            const double txy = s.tau_xy[i];
            s.von_mises[i] = std::sqrt((std::max)(0.0, sx * sx - sx * sy + sy * sy + 3.0 * txy * txy));
        }

        return s;
    }

    void write_results_txt(const std::string& filename,
        const Mesh& mesh,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& T,
        const NodalStressField& stress) {

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть файл результатов: " + filename);
        }

        out << std::setprecision(10);
        out << "id x y T ux uy ur sigma_xx sigma_yy tau_xy sigma_rr sigma_tt tau_rt von_mises\n";

        for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
            const Node& p = mesh.nodes[i];
            const int id = static_cast<int>(i);
            const vec2 ui = displacement_at_node(u, id);
            const double Ti = (T.size() == static_cast<int>(mesh.nodes.size()) ? T[id] : 0.0);

            out << id << ' '
                << p.x << ' ' << p.y << ' '
                << Ti << ' '
                << ui[0] << ' ' << ui[1] << ' '
                << radial_displacement_value(mesh, id, u) << ' '
                << stress.sigma_xx[id] << ' '
                << stress.sigma_yy[id] << ' '
                << stress.tau_xy[id] << ' '
                << stress.sigma_rr[id] << ' '
                << stress.sigma_tt[id] << ' '
                << stress.tau_rt[id] << ' '
                << stress.von_mises[id] << '\n';
        }
    }

    void write_tecplot_zone(std::ofstream& out,
        const std::string& zone_name,
        const Mesh& mesh,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& T,
        const NodalStressField* stress,
        double deformation_scale,
        double solution_time,
        bool use_deformed_coordinates,
        int strand_id = 0) {

        out << "ZONE T=\"" << zone_name << "\", N=" << mesh.nodes.size()
            << ", E=" << mesh.elements.size()
            << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE";

        // Для нестационарной анимации Tecplot должен знать,
        // какие зоны относятся к одному физическому телу во времени.
        // Поэтому для температурной истории явно задаются STRANDID и SOLUTIONTIME.
        if (strand_id > 0) {
            out << ", STRANDID=" << strand_id
                << ", SOLUTIONTIME=" << solution_time;
        }
        out << "\n";

        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const Node& p = mesh.nodes[i];
            const vec2 ui = displacement_at_node(u, i);
            const double x = use_deformed_coordinates ? p.x + deformation_scale * ui[0] : p.x;
            const double y = use_deformed_coordinates ? p.y + deformation_scale * ui[1] : p.y;
            const double Ti = (T.size() == static_cast<int>(mesh.nodes.size()) ? T[i] : 0.0);

            out << x << ' ' << y << ' '
                << p.x << ' ' << p.y << ' '
                << Ti << ' '
                << ui[0] << ' ' << ui[1] << ' '
                << radial_displacement_value(mesh, i, u);

            if (stress) {
                out << ' ' << stress->sigma_xx[i]
                    << ' ' << stress->sigma_yy[i]
                    << ' ' << stress->tau_xy[i]
                    << ' ' << stress->sigma_rr[i]
                    << ' ' << stress->sigma_tt[i]
                    << ' ' << stress->tau_rt[i]
                    << ' ' << stress->von_mises[i];
            }
            else {
                out << " 0 0 0 0 0 0 0";
            }
            out << '\n';
        }

        for (const auto& tri : mesh.elements) {
            out << tri.node[0] + 1 << ' '
                << tri.node[1] + 1 << ' '
                << tri.node[2] + 1 << '\n';
        }
    }

    void write_tecplot_final_fields(const std::string& filename,
        const BodyState& inner,
        const BodyState& outer,
        const CoupledResult& result,
        const NodalStressField& stress_inner,
        const NodalStressField& stress_outer) {

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть Tecplot-файл: " + filename);
        }

        out << std::setprecision(10);
        out << "TITLE=\"Thermoelastic contact results\"\n";
        out << "VARIABLES=\"X\" \"Y\" \"X0\" \"Y0\" \"T\" \"Ux\" \"Uy\" \"Ur\" "
            << "\"Sxx\" \"Syy\" \"Txy\" \"Srr\" \"Stt\" \"Srt\" \"VonMises\"\n";

        write_tecplot_zone(out, "inner_deformed", inner.mesh, result.inner.displacement,
            inner.thermal.temperature, &stress_inner, DEFORMATION_SCALE_FOR_TECPLOT, 0.0, true);
        write_tecplot_zone(out, "outer_deformed", outer.mesh, result.outer.displacement,
            outer.thermal.temperature, &stress_outer, DEFORMATION_SCALE_FOR_TECPLOT, 0.0, true);
    }

    const Eigen::VectorXd& temperature_at_time(const pipes2d::ThermalResult& thermal,
        double time,
        std::size_t& cached_index) {

        if (thermal.history.empty()) {
            return thermal.temperature;
        }

        // История температур обычно сохранена на каждом шаге времени.
        // Если одно тело достигло критерия раньше другого, для последующих кадров
        // оно должно оставаться видимым с последним найденным температурным полем.
        while (cached_index + 1 < thermal.history_time.size()
            && thermal.history_time[cached_index + 1] <= time + 1.0e-12) {
            ++cached_index;
        }
        return thermal.history[cached_index];
    }

    void append_temperature_times(std::vector<double>& times,
        const pipes2d::ThermalResult& thermal) {

        if (thermal.history_time.empty()) {
            times.push_back(thermal.final_time);
            return;
        }
        for (double t : thermal.history_time) {
            times.push_back(t);
        }
    }

    void write_tecplot_temperature_animation(const std::string& filename,
        const BodyState& inner,
        const BodyState& outer) {

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть Tecplot-файл: " + filename);
        }

        out << std::setprecision(10);
        out << "TITLE=\"Temperature animation\"\n";
        out << "VARIABLES=\"X\" \"Y\" \"X0\" \"Y0\" \"T\" \"Ux\" \"Uy\" \"Ur\" "
            << "\"Sxx\" \"Syy\" \"Txy\" \"Srr\" \"Stt\" \"Srt\" \"VonMises\"\n";

        Eigen::VectorXd zero_inner = Eigen::VectorXd::Zero(2 * static_cast<int>(inner.mesh.nodes.size()));
        Eigen::VectorXd zero_outer = Eigen::VectorXd::Zero(2 * static_cast<int>(outer.mesh.nodes.size()));

        std::vector<double> times;
        append_temperature_times(times, inner.thermal);
        append_temperature_times(times, outer.thermal);
        std::sort(times.begin(), times.end());
        times.erase(std::unique(times.begin(), times.end(),
            [](double a, double b) { return std::abs(a - b) < 1.0e-12; }),
            times.end());

        std::size_t inner_index = 0;
        std::size_t outer_index = 0;
        for (double time : times) {
            const Eigen::VectorXd& T_inner = temperature_at_time(inner.thermal, time, inner_index);
            const Eigen::VectorXd& T_outer = temperature_at_time(outer.thermal, time, outer_index);

            // Важно: в каждый момент времени записываются обе зоны.
            // STRANDID=1 соответствует внутренней трубе, STRANDID=2 — внешней.
            // Одинаковый SOLUTIONTIME заставляет Tecplot показывать обе трубы в одном кадре.
            write_tecplot_zone(out, "inner_T_t=" + std::to_string(time),
                inner.mesh, zero_inner, T_inner, nullptr, 1.0, time, false, 1);
            write_tecplot_zone(out, "outer_T_t=" + std::to_string(time),
                outer.mesh, zero_outer, T_outer, nullptr, 1.0, time, false, 2);
        }
    }

    void write_stress_profile_phi0(const std::string& filename,
        const BodyState& inner,
        const BodyState& outer,
        const CoupledResult& result,
        const NodalStressField& stress_inner,
        const NodalStressField& stress_outer) {

        struct Row {
            std::string body;
            double r = 0.0;
            double sigma_rr = 0.0;
            double sigma_tt = 0.0;
            double tau_rt = 0.0;
            double von_mises = 0.0;
            double ur = 0.0;
            double T = 0.0;
        };

        std::vector<Row> rows;

        auto add_body = [&](const std::string& name,
            const BodyState& body,
            const Eigen::VectorXd& u,
            const NodalStressField& s) {

                for (int i = 0; i < static_cast<int>(body.mesh.nodes.size()); ++i) {
                    const Node& p = body.mesh.nodes[i];
                    if (!is_on_phi0_ray(p)) {
                        continue;
                    }
                    Row row;
                    row.body = name;
                    row.r = std::hypot(p.x, p.y);
                    row.sigma_rr = s.sigma_rr[i];
                    row.sigma_tt = s.sigma_tt[i];
                    row.tau_rt = s.tau_rt[i];
                    row.von_mises = s.von_mises[i];
                    row.ur = radial_displacement_value(body.mesh, i, u);
                    row.T = body.thermal.temperature[i];
                    rows.push_back(row);
                }
            };

        add_body("inner", inner, result.inner.displacement, stress_inner);
        add_body("outer", outer, result.outer.displacement, stress_outer);

        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
            return a.r < b.r;
            });

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть файл профиля напряжений: " + filename);
        }

        out << std::setprecision(10);
        out << "body r sigma_rr sigma_tt tau_rt von_mises ur T\n";
        for (const Row& row : rows) {
            out << row.body << ' '
                << row.r << ' '
                << row.sigma_rr << ' '
                << row.sigma_tt << ' '
                << row.tau_rt << ' '
                << row.von_mises << ' '
                << row.ur << ' '
                << row.T << '\n';
        }
    }

    void write_contact_table(const std::string& filename, const CoupledResult& result) {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть файл контактного давления: " + filename);
        }

        out << std::setprecision(10);
        out << "id angle_rad length pressure_signed pressure_physical gap\n";

        const std::size_t n = (std::min)(result.contact_edges.size(), result.contact_pressure.size());
        for (std::size_t i = 0; i < n; ++i) {
            const ContactEdgeInfo& e = result.contact_edges[i];
            const double angle = std::atan2(e.normal.y(), e.normal.x());
            const double signed_p = result.contact_pressure[i];
            const double gap = (i < result.contact_gap.size() ? result.contact_gap[i] : 0.0);
            out << i << ' '
                << angle << ' '
                << e.length << ' '
                << signed_p << ' '
                << -signed_p << ' '
                << gap << '\n';
        }
    }

    void write_summary_table(const std::string& filename,
        double u0c1,
        double u0c2,
        double delta_h,
        double delta_an,
        const PressureStats& pressure,
        double p_formula_delta_h,
        double p_formula_delta_an,
        const ComparisonStats& cmp) {

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть таблицу результатов: " + filename);
        }

        out << std::setprecision(10);
        out << "parameter value unit\n";
        out << "u0c_inner " << u0c1 << " m\n";
        out << "u0c_outer " << u0c2 << " m\n";
        out << "Delta_h " << delta_h << " m\n";
        out << "Delta_an " << delta_an << " m\n";
        out << "p_contact_average " << pressure.average_physical << " Pa\n";
        out << "p_formula_Delta_h " << p_formula_delta_h << " Pa\n";
        out << "p_formula_Delta_an " << p_formula_delta_an << " Pa\n";
        out << "error_p_vs_Delta_h " << cmp.relative_to_delta_h_percent << " percent\n";
        out << "error_p_vs_Delta_an " << cmp.relative_to_delta_an_percent << " percent\n";
    }

    void write_report_txt(const std::string& filename,
        const std::string& grid_folder,
        double T0,
        double Tc,
        double dt,
        int max_thermal_steps,
        double thermal_tolerance,
        const BodyState& inner,
        const BodyState& outer,
        const char* mechanical_model_name,
        double u0c1,
        double u0c2,
        double delta_h,
        double delta_an,
        const PressureStats& pressure,
        double p_formula_delta_h,
        double p_formula_delta_an,
        const ComparisonStats& cmp,
        const GapStats& gap) {

        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Не удалось открыть текстовый отчет: " + filename);
        }

        const TemperatureStats Tin = temperature_stats(inner.thermal.temperature);
        const TemperatureStats Tout = temperature_stats(outer.thermal.temperature);
        const double delta_error_percent =
            100.0 * std::abs(delta_h - delta_an) / std::abs(delta_an);

        out << std::setprecision(10);
        out << "РЕЗУЛЬТАТЫ РАСЧЕТА\n";
        out << "=================\n\n";

        out << "Исходные данные\n";
        out << "---------------\n";
        out << "Каталог сетки: " << grid_folder << "\n";
        out << "T0, К: " << T0 << "\n";
        out << "Tc, К: " << Tc << "\n";
        out << "Механическая модель: " << mechanical_model_name << "\n";
        out << "Тепловая задача: нестационарная, до критерия останова\n";
        out << "dt, с: " << dt << "\n";
        out << "max_steps: " << max_thermal_steps << "\n";
        out << "Критерий max|T - Tc|, К: " << thermal_tolerance << "\n\n";

        out << "Контроль температуры\n";
        out << "--------------------\n";
        out << "Внутренняя труба: шагов = " << inner.thermal.performed_steps
            << ", max|T-Tc| = " << inner.thermal.final_max_error << " К"
            << ", Tmin/Tavg/Tmax = " << Tin.minimum << " / " << Tin.average << " / " << Tin.maximum << " К\n";
        out << "Внешняя труба: шагов = " << outer.thermal.performed_steps
            << ", max|T-Tc| = " << outer.thermal.final_max_error << " К"
            << ", Tmin/Tavg/Tmax = " << Tout.minimum << " / " << Tout.average << " / " << Tout.maximum << " К\n\n";

        out << "Свободное охлаждение\n";
        out << "--------------------\n";
        out << "u0c inner, м: " << u0c1 << "\n";
        out << "u0c outer, м: " << u0c2 << "\n";
        out << "Delta_h, м: " << delta_h << "\n";
        out << "Delta_an, м: " << delta_an << "\n";
        out << "Ошибка Delta_h относительно Delta_an, %: " << delta_error_percent << "\n\n";

        out << "Контактное давление\n";
        out << "-------------------\n";
        out << "Среднее контактное давление, Па: " << pressure.average_physical << "\n";
        out << "Давление по Delta_h, Па: " << p_formula_delta_h << "\n";
        out << "Давление по Delta_an, Па: " << p_formula_delta_an << "\n";
        out << "Расхождение с формулой по Delta_h, %: " << cmp.relative_to_delta_h_percent << "\n";
        out << "Расхождение с формулой по Delta_an, %: " << cmp.relative_to_delta_an_percent << "\n";
        out << "Средний остаточный зазор, м: " << gap.average << "\n";
        out << "Максимальный по модулю остаточный зазор, м: "
            << (std::max)(std::abs(gap.minimum), std::abs(gap.maximum)) << "\n\n";

        // Диагностический вывод. Раскомментировать при проверке контактного алгоритма.
        // out << "Диагностика контактного решения\n";
        // out << "-------------------------------\n";
        // out << "Средний множитель Лагранжа со знаком, Па: " << pressure.average_signed << "\n";
        // out << "Минимальный множитель Лагранжа со знаком, Па: " << pressure.minimum_signed << "\n";
        // out << "Максимальный множитель Лагранжа со знаком, Па: " << pressure.maximum_signed << "\n";
        // out << "Остаточный зазор min/avg/max, м: "
        //     << gap.minimum << " / " << gap.average << " / " << gap.maximum << "\n\n";

        out << "Файлы результатов\n";
        out << "-----------------\n";
        out << "result_report.txt\n";
        out << "summary_table.txt\n";
        out << "inner_results.txt\n";
        out << "outer_results.txt\n";
        out << "stress_profile_phi0.txt\n";
        out << "contact_pressure.txt\n";
        out << "tecplot_temperature_animation.dat\n";
        out << "tecplot_deformed_fields.dat\n";
    }

} // namespace

int main(int argc, char* argv[]) {
    try {
        const std::string grid_folder = select_grid_folder(argc, argv);
        const std::string inner_project_file = join_path(grid_folder, "project_inner.txt");
        const std::string outer_project_file = join_path(grid_folder, "project_outer.txt");

        if (!file_exists(inner_project_file)) {
            throw std::runtime_error("Не найден файл проекта внутренней трубы: " + inner_project_file);
        }
        if (!file_exists(outer_project_file)) {
            throw std::runtime_error("Не найден файл проекта внешней трубы: " + outer_project_file);
        }

        const double T0 = 400.0;
        const double Tc = 300.0;
        const double dt = 0.1;

        const int max_thermal_steps = 10000;
        const double thermal_tolerance = 1.0e-6;

        const double a = 1.0e-3;
        const double c = 4.0e-3;
        const double b = 5.0e-3;

        const int inner_contact_type = 4;
        const int outer_contact_type = 5;

        BodyState inner;
        inner.name = "inner";
        inner.material = Material{ 2.0e11, 0.33, 7.0e-6, 30.0, 1.0e4, 300.0, true };
        inner.mesh = read_gridder2d_project(inner_project_file);
        inner.reference_temperature = T0;

        BodyState outer;
        outer.name = "outer";
        outer.material = Material{ 2.0e11, 0.33, 7.0e-5, 30.0, 1.0e4, 300.0, true };
        outer.mesh = read_gridder2d_project(outer_project_file);
        outer.reference_temperature = T0;

        apply_quarter_symmetry_constraints(inner);
        apply_quarter_symmetry_constraints(outer);

        ThermalSolver2D thermal_inner(inner.mesh, inner.material);
        ThermalSolver2D thermal_outer(outer.mesh, outer.material);

        inner.thermal = thermal_inner.solve_unsteady_cooling_until_converged(
            inner_contact_type, Tc, dt, max_thermal_steps, T0, thermal_tolerance,
            THERMAL_HISTORY_STRIDE);

        outer.thermal = thermal_outer.solve_unsteady_cooling_until_converged(
            outer_contact_type, Tc, dt, max_thermal_steps, T0, thermal_tolerance,
            THERMAL_HISTORY_STRIDE);

        CoupledThermoMechanicalContact2D solver;

        const CoupledResult free_result = solver.solve(
            inner, outer,
            inner_contact_type, outer_contact_type,
            0.0);

        const double u0c1 = weighted_average_radial_displacement(
            inner.mesh, inner_contact_type, free_result.inner.displacement);
        const double u0c2 = weighted_average_radial_displacement(
            outer.mesh, outer_contact_type, free_result.outer.displacement);
        const double delta_h = u0c1 - u0c2;
        const double delta_an = c * (outer.material.alpha - inner.material.alpha) * (T0 - Tc);

        const CoupledResult result = solver.solve(
            inner, outer,
            inner_contact_type, outer_contact_type,
            1.0);

        const PressureStats pressure = pressure_stats_by_edges(result);
        const double p_formula_delta_h = explicit_contact_pressure_formula(delta_h, a, c, b, inner.material, outer.material);
        const double p_formula_delta_an = explicit_contact_pressure_formula(delta_an, a, c, b, inner.material, outer.material);
        const ComparisonStats cmp = make_comparison(pressure.average_physical, p_formula_delta_h, p_formula_delta_an);
        const GapStats gap = gap_stats(result.contact_gap);

        const NodalStressField stress_inner = compute_nodal_stress_field(
            inner.mesh, inner.material, result.inner.displacement,
            inner.thermal.temperature, inner.reference_temperature);
        const NodalStressField stress_outer = compute_nodal_stress_field(
            outer.mesh, outer.material, result.outer.displacement,
            outer.thermal.temperature, outer.reference_temperature);

        write_results_txt("inner_results.txt", inner.mesh, result.inner.displacement, inner.thermal.temperature, stress_inner);
        write_results_txt("outer_results.txt", outer.mesh, result.outer.displacement, outer.thermal.temperature, stress_outer);
        write_stress_profile_phi0("stress_profile_phi0.txt", inner, outer, result, stress_inner, stress_outer);
        write_contact_table("contact_pressure.txt", result);
        write_summary_table("summary_table.txt", u0c1, u0c2, delta_h, delta_an,
            pressure, p_formula_delta_h, p_formula_delta_an, cmp);
        write_tecplot_temperature_animation("tecplot_temperature_animation.dat", inner, outer);
        write_tecplot_final_fields("tecplot_deformed_fields.dat", inner, outer, result, stress_inner, stress_outer);

        write_report_txt("result_report.txt",
            grid_folder,
            T0, Tc,
            dt, max_thermal_steps, thermal_tolerance,
            inner, outer,
            plane_model_name(inner.material),
            u0c1, u0c2,
            delta_h, delta_an,
            pressure,
            p_formula_delta_h,
            p_formula_delta_an,
            cmp,
            gap);

        std::cout << std::setprecision(10);
        std::cout << "Расчет завершен\n";
        std::cout << "Сетка: " << grid_folder << "\n";
        std::cout << "Delta_h = " << delta_h << " м\n";
        std::cout << "p_avg = " << pressure.average_physical << " Па\n";
        std::cout << "Отчет: result_report.txt\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << '\n';
        return 1;
    }

    return 0;
}

