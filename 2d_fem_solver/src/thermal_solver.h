#pragma once
#include "fem2d_common.h"

namespace pipes2d {

    /// Результат теплового расчета.
    struct ThermalResult {
        /// Конечное температурное поле.
        Eigen::VectorXd temperature;

        /// История температур для Tecplot-анимации.
        /// Заполняется только если в решателе задан положительный history_stride.
        std::vector<Eigen::VectorXd> history;
        std::vector<double> history_time;

        /// Служебные параметры нестационарного расчета.
        int performed_steps = 0;
        double final_time = 0.0;
        double final_max_error = 0.0;
    };

    /// Решатель двумерной задачи теплопроводности.
    class ThermalSolver2D {
    public:
        ThermalSolver2D(const Mesh& mesh, const Material& material);

        /// Стационарная задача теплопроводности.
        ThermalResult solve_steady(
            const BoundaryCondition2D& bc,
            const std::function<double(double, double)>& source = nullptr) const;

        /// Нестационарная задача, неявная схема Эйлера с диагональной матрицей теплоемкости.
        ThermalResult solve_transient_implicit(
            const BoundaryCondition2D& bc,
            double initial_temperature,
            double dt,
            int steps = 1,
            const std::function<double(double, double, double)>& source = nullptr) const;

        /// Охлаждение за заданное число шагов.
        /// На границе cold_boundary_type задается T = cold_temperature.
        ThermalResult solve_unsteady_cooling(
            int cold_boundary_type,
            double cold_temperature,
            double dt,
            int steps,
            double initial_temperature) const;

        /// Охлаждение до критерия max_i |T_i - cold_temperature| <= temperature_tolerance.
        /// history_stride > 0 включает сохранение истории температур для Tecplot.
        ThermalResult solve_unsteady_cooling_until_converged(
            int cold_boundary_type,
            double cold_temperature,
            double dt,
            int max_steps,
            double initial_temperature,
            double temperature_tolerance,
            int history_stride = 0) const;

    private:
        const Mesh& mesh_;
        const Material& material_;
    };

} // namespace pipes2d
