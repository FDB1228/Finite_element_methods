#pragma once
#include "fem2d_common.h"
#include "thermal_solver.h"

namespace pipes2d {

    struct FixedDof2D {
        int node = -1;       // номер узла в локальной сетке тела
        int component = 0;   // 0 -> ux, 1 -> uy
        double value = 0.0;
    };

    struct BodyState {
        Mesh mesh;
        Material material;
        ThermalResult thermal;
        BoundaryCondition2D bc;
        std::vector<FixedDof2D> fixed_dofs; // точечные закрепления, не зависящие от type границы
        double reference_temperature = 400.0;
        std::string name;
    };

    struct MechanicalBodyResult {
        Eigen::VectorXd displacement; // [ux0, uy0, ux1, uy1, ...]
    };

    struct ContactEdgeInfo {
        int slave_n1 = -1;
        int slave_n2 = -1;
        int master_n1 = -1;
        int master_n2 = -1;
        double master_t = 0.5;
        double length = 0.0;
        vec2 normal = vec2::Zero();
    };

    struct CoupledResult {
        MechanicalBodyResult inner;
        MechanicalBodyResult outer;

        // Старое поле оставлено для совместимости с прежним diplom.cpp.
        std::vector<ContactPair> pairs;

        // Новая контактная геометрия для segment-to-segment/mortar-постобработки.
        std::vector<ContactEdgeInfo> contact_edges;

        // Контактное давление на ребрах Γc. Индексы совпадают с contact_edges.
        std::vector<double> contact_pressure;

        // Остаточный нормальный зазор на контактных ребрах.
        std::vector<double> contact_gap;
    };

    class CoupledThermoMechanicalContact2D {
    public:
        // contact_mode_parameter <= 0: механическая задача без контакта.
        // contact_mode_parameter > 0: контактная задача с рёберными множителями Лагранжа.
        CoupledResult solve(const BodyState& inner,
            const BodyState& outer,
            int inner_contact_type,
            int outer_contact_type,
            double contact_mode_parameter) const;

        double radial_displacement(const Mesh& mesh,
            int node,
            const Eigen::VectorXd& u) const;

    private:
        std::vector<ContactPair> build_contact_pairs(const Mesh& inner,
            int inner_type,
            const Mesh& outer,
            int outer_type) const;

        std::vector<triplet> assemble_body(const BodyState& body,
            Eigen::VectorXd& rhs) const;
    };

} // namespace pipes2d
