#pragma once

#include <vector>
#include <array>
#include <string>
#include <cstddef>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace fem2d
{
    // узел сетки
    using node_t = std::array<double, 2>;
    
    // L_2 и С нормы
    struct error_norms_t
    {
        double l2 = 0.0;
        double c = 0.0;
    };

    // треугольный конечный элемент 
    struct tri_t
    {
        std::size_t n1 = 0;
        std::size_t n2 = 0;
        std::size_t n3 = 0;
    };

    // граничное ребро внешней границы
    struct boundary_edge_t
    {
        std::size_t n1 = 0;
        std::size_t n2 = 0;
        std::size_t elem = 0; // элемент, которому принадлежит ребро 
        int type = 0;         // тип ГУ
    };

    // ребро внутреннего ограничения 
    struct constraint_edge_t
    {
        std::size_t n1 = 0;
        std::size_t n2 = 0;
        std::size_t elem1 = 0;
        std::size_t elem2 = 0;
        int type = 0;
    };

    // структура сетки
    struct mesh_t
    {
        std::vector<node_t> nodes;                 // список узлов
        std::vector<tri_t>  elems;                 // список элементов
        std::vector<int>    elem_region;           // номер подобласти для каждого элемента 

        std::vector<boundary_edge_t>   boundary;   // внешняя граница 
        std::vector<constraint_edge_t> constraints;// внутренняя грпница

        std::size_t boundary_node_count = 0;       
    };

    // читаем project-файл
    mesh_t read_gridder2d_project(const std::string& project_file);


    // читает сетку из трёх файлов Gridder2D: nodes/elems/border
    mesh_t read_gridder2d_mesh(const std::string& nodes_file,
        const std::string& elems_file,
        const std::string& border_file);


    // параметры материала 
    struct material_t
    {
        double E = 2e11;           // модуль Юнга
        double nu = 0.3;           // коэффициент Пуассона
        bool plane_stress = false; // 0 = ПДС, true = ПНС
    };

    // создаем вектор перемещений
    using vec2 = Eigen::Matrix<double, 2, 1>;

    // точное решение
    using exact_u_func_t = vec2(*)(double, double);

    // ГУ 2-го рода
    struct traction_bc_t
    {
        int type = 0;                 // тип ГУ
        bool is_const = true;         // const / функциональное ГУ
        vec2 cval = vec2(0.0, 0.0);   // значение постоянного ГУ
        vec2(*func)(double, double) = nullptr; // функция ГУ
    };

    // ГУ Дирихле
    struct dirichlet_bc_t
    {
        std::size_t node = 0;
        int comp = 0; // компонента x, компонента y
        double value = 0.0;
    };

    // решатель 
    class elasticity2d_solver_t
    {
    public:
        // указатель на сетку
        void set_mesh(const mesh_t* mesh);

        // материал
        void set_material(const material_t& mat);

       // объемная сила
        void set_body_force(vec2(*f)(double, double));

        // ГУ 2 рода константное
        void add_traction_const(int type, const vec2& p);

        // ГУ 2 рода функциональное
        void add_traction_func(int type, vec2(*p)(double, double));

        // ГУ Дирихле u_x
        void fix_ux(std::size_t node, double ux);

        // ГУ Дирихле u_y
        void fix_uy(std::size_t node, double uy);

        // удобные функции: применить одинаковое значение на узлы,
        // которые встречаются на ребрах boundary с данным type
        void fix_ux_on_boundary_type(int type, double ux);
        void fix_uy_on_boundary_type(int type, double uy);

        // инициализация массивов
        void init();

        // сборка глобальной матрицы жёсткости и правой части 
        void assemble();

        // решаем систему с учётом ГУ Дирихле 
        Eigen::VectorXd solve() const;

    private:
        //  нулевой вектор 
        static vec2 zero_vec(double, double);

        //  матрица упругости
        static Eigen::Matrix3d make_C(const material_t& m);

        // вклад массовых сил
        void apply_tractions(Eigen::VectorXd& rhs) const;

        // строим систему kbc * u = rhs_bc с учётом Дирихле 
        void build_elimination_system(const Eigen::SparseMatrix<double>& Kraw,
            const Eigen::VectorXd& rhs_raw,
            Eigen::SparseMatrix<double>& Kbc,
            Eigen::VectorXd& rhs_bc) const;

    private:
        const mesh_t* mesh_ = nullptr; // сетка
        material_t mat_;               // материал
        std::size_t Nd_ = 0;           // число степеней свободы (2 * число узлов)

        vec2(*body_force_)(double, double) = &elasticity2d_solver_t::zero_vec;

        std::vector<traction_bc_t>  tractions_; // список ГУ 2 рода
        std::vector<dirichlet_bc_t> dirichlet_; // список ГУ Дирихле

        //  матрица жёсткости 
        std::vector<Eigen::Triplet<double>> Ktrip_;
        Eigen::VectorXd rhs_;
    };

    // вычисляем норму ошибки
    error_norms_t compute_error_norms(const mesh_t& mesh,
        const Eigen::VectorXd& u_h,
        exact_u_func_t u_exact);


    // экспорт результатов
    void export_nodal_results(const std::string& file_name,
        const mesh_t& mesh,
        const Eigen::VectorXd& u_h,
        exact_u_func_t u_exact = nullptr);
} 
