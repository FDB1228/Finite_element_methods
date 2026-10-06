#include "fem2d.h"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <locale>
#include <cstdlib>

namespace fem2d
{
    // читаем project-файл:
    // файл координат узлов
    // файл элементов
    // файл граничных ребер
    mesh_t read_gridder2d_project(const std::string& project_file)
    {
        std::ifstream in(project_file.c_str());
        if (!in)
            throw std::logic_error("cannot open project file: " + project_file);

        std::string nodes_file;
        std::string elems_file;
        std::string border_file;

        in >> nodes_file >> elems_file >> border_file;

        if (!in)
            throw std::logic_error("bad project file format (need 3 names): " + project_file);

        return read_gridder2d_mesh(nodes_file, elems_file, border_file);
    }


  
    // Импортируем узлы
    // общее количество узлов NN 
    // узлы на границе NBN  
    // NN строк с координатами (x y)
    static void read_gridder2d_nodes(const std::string& file, mesh_t& mesh)
    {
        std::ifstream in(file.c_str());
        if (!in) throw std::logic_error("cannot open nodes file: " + file);

        long long NN = 0, NBN = 0;
        in >> NN >> NBN;

        mesh.nodes.resize((std::size_t)NN);
        mesh.boundary_node_count = (std::size_t)NBN;

        for (long long i = 0; i < NN; ++i)
        {
            in >> mesh.nodes[(std::size_t)i][0] >> mesh.nodes[(std::size_t)i][1];
        }
    }


    // переводим индекс с 1 -> 0 
    static std::size_t to_index_0based(long long v)
    {
        return (std::size_t)(v - 1);
    }


    // импортируем элементы
    // количество элементов NE 
    // NE строк n1 n2 n3 e1 e2 e3 v
    static void read_gridder2d_elems(const std::string& file, mesh_t& mesh)
    {
        std::ifstream in(file.c_str());
        if (!in) throw std::logic_error("cannot open elems file: " + file);

        long long NE = 0;
        in >> NE;
        if (!in) throw std::logic_error("bad elems file header: " + file);

        mesh.elems.resize((std::size_t)NE);
        mesh.elem_region.resize((std::size_t)NE);

        for (long long e = 0; e < NE; ++e)
        {
            long long n1 = 0, n2 = 0, n3 = 0;
            long long e1 = 0, e2 = 0, e3 = 0;
            long long v = 0;

            in >> n1 >> n2 >> n3 >> e1 >> e2 >> e3 >> v;

            if (!in) throw std::logic_error("bad elems file data: " + file);

            tri_t t;
            t.n1 = to_index_0based(n1);
            t.n2 = to_index_0based(n2);
            t.n3 = to_index_0based(n3);

            mesh.elems[(std::size_t)e] = t;
            mesh.elem_region[(std::size_t)e] = (int)v;
        }
    }


    // импортируем граничные ребра
    // количество граничных ребер NBR
    // NBR строк n1 n2 E Type
    // строка NBR+2:
    // количество ребер на ограничениях NRR 
    // далее NRR строк n1 n2 E1 E2 Type
    static void read_gridder2d_border(const std::string& file, mesh_t& mesh)
    {
        std::ifstream in(file.c_str());
        if (!in) throw std::logic_error("cannot open border file: " + file);

        long long NBR = 0;
        in >> NBR;
        if (!in) throw std::logic_error("bad border file header: " + file);

        mesh.boundary.resize((std::size_t)NBR);

        for (long long i = 0; i < NBR; ++i)
        {
            long long n1 = 0, n2 = 0, E = 0, Type = 0;

            in >> n1 >> n2 >> E >> Type;

            if (!in) throw std::logic_error("bad border file data (boundary): " + file);

            boundary_edge_t be;
            be.n1 = to_index_0based(n1);
            be.n2 = to_index_0based(n2);
            be.elem = to_index_0based(E);
            be.type = (int)Type;

            mesh.boundary[(std::size_t)i] = be;
        }


        long long NRR = 0;
        if (!(in >> NRR))
        {
            mesh.constraints.clear();
            return;
        }

        mesh.constraints.resize((std::size_t)NRR);

        for (long long i = 0; i < NRR; ++i)
        {
            long long n1 = 0, n2 = 0, E1 = 0, E2 = 0, Type = 0;

            in >> n1 >> n2 >> E1 >> E2 >> Type;

            if (!in) throw std::logic_error("bad border file data (constraints): " + file);

            constraint_edge_t ce;
            ce.n1 = to_index_0based(n1);
            ce.n2 = to_index_0based(n2);
            ce.elem1 = to_index_0based(E1);
            ce.elem2 = to_index_0based(E2);
            ce.type = (int)Type;

            mesh.constraints[(std::size_t)i] = ce;
        }
    }

    // создаём и заполняем сетку
    mesh_t read_gridder2d_mesh(const std::string& nodes_file,
        const std::string& elems_file,
        const std::string& border_file)
    {
        mesh_t mesh;
        read_gridder2d_nodes(nodes_file, mesh);
        read_gridder2d_elems(elems_file, mesh);
        read_gridder2d_border(border_file, mesh);

        return mesh;
    }

   // решатель

    vec2 elasticity2d_solver_t::zero_vec(double, double)
    {
        vec2 z;
        z.setZero();
        return z;
    }

    void elasticity2d_solver_t::set_mesh(const mesh_t* mesh)
    {
        mesh_ = mesh;
    }

    void elasticity2d_solver_t::set_material(const material_t& mat)
    {
        mat_ = mat;
    }

    void elasticity2d_solver_t::set_body_force(vec2(*f)(double, double))
    {
        body_force_ = (f ? f : &zero_vec);
    }

    void elasticity2d_solver_t::add_traction_const(int type, const vec2& p)
    {
        traction_bc_t bc;
        bc.type = type;
        bc.is_const = true;
        bc.cval = p;
        bc.func = nullptr;
        tractions_.push_back(bc);
    }

    void elasticity2d_solver_t::add_traction_func(int type, vec2(*p)(double, double))
    {
        traction_bc_t bc;
        bc.type = type;
        bc.is_const = false;
        bc.cval.setZero();
        bc.func = (p ? p : &zero_vec);
        tractions_.push_back(bc);
    }

    void elasticity2d_solver_t::fix_ux(std::size_t node, double ux)
    {
        dirichlet_bc_t bc;
        bc.node = node;
        bc.comp = 0;
        bc.value = ux;
        dirichlet_.push_back(bc);
    }

    void elasticity2d_solver_t::fix_uy(std::size_t node, double uy)
    {
        dirichlet_bc_t bc;
        bc.node = node;
        bc.comp = 1;
        bc.value = uy;
        dirichlet_.push_back(bc);
    }

    void elasticity2d_solver_t::fix_ux_on_boundary_type(int type, double ux)
    {
        std::vector<char> mark(mesh_->nodes.size(), 0);
        for (std::size_t i = 0; i < mesh_->boundary.size(); ++i)
        {
            const boundary_edge_t& e = mesh_->boundary[i];
            if (e.type == type)
            {
                mark[e.n1] = 1;
                mark[e.n2] = 1;
            }
        }

        for (std::size_t n = 0; n < mark.size(); ++n)
        {
            if (mark[n]) fix_ux(n, ux);
        }
    }

    void elasticity2d_solver_t::fix_uy_on_boundary_type(int type, double uy)
    {
        std::vector<char> mark(mesh_->nodes.size(), 0);
        for (std::size_t i = 0; i < mesh_->boundary.size(); ++i)
        {
            const boundary_edge_t& e = mesh_->boundary[i];
            if (e.type == type)
            {
                mark[e.n1] = 1;
                mark[e.n2] = 1;
            }
        }

        for (std::size_t n = 0; n < mark.size(); ++n)
        {
            if (mark[n]) fix_uy(n, uy);
        }
    }

    Eigen::Matrix3d elasticity2d_solver_t::make_C(const material_t& m)
    {
        // матрица упругости 
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();

        if (!m.plane_stress)
        {
            // ПДС
            const double mu = m.E / (2.0 * (1.0 + m.nu));
            const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));

            C(0, 0) = lam + 2.0 * mu;  C(0, 1) = lam;
            C(1, 0) = lam;             C(1, 1) = lam + 2.0 * mu;
            C(2, 2) = mu;
        }
        else
        {
            // ПНС
            const double k = m.E / (1.0 - m.nu * m.nu);

            C(0, 0) = k;        C(0, 1) = k * m.nu;
            C(1, 0) = k * m.nu; C(1, 1) = k;
            C(2, 2) = k * (1.0 - m.nu) / 2.0;
        }

        return C;
    }

    void elasticity2d_solver_t::init()
    {
        // число степеней свободы: u_x и u_y в каждом узле
        Nd_ = 2 * mesh_->nodes.size();

        // правая часть 
        rhs_ = Eigen::VectorXd::Zero((Eigen::Index)Nd_);

      
        // на треугольнике 6 степеней свободы, локальная матрица 6x6 
        Ktrip_.clear();
        Ktrip_.reserve(mesh_->elems.size() * 36 + 64);
    }

    void elasticity2d_solver_t::assemble()
    {

        rhs_.setZero();
        Ktrip_.clear();
        Ktrip_.reserve(mesh_->elems.size() * 36 + 64);

        const Eigen::Matrix3d C = make_C(mat_);
        const std::vector<node_t>& nodes = mesh_->nodes;

        // сборка по элементам треугольника
        // для каждого треугольника строим локальную матрицу ke и добавляем в Ktrip_
        for (std::size_t e = 0; e < mesh_->elems.size(); ++e)
        {
            const tri_t& t = mesh_->elems[e];
            const std::size_t I = t.n1;
            const std::size_t J = t.n2;
            const std::size_t K = t.n3;

            const double x1 = nodes[I][0], y1 = nodes[I][1];
            const double x2 = nodes[J][0], y2 = nodes[J][1];
            const double x3 = nodes[K][0], y3 = nodes[K][1];

            const double A2 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            const double absA2 = std::abs(A2);
            if (absA2 == 0.0) continue;

            const double area = 0.5 * absA2;

            // коэффициенты для градиентов базисных функций 
            const double b1 = y2 - y3;
            const double b2 = y3 - y1;
            const double b3 = y1 - y2;

            const double c1 = x3 - x2;
            const double c2 = x1 - x3;
            const double c3 = x2 - x1;

            // производные базисных функций
            const double dphix1 = b1 / A2;
            const double dphiy1 = c1 / A2;
            const double dphix2 = b2 / A2;
            const double dphiy2 = c2 / A2;
            const double dphix3 = b3 / A2;
            const double dphiy3 = c3 / A2;

            Eigen::Matrix<double, 3, 6> R;
            R <<
                dphix1, 0.0, dphix2, 0.0, dphix3, 0.0,
                0.0, dphiy1, 0.0, dphiy2, 0.0, dphiy3,
                dphiy1, dphix1, dphiy2, dphix2, dphiy3, dphix3;

            // локальная матрица жёсткости
            Eigen::Matrix<double, 6, 6> Ke = area * (R.transpose() * (C * R));

            // глобальные номера степеней свободы для трёх узлов
            const std::size_t dof[6] =
            {
                2 * I, 2 * I + 1,
                2 * J, 2 * J + 1,
                2 * K, 2 * K + 1
            };

            for (int a = 0; a < 6; ++a)
            {
                for (int b = 0; b < 6; ++b)
                {
                    Ktrip_.push_back(Eigen::Triplet<double>((int)dof[a], (int)dof[b], Ke(a, b)));
                }
            }

            // вклад объёмной силы в rhs
            const double xm = (x1 + x2 + x3) / 3.0;
            const double ym = (y1 + y2 + y3) / 3.0;
            const vec2 f = body_force_(xm, ym);

            const double w = absA2 / 6.0;

            rhs_[(Eigen::Index)dof[0]] += f[0] * w;
            rhs_[(Eigen::Index)dof[1]] += f[1] * w;
            rhs_[(Eigen::Index)dof[2]] += f[0] * w;
            rhs_[(Eigen::Index)dof[3]] += f[1] * w;
            rhs_[(Eigen::Index)dof[4]] += f[0] * w;
            rhs_[(Eigen::Index)dof[5]] += f[1] * w;
        }

        // добавляем в rhs ГУ 2 рода 
        apply_tractions(rhs_);
    }

    void elasticity2d_solver_t::apply_tractions(Eigen::VectorXd& rhs) const
    {
        if (tractions_.empty()) return;

        const std::vector<node_t>& nodes = mesh_->nodes;

        for (std::size_t i = 0; i < mesh_->boundary.size(); ++i)
        {
            const boundary_edge_t& e = mesh_->boundary[i];

            const traction_bc_t* bc = nullptr;
            for (std::size_t k = 0; k < tractions_.size(); ++k)
            {
                if (tractions_[k].type == e.type)
                {
                    bc = &tractions_[k];
                    break;
                }
            }
            if (!bc) continue;

            const std::size_t I = e.n1;
            const std::size_t J = e.n2;

            const double x1 = nodes[I][0], y1 = nodes[I][1];
            const double x2 = nodes[J][0], y2 = nodes[J][1];

            const double xm = 0.5 * (x1 + x2);
            const double ym = 0.5 * (y1 + y2);

            const double dx = x2 - x1;
            const double dy = y2 - y1;
            const double len = std::sqrt(dx * dx + dy * dy);
            if (len == 0.0) continue;

            vec2 p;
            if (bc->is_const) p = bc->cval;
            else p = (bc->func ? bc->func(xm, ym) : zero_vec(xm, ym));
            const double w = 0.5 * len;

            rhs[(Eigen::Index)(2 * I)] += p[0] * w;
            rhs[(Eigen::Index)(2 * I + 1)] += p[1] * w;
            rhs[(Eigen::Index)(2 * J)] += p[0] * w;
            rhs[(Eigen::Index)(2 * J + 1)] += p[1] * w;
        }
    }

    void elasticity2d_solver_t::build_elimination_system(const Eigen::SparseMatrix<double>& Kraw,
        const Eigen::VectorXd& rhs_raw,
        Eigen::SparseMatrix<double>& Kbc,
        Eigen::VectorXd& rhs_bc) const
    {
        // учитываем ГУ Дирихле

        std::vector<char> fixed(Nd_, 0);
        std::vector<double> u0(Nd_, 0.0);


        for (std::size_t k = 0; k < dirichlet_.size(); ++k)
        {
            const dirichlet_bc_t& bc = dirichlet_[k];
            const std::size_t dof = 2 * bc.node + (std::size_t)bc.comp;
            fixed[dof] = 1;
            u0[dof] = bc.value;
        }

        rhs_bc = rhs_raw;

        std::vector<Eigen::Triplet<double>> t;
        t.reserve((std::size_t)Kraw.nonZeros() + dirichlet_.size() + 16);

        // пробегаем все ненулевые элементы матрицы
        for (int col = 0; col < Kraw.outerSize(); ++col)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Kraw, col); it; ++it)
            {
                const int i = it.row();
                const int j = it.col();
                const double aij = it.value();

                const bool fi = fixed[(std::size_t)i] != 0;
                const bool fj = fixed[(std::size_t)j] != 0;

                if (!fi && !fj)
                {
                    t.push_back(Eigen::Triplet<double>(i, j, aij));
                }
                else if (!fi && fj)
                {
                    rhs_bc[i] -= aij * u0[(std::size_t)j];
                }
            }
        }

        for (std::size_t d = 0; d < Nd_; ++d)
        {
            if (fixed[d])
            {
                t.push_back(Eigen::Triplet<double>((int)d, (int)d, 1.0));
                rhs_bc[(Eigen::Index)d] = u0[d];
            }
        }

        Kbc.resize((Eigen::Index)Nd_, (Eigen::Index)Nd_);
        Kbc.setFromTriplets(t.begin(), t.end());
        Kbc.makeCompressed();
    }

    Eigen::VectorXd elasticity2d_solver_t::solve() const
    {
        // собираем разреженную матрицу Kraw
        Eigen::SparseMatrix<double> Kraw((Eigen::Index)Nd_, (Eigen::Index)Nd_);
        Kraw.setFromTriplets(Ktrip_.begin(), Ktrip_.end());
        Kraw.makeCompressed();

        // применяем ГУ Дирихле
        Eigen::SparseMatrix<double> Kbc;
        Eigen::VectorXd rhs_bc;
        build_elimination_system(Kraw, rhs_, Kbc, rhs_bc);

        // решаем систему
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(Kbc);
        if (ldlt.info() != Eigen::Success)
            throw std::logic_error("ldlt factorization failed (not enough constraints or bad matrix)");

        Eigen::VectorXd u = ldlt.solve(rhs_bc);
        if (ldlt.info() != Eigen::Success)
            throw std::logic_error("ldlt solve failed");

        return u;
    }

    //вычисляем норму ошибку
    error_norms_t compute_error_norms(const mesh_t& mesh,
        const Eigen::VectorXd& u_h,
        exact_u_func_t u_exact)
    {
        if (!u_exact)
            throw std::logic_error("compute_error_norms: u_exact == nullptr");

        const std::size_t NN = mesh.nodes.size();
        if ((std::size_t)u_h.size() != 2 * NN)
            throw std::logic_error("compute_error_norms: размер u_h должен быть 2*NN");

        error_norms_t out;

        // норма C
        double cmax = 0.0;
        for (std::size_t i = 0; i < NN; ++i)
        {
            const double x = mesh.nodes[i][0];
            const double y = mesh.nodes[i][1];

            const Eigen::Matrix<double, 2, 1> ue = u_exact(x, y);

            const double ux_h = u_h[(Eigen::Index)(2 * i)];
            const double uy_h = u_h[(Eigen::Index)(2 * i + 1)];

            const double dux = ux_h - ue[0];
            const double duy = uy_h - ue[1];

            const double e = std::sqrt(dux * dux + duy * duy);
            if (e > cmax) cmax = e;
        }
        out.c = cmax;

        // норма L2:

        double sum = 0.0;

        for (std::size_t e = 0; e < mesh.elems.size(); ++e)
        {
            const auto& t = mesh.elems[e];

            const std::size_t I = t.n1;
            const std::size_t J = t.n2;
            const std::size_t K = t.n3;

            const double x1 = mesh.nodes[I][0], y1 = mesh.nodes[I][1];
            const double x2 = mesh.nodes[J][0], y2 = mesh.nodes[J][1];
            const double x3 = mesh.nodes[K][0], y3 = mesh.nodes[K][1];


            const double A2 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            const double area = 0.5 * std::abs(A2);
            if (area == 0.0) continue;


            const double xc = (x1 + x2 + x3) / 3.0;
            const double yc = (y1 + y2 + y3) / 3.0;


            const double uxI = u_h[(Eigen::Index)(2 * I)];
            const double uyI = u_h[(Eigen::Index)(2 * I + 1)];
            const double uxJ = u_h[(Eigen::Index)(2 * J)];
            const double uyJ = u_h[(Eigen::Index)(2 * J + 1)];
            const double uxK = u_h[(Eigen::Index)(2 * K)];
            const double uyK = u_h[(Eigen::Index)(2 * K + 1)];

            const double ux_c = (uxI + uxJ + uxK) / 3.0;
            const double uy_c = (uyI + uyJ + uyK) / 3.0;

            const Eigen::Matrix<double, 2, 1> ue = u_exact(xc, yc);

            const double dux = ux_c - ue[0];
            const double duy = uy_c - ue[1];

            const double e2 = dux * dux + duy * duy;

            sum += area * e2;
        }

        out.l2 = std::sqrt(sum);
        return out;
    }

    // экспорт результатов
    void export_nodal_results(const std::string& file_name,
        const mesh_t& mesh,
        const Eigen::VectorXd& u_h,
        exact_u_func_t u_exact)
    {
        std::ofstream out(file_name.c_str());
        if (!out) throw std::logic_error("unable to write file: " + file_name);
        out.imbue(std::locale::classic());

        out << std::setprecision(16);

        const std::size_t NN = mesh.nodes.size();

        out << NN << "\n";
        if (u_exact)
            out << "# id x y ux uy ux_exact uy_exact err\n";
        else
            out << "# id x y ux uy\n";

        for (std::size_t i = 0; i < NN; ++i)
        {
            const double x = mesh.nodes[i][0];
            const double y = mesh.nodes[i][1];

            const double ux_h = u_h[(Eigen::Index)(2 * i)];
            const double uy_h = u_h[(Eigen::Index)(2 * i + 1)];

            out << i << " " << x << " " << y << " " << ux_h << " " << uy_h;

            if (u_exact)
            {
                const vec2 ue = u_exact(x, y);

                const double dux = ux_h - ue[0];
                const double duy = uy_h - ue[1];
                const double err = std::sqrt(dux * dux + duy * duy);

                out << " " << ue[0] << " " << ue[1] << " " << err;
            }

            out << "\n";
        }
    }
} 
