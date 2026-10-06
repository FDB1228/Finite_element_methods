#define _USE_MATH_DEFINES
#include "fem2d.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace fem2d;


static constexpr double E0 = 2e11;
static constexpr double nu0 = 0.3;
static constexpr bool plane_stress0 = false; // false=ПДС, true=ПНС

static inline void lame_params(double& lam, double& mu)
{
    mu = E0 / (2.0 * (1.0 + nu0));
    if (!plane_stress0)
    {
        // ПДС
        lam = E0 * nu0 / ((1.0 + nu0) * (1.0 - 2.0 * nu0));
    }
    else
    {
        // ПНС
        lam = E0 * nu0 / (1.0 - nu0 * nu0);
    }
}


// тесты 

// линейный тест:
static vec2 u_exact_0(double x, double y)
{
    vec2 u;
    u[0] = 1e-3 + 2e-3 * x - 1.5e-3 * y;
    u[1] = -5e-4 + 1e-3 * x + 3e-3 * y;
    return u;
}

static vec2 body_force_0(double, double)
{
    vec2 f;
    f.setZero();
    return f;
}


// нелинейный тесты

// синусы
static vec2 u_exact_1(double x, double y)
{
    const double s = std::sin(M_PI * x) * std::sin(M_PI * y);
    vec2 u;
    u[0] = s;
    u[1] = s;
    return u;
}

static vec2 body_force_1(double x, double y)
{
    double lam, mu;
    lame_params(lam, mu);

    const double s = std::sin(M_PI * x) * std::sin(M_PI * y);
    const double cxy = std::cos(M_PI * (x + y));

    const double common = M_PI * M_PI * (2.0 * mu * s - (lam + mu) * cxy);

    vec2 f;
    f[0] = common;
    f[1] = common;
    return f;
}

// параболы

static constexpr double A2 = 1e-3;

static vec2 u_exact_2(double x, double)
{
    vec2 u;
    u[0] = A2 * x * x;
    u[1] = 0.0;
    return u;
}

static vec2 body_force_2(double, double)
{
    double lam, mu;
    lame_params(lam, mu);

    vec2 f;
    f[0] = -2.0 * A2 * (lam + 2.0 * mu);
    f[1] = 0.0;
    return f;
}


// билиненый тест

static constexpr double A3 = 1e-3;

static vec2 u_exact_3(double x, double y)
{
    vec2 u;
    u[0] = A3 * x * y;
    u[1] = 0.0;
    return u;
}

static vec2 body_force_3(double, double)
{
    double lam, mu;
    lame_params(lam, mu);

    vec2 f;
    f[0] = 0.0;
    f[1] = -(lam + mu) * A3;
    return f;
}




int main()
{
    try
    {
        // читаем сетку из project-файла
        mesh_t mesh = read_gridder2d_project("Project3.txt");


        // решатель
        elasticity2d_solver_t solver;
        solver.set_mesh(&mesh);

        material_t mat;
        mat.E = 2e11;
        mat.nu = 0.3;
        mat.plane_stress = plane_stress0; // ПНС
        solver.set_material(mat);

        solver.set_body_force(&body_force_2);

        // условие Дирихле 
        std::vector<char> on_bnd(mesh.nodes.size(), 0);
        for (std::size_t i = 0; i < mesh.boundary.size(); ++i)
        {
            on_bnd[mesh.boundary[i].n1] = 1;
            on_bnd[mesh.boundary[i].n2] = 1;
        }

        for (std::size_t n = 0; n < on_bnd.size(); ++n)
        {
            if (on_bnd[n])
            {
                const double x = mesh.nodes[n][0];
                const double y = mesh.nodes[n][1];
                const vec2 ue = u_exact_2(x, y);

                // задаём значения узловых степеней свободы напрямую

                //для теста 0 2 3  поменять функцию выше
                solver.fix_ux(n, ue[0]);
                solver.fix_uy(n, ue[1]);
                
                //для теста 1
//                solver.fix_ux(n, 0.0);
//                solver.fix_uy(n, 0.0);
            }
        }

        // собираем систему и решаем
        solver.init();
        solver.assemble();
        Eigen::VectorXd u_h = solver.solve();

        // считаем ошибку 
        const fem2d::error_norms_t err = compute_error_norms(mesh, u_h, &u_exact_2);

        std::cout << std::setprecision(18);
        std::cout << "L2 error = " << err.l2 << "\n";
        std::cout << "C  error = " << err.c << "\n";


        // экспортируем значения в узлах
        export_nodal_results("results.txt", mesh, u_h, &u_exact_2);
        std::cout << "exported: results.txt\n";
    }
    catch (const std::exception& e)
    {
        std::cout << "error: " << e.what() << "\n";
    }
    return 0;
}
