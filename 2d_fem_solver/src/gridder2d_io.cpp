#include "gridder2d_io.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace pipes2d {
    namespace fs = std::filesystem;

    // Коэффициент перевода координат сетки Gridder2D в метры.
    // Новые описатели области построены в миллиметровом масштабе:
    // 1, 4, 5 соответствуют 0.001, 0.004, 0.005 м.
    // Поэтому сразу после чтения узлов переводим координаты в СИ.
    static constexpr double GRIDDER_LENGTH_SCALE = 1.0e-3;

    // Выводит диапазон радиусов узлов сетки.
    // Это контроль того, что масштабирование координат выполнено корректно.
    static void print_radius_range(const Mesh& mesh, const std::string& name) {
        double r_min = 1.0e100;
        double r_max = -1.0e100;

        for (const auto& node : mesh.nodes) {
            const double r = std::sqrt(node.x * node.x + node.y * node.y);
            r_min = std::min(r_min, r);
            r_max = std::max(r_max, r);
        }

        std::cout << "\nПроверка радиусов сетки: " << name << "\n";
        std::cout << "r_min = " << r_min << "\n";
        std::cout << "r_max = " << r_max << "\n";
    }

    // Удаляет пробелы и служебные символы в начале и в конце строки.
    // Нужна для аккуратного чтения project-файла Gridder2D.
    static std::string trim(const std::string& s) {
        const auto a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        const auto b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

    // Читает все непустые строки project-файла.
    // В Gridder2D значимыми являются первые 6 строк, поэтому здесь
    // дополнительно проверяется, что их не меньше шести.
    static std::vector<std::string> read_project_lines(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Cannot open Gridder2D project file: " + filename);

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (!line.empty()) lines.push_back(line);
        }

        if (lines.size() < 6) {
            throw std::runtime_error(
                "Gridder2D project file must contain at least 6 non-empty lines");
        }
        return lines;
    }

    // Преобразует относительный путь из project-файла в абсолютный путь,
    // используя каталог, в котором лежит сам project-файл.
    static fs::path resolve_against(const fs::path& base_dir, const std::string& raw) {
        fs::path p(raw);
        if (p.is_relative()) p = base_dir / p;
        return p;
    }

    // Читает сетку по project-файлу Gridder2D.
    // Используются только файлы узлов, элементов и граничных ребер.
    // Геометрический файл, шаг сетки и число итераций сглаживания здесь
    // напрямую не нужны, потому что к этому моменту сетка уже построена Gridder2D.
    Mesh read_gridder2d_project(const std::string& project_file) {
        const auto lines = read_project_lines(project_file);
        const fs::path base = fs::absolute(fs::path(project_file)).parent_path();

        // В формате project-файла Gridder2D:
        // lines[0] — файл области,
        // lines[1] — файл узлов,
        // lines[2] — файл элементов,
        // lines[3] — файл граничных ребер.
        const fs::path nodes_path = resolve_against(base, lines[1]);
        const fs::path elems_path = resolve_against(base, lines[2]);
        const fs::path bnd_path = resolve_against(base, lines[3]);

        Mesh mesh;

        {
            // Чтение файла узлов.
            // Первая строка содержит:
            // nn  — число всех узлов,
            // nbn — число узлов на границе и ограничениях.
            // Поле nbn пока не используется, но считывается для корректного разбора формата.
            std::ifstream in(nodes_path);
            if (!in) throw std::runtime_error("Cannot open nodes file: " + nodes_path.string());

            int nn = 0;
            int nbn = 0;
            in >> nn >> nbn;
            if (!in || nn <= 0) throw std::runtime_error("Invalid nodes file header");

            mesh.nodes.resize(nn);
            for (int i = 0; i < nn; ++i) {
                in >> mesh.nodes[i].x >> mesh.nodes[i].y;

                // Gridder2D построил сетку по описателю в миллиметрах.
                // Вся расчетная часть программы работает в СИ, поэтому
                // координаты узлов сразу переводятся в метры.
                mesh.nodes[i].x *= GRIDDER_LENGTH_SCALE;
                mesh.nodes[i].y *= GRIDDER_LENGTH_SCALE;
            }

            if (!in) throw std::runtime_error("Failed while reading nodes file");
        }

        {
            // Чтение файла элементов.
            // Gridder2D хранит индексацию с 1, внутри программы переводим ее к 0-based.
            std::ifstream in(elems_path);
            if (!in) throw std::runtime_error("Cannot open elements file: " + elems_path.string());

            int ne = 0;
            in >> ne;
            if (!in || ne <= 0) throw std::runtime_error("Invalid elements file header");

            mesh.elements.resize(ne);
            for (int e = 0; e < ne; ++e) {
                Triangle tri;
                in >> tri.node[0] >> tri.node[1] >> tri.node[2]
                    >> tri.neigh[0] >> tri.neigh[1] >> tri.neigh[2] >> tri.subdomain;

                for (int i = 0; i < 3; ++i) {
                    tri.node[i] -= 1;
                }
                mesh.elements[e] = tri;
            }

            if (!in) throw std::runtime_error("Failed while reading elements file");
        }

        {
            // Чтение файла граничных ребер.
            // Первая секция — внешняя граница,
            // вторая секция — внутренние ограничения (constraints).
            std::ifstream in(bnd_path);
            if (!in) throw std::runtime_error("Cannot open boundary file: " + bnd_path.string());

            int nbr = 0;
            in >> nbr;
            if (!in || nbr < 0) throw std::runtime_error("Invalid boundary file header");

            mesh.boundary.resize(nbr);
            for (int i = 0; i < nbr; ++i) {
                auto& e = mesh.boundary[i];
                in >> e.n1 >> e.n2 >> e.owner >> e.type;
                e.n1 -= 1;
                e.n2 -= 1;
                e.owner -= 1;
            }

            int nrr = 0;
            in >> nrr;
            if (!in || nrr < 0) throw std::runtime_error("Invalid restriction section in boundary file");

            mesh.constraints.resize(nrr);
            for (int i = 0; i < nrr; ++i) {
                auto& e = mesh.constraints[i];
                int second_owner = -1;
                in >> e.n1 >> e.n2 >> e.owner >> second_owner >> e.type;
                e.n1 -= 1;
                e.n2 -= 1;
                e.owner -= 1;

                // В текущей структуре BoundaryEdge хранится только один владеющий элемент.
                // Поэтому второй соседний элемент для внутреннего ограничения пока
                // считывается, но не сохраняется.
            }

            if (!in) throw std::runtime_error("Failed while reading boundary file");
        }

        // Диагностический вывод масштаба сетки.
        // Раскомментировать при проблемах с единицами измерения.
        // print_radius_range(mesh, project_file);

        return mesh;
    }

} // namespace pipes2d
