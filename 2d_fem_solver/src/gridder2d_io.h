#pragma once

#include "fem2d_common.h"

namespace pipes2d {

// Читает project-файл Gridder2D и по нему загружает сетку:
// узлы, элементы, граничные ребра и внутренние ограничения.
// Ожидается стандартный формат Gridder2D:
// 1) файл описания области
// 2) файл узлов
// 3) файл элементов
// 4) файл граничных ребер
// 5) шаг сетки
// 6) число циклов оптимизации
Mesh read_gridder2d_project(const std::string& project_file);

} // namespace pipes2d
