"""Генератор сеток для программы Kyrsovaya_7sem (упругость, МКЭ).

Строит структурированную сетку в единичном квадрате [0, 1] x [0, 1]: N x N ячеек,
каждая делится диагональю на два прямоугольных треугольника (катет h = 1/N).
Сетки пишутся в текстовом формате, который читает fem2d.cpp (read_gridder2d_project), поэтому
программа работает с ними без изменений.

Создаются три файла и файл проекта:
    nodes<suffix>.dat    <число узлов NN> <число граничных узлов NBN>, затем строки "x y"
    elems<suffix>.dat    <число элементов NE>, затем строки "n1 n2 n3 e1 e2 e3 v"
    border<suffix>.dat   <число граничных рёбер NBR>, затем строки "n1 n2 E Type", затем 0
    Project<suffix>.txt  имена трёх файлов подряд

Нумерация узлов и элементов с 1. Поля e1 e2 e3 (соседние элементы) программа не использует
и заполняются нулями. Тип границы Type: 1 снизу, 2 справа, 3 сверху, 4 слева.
Граничные узлы нумеруются первыми (как в Gridder2D).

Использование:
    python tools/gen_mesh.py 10 20 40                  # h = 0.1, 0.05, 0.025
    python tools/gen_mesh.py 20 --out src/meshes       # в указанную папку
"""
import argparse
import os


def build(n):
    h = 1.0 / n
    # граничные узлы идут первыми
    coords = {}
    order = []
    for j in range(n + 1):
        for i in range(n + 1):
            if i in (0, n) or j in (0, n):
                coords[(i, j)] = None
                order.append((i, j))
    nbn = len(order)
    for j in range(1, n):
        for i in range(1, n):
            order.append((i, j))
    ident = {ij: k + 1 for k, ij in enumerate(order)}
    nodes = [(ij[0] * h, ij[1] * h) for ij in order]

    elems = []
    elem_of_cell = {}
    for j in range(n):
        for i in range(n):
            a, b = ident[(i, j)], ident[(i + 1, j)]
            c, d = ident[(i + 1, j + 1)], ident[(i, j + 1)]
            elems.append((a, b, c))
            elems.append((a, c, d))
            elem_of_cell[(i, j)] = (len(elems) - 1, len(elems))     # номера с 1

    border = []
    for i in range(n):                       # снизу, элемент (a, b, c)
        border.append((ident[(i, 0)], ident[(i + 1, 0)], elem_of_cell[(i, 0)][0], 1))
    for j in range(n):                       # справа, элемент (a, b, c)
        border.append((ident[(n, j)], ident[(n, j + 1)], elem_of_cell[(n - 1, j)][0], 2))
    for i in range(n - 1, -1, -1):           # сверху, элемент (a, c, d)
        border.append((ident[(i + 1, n)], ident[(i, n)], elem_of_cell[(i, n - 1)][1], 3))
    for j in range(n - 1, -1, -1):           # слева, элемент (a, c, d)
        border.append((ident[(0, j + 1)], ident[(0, j)], elem_of_cell[(0, j)][1], 4))
    return nodes, nbn, elems, border


def write_mesh(n, out_dir, suffix):
    nodes, nbn, elems, border = build(n)
    names = (f"nodes{suffix}.dat", f"elems{suffix}.dat", f"border{suffix}.dat")
    with open(os.path.join(out_dir, names[0]), "w") as f:
        f.write(f"{len(nodes)} {nbn}\n")
        for x, y in nodes:
            f.write(f"{x:.12g} {y:.12g}\n")
    with open(os.path.join(out_dir, names[1]), "w") as f:
        f.write(f"{len(elems)}\n")
        for a, b, c in elems:
            f.write(f"{a} {b} {c} 0 0 0 1\n")
    with open(os.path.join(out_dir, names[2]), "w") as f:
        f.write(f"{len(border)}\n")
        for a, b, e, t in border:
            f.write(f"{a} {b} {e} {t}\n")
        f.write("0\n")
    project = os.path.join(out_dir, f"Project{suffix}.txt")
    with open(project, "w") as f:
        f.write("\n".join(names) + "\n")
    return len(nodes), len(elems), project


def main():
    parser = argparse.ArgumentParser(description="Генератор сеток единичного квадрата")
    parser.add_argument("n", type=int, nargs="+", help="число ячеек по стороне, например 10 20 40")
    parser.add_argument("--out", default=".", help="папка для файлов (по умолчанию текущая)")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for n in args.n:
        nn, ne, project = write_mesh(n, args.out, f"_n{n}")
        print(f"{project}: h = {1.0 / n:g}, узлов {nn}, элементов {ne}")


if __name__ == "__main__":
    main()
