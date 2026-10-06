"""График решения по файлу results.txt (замена Tecplot / ParaView для быстрой проверки).

Строит две карты на сетке из файла проекта: перемещение u_x (численное) и ошибку
|u_h - u| в узлах. Программа elasticity_fem записывает results.txt в рабочую папку, то есть
в папку с файлами сетки.

Использование:
    cd meshes/generated
    ../../src/build/elasticity_fem Project_n20.txt 1
    python ../../tools/plot_solution.py Project_n20.txt results.txt ../../reports/figures/solution_sin.png

Нужны numpy и matplotlib.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def read_project(project):
    base = os.path.dirname(os.path.abspath(project))
    with open(project, encoding="utf-8") as f:
        nodes_f, elems_f, border_f = f.read().split()[:3]
    with open(os.path.join(base, nodes_f)) as f:
        nn = int(f.readline().split()[0])
        xy = np.loadtxt(f, max_rows=nn)
    with open(os.path.join(base, elems_f)) as f:
        ne = int(f.readline().split()[0])
        el = np.loadtxt(f, max_rows=ne, dtype=int)[:, :3] - 1       # нумерация с 1
    return xy, el


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    project, results, out = sys.argv[1:]
    xy, el = read_project(project)
    data = np.loadtxt(results, skiprows=2)
    ux, uy, ux_e, uy_e, err = data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1], el)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    m = axes[0].tripcolor(tri, ux, shading="gouraud", cmap="viridis")
    axes[0].triplot(tri, color="k", lw=0.2, alpha=0.4)
    axes[0].set_title("Перемещение $u_x$ (МКЭ)")
    fig.colorbar(m, ax=axes[0])
    m = axes[1].tripcolor(tri, err, shading="gouraud", cmap="magma")
    axes[1].set_title(r"Ошибка в узлах $|u_h-u|$")
    fig.colorbar(m, ax=axes[1], format="%.1e")
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"сохранено: {out}; max|u_h - u| = {err.max():.3e}")


if __name__ == "__main__":
    main()
