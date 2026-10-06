"""Исследование сходимости программы elasticity_fem на структурированных сетках единичного квадрата.

Для каждого N строится сетка (tools/gen_mesh.py), запускается программа
    elasticity_fem Project_nN.txt <номер теста>
и из её вывода берутся нормы ошибки L2 и C. Результат печатается в виде таблицы
Markdown с отношениями ошибок на соседних сетках (для сетки, измельчённой вдвое,
отношение около 4 означает второй порядок).

Использование (из папки Elasticity_fem):
    python tools/run_convergence.py --exe src/build/elasticity_fem
    python tools/run_convergence.py --exe src/build/elasticity_fem --tests 1 --n 10 20 40 80 160 \\
        --plot reports/figures/convergence.png

Номера тестов: 0 линейная функция, 1 sin(pi x) sin(pi y), 2 парабола u1 = 1e-3 x^2,
3 билинейная u1 = 1e-3 x y (в отчёте это тесты 1-4).
"""
import argparse
import math
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_mesh  # noqa: E402

NAMES = {
    0: "линейная функция",
    1: "sin(πx)·sin(πy)",
    2: "парабола u₁ = 10⁻³·x²",
    3: "билинейная u₁ = 10⁻³·x·y",
}


def run_case(exe, mesh_dir, n, test):
    out = subprocess.run([os.path.abspath(exe), f"Project_n{n}.txt", str(test)],
                         cwd=mesh_dir, capture_output=True, text=True, check=True).stdout
    l2 = float(re.search(r"L2 error\s*=\s*(\S+)", out).group(1))
    c = float(re.search(r"C\s+error\s*=\s*(\S+)", out).group(1))
    return l2, c


def fmt(x):
    return "0" if x == 0 else f"{x:.3e}"


def ratio(prev, cur):
    if prev is None or cur == 0 or prev == 0:
        return ""
    return f"{prev / cur:.2f}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", required=True, help="путь к собранной программе elasticity_fem")
    p.add_argument("--n", type=int, nargs="+", default=[10, 20, 40, 80], help="число ячеек по стороне")
    p.add_argument("--tests", type=int, nargs="+", default=[0, 1, 2, 3], help="номера тестов 0..3")
    p.add_argument("--mesh-dir", default=os.path.join(os.path.dirname(__file__), "..", "meshes", "generated"),
                   help="куда писать сетки (по умолчанию meshes/generated)")
    p.add_argument("--plot", help="сохранить график сходимости в этот файл (нужен matplotlib)")
    args = p.parse_args()

    mesh_dir = os.path.abspath(args.mesh_dir)
    os.makedirs(mesh_dir, exist_ok=True)
    sizes = {}
    for n in args.n:
        nn, ne, _ = gen_mesh.write_mesh(n, mesh_dir, f"_n{n}")
        sizes[n] = (nn, ne)

    results = {}
    for t in args.tests:
        print(f"\nТест {t}: {NAMES[t]}\n")
        print("| h        | Узлов | Элементов | Ошибка L₂  | Отношение | Ошибка C   | Отношение |")
        print("| -------- | ----- | --------- | ---------- | --------- | ---------- | --------- |")
        prev = None
        results[t] = []
        for n in args.n:
            l2, c = run_case(args.exe, mesh_dir, n, t)
            results[t].append((n, l2, c))
            r2 = ratio(prev[0], l2) if prev else ""
            rc = ratio(prev[1], c) if prev else ""
            print(f"| {1.0 / n:<8g} | {sizes[n][0]:<5} | {sizes[n][1]:<9} | {fmt(l2):<10} | {r2:<9} | {fmt(c):<10} | {rc:<9} |")
            prev = (l2, c)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.2, 4.4))
        for t, rows in results.items():
            hs = [1.0 / n for n, l2, c in rows if l2 > 1e-14]
            if len(hs) < 2:
                continue
            ax.loglog(hs, [l2 for n, l2, c in rows if l2 > 1e-14], "o-", label=f"L₂, тест {t}: {NAMES[t]}")
            if t == 1:
                ax.loglog(hs, [c for n, l2, c in rows], "s--", label=f"C, тест {t}: {NAMES[t]}")
        h0 = 1.0 / args.n[0]
        ref = results[1][0][1] if 1 in results else 1.0
        ax.loglog([h0, 1.0 / args.n[-1]], [ref, ref * (1.0 / args.n[-1] / h0) ** 2], "k:", label="наклон 2 (h²)")
        ax.set_xlabel("h")
        ax.set_ylabel("ошибка")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nграфик: {args.plot}")


if __name__ == "__main__":
    main()
