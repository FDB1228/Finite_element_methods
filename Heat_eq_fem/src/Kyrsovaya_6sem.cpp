#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <locale>
#include <map>
#include <algorithm>
#include <numeric>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include <Eigen\Dense>


//РАСКОММЕНТИРОВАТЬ ДЛЯ ВЕРСИИ EIGEN МЕНЬШЕ 3.4
//double* begin(Eigen::VectorXd& vec)
// {
// return vec.data();
// }
//double* end(Eigen::VectorXd& vec)
// {
// return vec.data() + vec.size();
// }

template<typename S, typename VV, typename N = std::initializer_list<const char*>>
void write_vec_vec(S& out, const VV& vv, const N& names)
{
	if (vv.empty()) return;

	assert(vv[0].size() >= (int)names.size());

	const auto endl = '\n';

	size_t i = 0;
	for (auto& name : names)
	{
		out << name << endl;
		for (auto& v : vv) out << v[i] << endl;
		++i;
	}
}

inline const auto endl = '\n';

void VAR(double& var, const char* name);
void read_vars(std::istream& in);
void write_vars(std::ostream& out);

void read_input(const char* file_name = "input.txt")
{
	std::ifstream in(file_name);
	if (!in) throw std::logic_error("input file not found");
	read_vars(in);
}


using node_t = std::array<double, 2>;
using elem_t = std::array<size_t, 3>;
using bound_t = std::vector<size_t>; //номера узлов вдоль границы

struct mesh_t
{
	std::vector<node_t> nodes;
	std::vector<elem_t> elems;
	std::vector<bound_t> bounds;
};

const double SMALL_JAC = 1e-8;

mesh_t read_mesh(std::istream& in)
{
	mesh_t mesh;

	//READING NODES
	size_t num_nodes = 0;
	in >> num_nodes;

	auto& nodes = mesh.nodes;

	nodes.resize(num_nodes);

	for (size_t i = 0; i < num_nodes; ++i)
	{
		size_t n;
		in >> n;
		in >> nodes[i][0] >> nodes[i][1]; //координаты узлов
	}

	if (!in) throw std::logic_error("reading nodes error");

	std::cout << "num nodes = " << num_nodes << std::endl;

	//READING ELEMS
	size_t num_elems = 0;
	in >> num_elems;

	auto& elems = mesh.elems;

	elems.resize(num_elems);

	for (size_t i = 0; i < num_elems; ++i)
	{
		//size_t n;
		//in >> n;
		in >> elems[i][0] >> elems[i][1] >> elems[i][2]; //номера узлов
		--elems[i][0]; --elems[i][1]; --elems[i][2]; //в файле нумерация от единицы
		//
		//геометрические характеристики
		auto [I, J, K] = elems[i];
		const double x[3] = { nodes[I][0], nodes[J][0], nodes[K][0] };
		const double y[3] = { nodes[I][1], nodes[J][1], nodes[K][1] };

		const double b[3] = { y[1] - y[2], y[2] - y[0], y[0] - y[1] };
		const double c[3] = { x[2] - x[1], x[0] - x[2], x[1] - x[0] };
		const double A2 = c[2] * b[1] - b[2] * c[1];

		if (A2 < 0) std::cout << "WARNING: negative jacobian" << i + 1 << ' ' << I + 1 << ' ' << J + 1 << ' ' << K + 1 << ' ' << "jac = " << A2 << std::endl;
		if (abs(A2) < SMALL_JAC) std::cout << "WARNING: small jacobian" << i + 1 << ' ' << I + 1 << ' ' << J + 1 << ' ' << K + 1 << ' ' << "jac = " << A2 << std::endl;
	}

	if (!in) throw std::logic_error("reading elements error");

	std::cout << "num elems = " << num_elems << std::endl;

	//READING BOUNDS
	size_t num_bounds = 0; //количество границ
	in >> num_bounds;

	auto& bounds = mesh.bounds;

	bounds.resize(num_bounds);

	for (size_t i = 0; i < num_bounds; ++i)
	{
		size_t num_bound_nodes = 0;
		in >> num_bound_nodes; //количество узлов на каждой из границ

		bounds[i].resize(num_bound_nodes);

		for (size_t j = 0; j < num_bound_nodes; ++j)
		{
			in >> bounds[i][j];
			--bounds[i][j]; //в файле нумерация от единицы
		}
	}

	if (!in) throw std::logic_error("reading bounds error");

	std::cout << "num bounds = " << num_bounds << std::endl;

	if (0) for (auto& bound : bounds)
		std::cout << '\t' << "num = " << bound.size() << std::endl;

	return mesh;
}

mesh_t read_mesh(const char* file_name = "mesh.txt")
{
	std::ifstream in(file_name);
	if (!in) throw std::logic_error("mesh file not found");
	return read_mesh(in);
}

template<typename T>
T zero()
{
	return T::Zero();
}

template<>
double zero<double>()
{
	return 0.0;
}

template<>
std::array<double, 2> zero<std::array<double, 2>>()
{
	return { 0,0 };
}

void operator+=(std::array<double, 2>& a, const std::array<double, 2>& b)
{
	a[0] += b[0];
	a[1] += b[1];
}

void operator/=(std::array<double, 2>& a, double val)
{
	a[0] /= val;
	a[1] /= val;
}

template<typename T>
auto to_nodes(const mesh_t& mesh, const std::vector<T>& edata, std::vector<T>& ndata)
{
	size_t N = mesh.nodes.size();

	ndata.resize(N, zero<T>());

	std::vector<size_t> count(N, 0);

	size_t el = 0;
	for (auto [I, J, K] : mesh.elems)
	{
		ndata[I] += edata[el];
		++count[I];
		ndata[J] += edata[el];
		++count[J];
		ndata[K] += edata[el];
		++count[K];

		++el;
	}

	size_t i = 0;
	for (auto& val : ndata)
	{
		val /= double(count[i]);
		++i;
	}
}

struct material_t
{
	double E = 2e11, mu = 0.3, density = 7800;
	double alpha = 1e-5; //температурное расширение
	double K = 1.0;
	double C = 1.0;
};

struct problem_t
{
	mesh_t mesh;
	material_t material;
};

using g_t = std::array<double, 2>; //grad для теплопроводности, можно заменить на Eigen

namespace thermal2dfullmatrix
{
	struct solver_t
	{
		Eigen::MatrixXd K;
		Eigen::VectorXd rhs;


		template<typename KL, typename DS>
		void to_K(const KL& k_loc, const DS& dofs)
		{
			for (size_t j = 0; j < k_loc.cols(); ++j)
				for (size_t i = 0; i < k_loc.rows(); ++i)
					K(dofs[i], dofs[j]) += k_loc(i, j);
		}


		template<typename RL, typename DS>
		void to_rhs(const RL& r_loc, const DS& dofs)
		{
			for (size_t i = 0; i < r_loc.size(); ++i)
				rhs(dofs[i]) += r_loc(i);
		}


		void to_K(double K_loc, size_t I)
		{
			K(I, I) += K_loc;
		}
		void to_rhs(double r_loc, size_t I)
		{
			rhs(I) += r_loc;
		}
	};

	double heatGen(double x, double y)
	{
		//const double ta = 0.01;
		//if (g_time < ta)
		//	return g_time / ta*heat_gen;

		const double heat_gen = 1;

		return heat_gen;
	}

	void assemble(const problem_t& problem, solver_t& solver)
	{
		auto& nodes = problem.mesh.nodes;

		size_t el = 0; //номер элемента
		for (auto [I, J, K] : problem.mesh.elems)
		{
			//геометрические характеристики
			const double x[3] = { nodes[I][0], nodes[J][0], nodes[K][0] };
			const double y[3] = { nodes[I][1], nodes[J][1], nodes[K][1] };

			const double b[3] = { y[1] - y[2], y[2] - y[0], y[0] - y[1] };
			const double c[3] = { x[2] - x[1], x[0] - x[2], x[1] - x[0] };
			const double A2 = c[2] * b[1] - b[2] * c[1];

			//матрица градиентов
			Eigen::Matrix<double, 2, 3> B;
			B <<
				b[0], b[1], b[2],
				c[0], c[1], c[2];

			//свойства материала

			const double D = problem.material.K;

			Eigen::Matrix<double, 3, 2> hBD = 0.5 * B.transpose() * D;

			//локальная матрица 
			Eigen::Matrix<double, 3, 3> k_loc;
			k_loc = hBD * B / A2; // 0.5*A2/A2/A2

			//сборка
			const size_t dofs[] = { I, J, K };

			solver.to_K(k_loc, dofs);

			//правая часть
			Eigen::Matrix<double, 3, 1> f_loc;

			const double xm = (x[0] + x[1] + x[2]) / 3;
			const double ym = (y[0] + y[1] + y[2]) / 3;
			const double Q = heatGen(xm, ym) * A2 / 6.;

			f_loc.array() = Q;

			//сборка правой части
			solver.to_rhs(f_loc, dofs);

			++el;
		}
	}

	//POSTPROCESS GRADS AND FLOWS
	void postprocess(const problem_t& problem, const Eigen::VectorXd& us, std::vector<g_t>& qs, std::vector<g_t>& grads)
	{
		size_t Ne = problem.mesh.elems.size();

		qs.resize(Ne);
		grads.resize(Ne);
		//
		auto& nodes = problem.mesh.nodes;
		//
		size_t el = 0; //номер элемента
		for (auto [I, J, K] : problem.mesh.elems)
		{
			const double x[3] = { nodes[I][0], nodes[J][0], nodes[K][0] };
			const double y[3] = { nodes[I][1], nodes[J][1], nodes[K][1] };

			const double b[3] = { y[1] - y[2], y[2] - y[0], y[0] - y[1] };
			const double c[3] = { x[2] - x[1], x[0] - x[2], x[1] - x[0] };
			const double A2 = c[2] * b[1] - b[2] * c[1];

			const double u[] = { us[I], us[J], us[K] };

			auto& gs = grads[el];

			gs[0] = (b[0] * u[0] + b[1] * u[1] + b[2] * u[2]) / A2;
			gs[1] = (c[0] * u[0] + c[1] * u[1] + c[2] * u[2]) / A2;

			const double D = problem.material.K;

			auto& q = qs[el];

			q[0] = -D * gs[0];
			q[1] = -D * gs[1];

			//
			++el;
		}

	}

	//BOUNDARY CONDITIONS

	void apply_flow(double Q, const std::vector<size_t>& bound, const mesh_t& mesh, solver_t& solver)
	{
		auto& nodes = mesh.nodes;
		for (size_t i = 0, e = bound.size() - 1; i < e; ++i)
		{
			const size_t I = bound[i];
			const size_t J = bound[i + 1];
			const auto [x1, y1] = nodes[I];
			const auto [x2, y2] = nodes[J];

			const double dx = x2 - x1; //ny
			const double dy = y2 - y1; //-nx
			const double len = sqrt(dx * dx + dy * dy);

			const double r_loc = 0.5 * Q * len;

			solver.to_rhs(r_loc, I);
			solver.to_rhs(r_loc, J);
		}
	}

	void apply_conv(double alpha, double Tinf, const std::vector<size_t>& bound, const mesh_t& mesh, solver_t& solver)
	{
		auto& nodes = mesh.nodes;

		for (size_t i = 0, e = bound.size() - 1; i < e; ++i)
		{
			const size_t I = bound[i];
			const size_t J = bound[i + 1];
			const auto [x1, y1] = nodes[I];
			const auto [x2, y2] = nodes[J];

			const double dx = x2 - x1; //ny*len == c[]
			const double dy = y2 - y1; //-nx*len == -b[]
			const double len = sqrt(dx * dx + dy * dy);

			const double r3 = alpha * len / 3.;
			const double r6 = r3 / 2;

			Eigen::Matrix2d k_loc;
			k_loc << r3, r6,
				r6, r3;

			const size_t dofs[] = { I, J };

			solver.to_K(k_loc, dofs);

			const double r_loc = 0.5 * alpha * Tinf * len;

			solver.to_rhs(r_loc, I);
			solver.to_rhs(r_loc, J);
		}
	}

	const double M = 1e20; //параметр штрафа

	void fix_T(solver_t& solver, size_t I, double u0 = 0.0)
	{
		const size_t N = I;
		solver.to_K(M, N);
		solver.to_rhs(M * u0, N);
	}

	// ---------------- Новый класс: нестационарный решатель ----------------
	struct unsteady_solver_t : solver_t
	{
		Eigen::VectorXd Mdiag;  // диагонализированная матрица масс (вектор)

		void init(size_t Nd)
		{
			K = Eigen::MatrixXd::Zero(Nd, Nd);
			rhs = Eigen::VectorXd::Zero(Nd);
			Mdiag = Eigen::VectorXd::Zero(Nd);
		}

		// Сборка диагональной матрицы масс (lumped mass)
		void assemble_mass(const problem_t& problem)
		{
			auto& nodes = problem.mesh.nodes;

			for (auto [I, J, Knode] : problem.mesh.elems)
			{
				const double x[3] = { nodes[I][0], nodes[J][0], nodes[Knode][0] };
				const double y[3] = { nodes[I][1], nodes[J][1], nodes[Knode][1] };

				const double b[3] = { y[1] - y[2], y[2] - y[0], y[0] - y[1] };
				const double c[3] = { x[2] - x[1], x[0] - x[2], x[1] - x[0] };
				const double A2 = c[2] * b[1] - b[2] * c[1];
				double area = 0.5 * std::abs(A2);

				double rho = problem.material.density;
				double Cp = problem.material.C;

				double val = rho * Cp * area / 3.0; // каждая вершина получает 1/3 от площади

				Mdiag[I] += val;
				Mdiag[J] += val;
				Mdiag[Knode] += val;
			}
		}

		// Добавление локальной матрицы в глобальную M
		template<typename ML, typename DS>
		void add_to_matrix(const ML& m_loc, const DS& dofs, Eigen::MatrixXd& Mglob)
		{
			for (size_t j = 0; j < m_loc.cols(); ++j)
				for (size_t i = 0; i < m_loc.rows(); ++i)
					Mglob(dofs[i], dofs[j]) += m_loc(i, j);
		}

		// Один шаг по времени (неявный Эйлер с диагонализированным [C])
		Eigen::VectorXd time_step(
			const Eigen::VectorXd& T_prev,
			double dt)
		{
			size_t Nd = T_prev.size();
			Eigen::MatrixXd A = K;   // начнём с K
			Eigen::VectorXd b = rhs; // начнём с F

			// [A] = (1/dt)*Cdiag + K
			// {b} = F + (1/dt)*Cdiag*T^n
			for (size_t i = 0; i < Nd; ++i)
			{
				A(i, i) += Mdiag[i] / dt;
				b[i] += (Mdiag[i] / dt) * T_prev[i];
			}

			Eigen::VectorXd T_new = A.lu().solve(b);
			return T_new;
		}
	};
}

// ---------- МОДУЛЬ ДЛЯ ОЦЕНКИ ОШИБОК ----------------

namespace error_analysis
{
	// Аналитическое решение
	double exact_solution(double x, double y, double t)
	{
		return std::exp(-2.0 * M_PI * M_PI * t) * std::cos(M_PI * x) * std::cos(M_PI * y);
	}

	// Норма L2 и C для численного решения
	struct error_norms_t {
		double L2_abs;
		double L2_rel;
		double C_abs;
		double C_rel;
	};

	// Вычисление норм ошибок
	error_norms_t compute_error(
		const problem_t& problem,
		const Eigen::VectorXd& Th,  // численное решение в узлах
		double time)
	{
		size_t Nd = problem.mesh.nodes.size();
		size_t Ne = problem.mesh.elems.size();

		double errL2_sq = 0.0;
		double normL2_sq = 0.0;
		double errC = 0.0;
		double normC = 0.0;

		// ---- L∞ норма (узловые значения) ----
		for (size_t i = 0; i < Nd; ++i) {
			double x = problem.mesh.nodes[i][0];
			double y = problem.mesh.nodes[i][1];
			double T_exact = exact_solution(x, y, time);
			double diff = std::abs(Th[i] - T_exact);

			errC = std::max(errC, diff);
			normC = std::max(normC, std::abs(T_exact));
		}

		// ---- L2 норма (по элементам, через центр тяжести) ----
		for (size_t e = 0; e < Ne; ++e) {
			auto [I, J, K] = problem.mesh.elems[e];
			auto& nI = problem.mesh.nodes[I];
			auto& nJ = problem.mesh.nodes[J];
			auto& nK = problem.mesh.nodes[K];

			// площадь треугольника
			double A = 0.5 * std::abs(
				(nJ[0] - nI[0]) * (nK[1] - nI[1]) -
				(nK[0] - nI[0]) * (nJ[1] - nI[1]));

			// центр тяжести
			double xm = (nI[0] + nJ[0] + nK[0]) / 3.0;
			double ym = (nI[1] + nJ[1] + nK[1]) / 3.0;

			// численное и точное в центре
			double Th_val = (Th[I] + Th[J] + Th[K]) / 3.0;
			double Te_val = exact_solution(xm, ym, time);

			double diff = Th_val - Te_val;

			errL2_sq += diff * diff * A;
			normL2_sq += Te_val * Te_val * A;
		}

		error_norms_t norms;
		norms.L2_abs = std::sqrt(errL2_sq);
		norms.L2_rel = norms.L2_abs / std::sqrt(normL2_sq);
		norms.C_abs = errC;
		norms.C_rel = errC / normC;

		return norms;
	}
}

namespace thermal_solver
{
	using namespace thermal2dfullmatrix;

	void apply_bc(problem_t& problem, solver_t& solver)
	{
		/*size_t fixed_node = 0;
		double T = 30;
		double alpha = 100;

		apply_conv(alpha, T, problem.mesh.bounds[1], problem.mesh, solver);

		double T0 = 20;

		for (size_t node : problem.mesh.bounds[3])
		{
			fix_T(solver, node, T0);
		}*/
	}

	void run()
	{
		std::cout << "run stud_fem thermal" << std::endl;
		problem_t problem;

		VAR(problem.material.K, "K");
		VAR(problem.material.C, "C");

		read_input();
		write_vars(std::cout);

		problem.mesh = read_mesh("mesh1.txt");

		size_t Nd = problem.mesh.nodes.size();

		solver_t s;
		s.K = Eigen::MatrixXd::Zero(Nd, Nd);
		s.rhs = Eigen::VectorXd::Zero(Nd);

		const bool PRINT = true;
#define print if(PRINT) std::cout

		assemble(problem, s);
		print << s.K << endl << endl;

		apply_bc(problem, s);
		print << s.K << endl;
		print << s.rhs << endl << endl;

		Eigen::VectorXd us;
		us = s.K.lu().solve(s.rhs);
		print << us << endl;

#undef print
		//результаты по элементам
		std::vector<g_t> qs; //flows
		std::vector<g_t> gs; //grads

		postprocess(problem, us, qs, gs);

		//результаты по узлам
		std::vector<g_t> qsn;
		std::vector<g_t> gsn;

		auto write_nodes_results = [&](std::ofstream& out)
			{
				const auto endl = '\n';

				size_t NN = problem.mesh.nodes.size(); //us.size()
				out << NN << endl;

				out << "T" << endl;
				for (size_t i = 0, e = NN; i < e; i += 1) out << us[i] << endl;

				write_vec_vec(out, gsn, { "Gx", "Gy" });
				write_vec_vec(out, qsn, { "Qx", "Qy" });

			};

		auto write_elems_results = [&](std::ofstream& out)
			{
				const auto endl = '\n';

				out << gs.size() << endl;

				write_vec_vec(out, gs, { "el_Gx", "el_Gy" });

				write_vec_vec(out, qs, { "el_Qx", "el_Qy" });

			};

		auto postproces_and_write = [&](auto problem_name)
			{
				std::ofstream out(std::string("results") + problem_name + ".txt");
				write_nodes_results(out);

				to_nodes(problem.mesh, qs, qsn);
				to_nodes(problem.mesh, gs, gsn);

				std::ofstream e_out(std::string("elem_results") + problem_name + ".txt");
				write_elems_results(e_out);
			};

		postproces_and_write("_therm");
	}

	// нестационарный случай
	void run_unsteady()
	{
		std::cout << "run unsteady thermal" << std::endl;
		problem_t problem;

		VAR(problem.material.K, "K");
		VAR(problem.material.C, "C");
		VAR(problem.material.density, "rho");

		read_input();
		write_vars(std::cout);

		problem.mesh = read_mesh("mesh1.txt");
		size_t Nd = problem.mesh.nodes.size();
		size_t Ne = problem.mesh.elems.size();

		unsteady_solver_t s;
		s.init(Nd);

		assemble(problem, s);
		s.assemble_mass(problem);
		apply_bc(problem, s);

		Eigen::VectorXd T(Nd);
		for (size_t i = 0; i < Nd; ++i) {
			double x = problem.mesh.nodes[i][0];
			double y = problem.mesh.nodes[i][1];
			T[i] = cos(M_PI * x) * cos(M_PI * y);  // 2D версия условия
		};

		double dt = 0.01;
		int Nt = 24;

		/*double dt = 0.000025;
		int Nt = 100;*/

		/*double dt = 0.00000625;
		int Nt = 400;*/

		// Подготовка файлов для Tecplot
		std::ofstream tec_out("results_unsteady.dat");
		if (!tec_out) {
			throw std::logic_error("Cannot open results_unsteady.dat for writing");
		}

		// Заголовок Tecplot
		tec_out << "TITLE = \"Unsteady thermal results\"\n";
		tec_out << "VARIABLES = \"X\",\"Y\",\"T\",\"Gx\",\"Gy\",\"Qx\",\"Qy\"\n";
		tec_out << std::flush;

		// Вспомогательные контейнеры для постпроцессинга
		std::vector<g_t> elem_qs; // по элементам
		std::vector<g_t> elem_gs; // по элементам
		std::vector<g_t> node_qs; // по узлам
		std::vector<g_t> node_gs; // по узлам

		// Для каждого шага времени: делаем шаг, постпроцессим и пишем ZONE
		for (int n = 0; n < Nt; ++n)
		{
			T = s.time_step(T, dt);
			double current_time = (n + 1) * dt; // время после шага

			std::cout << "step " << n
				<< ", T_max = " << T.maxCoeff()
				<< ", time = " << current_time
				<< std::endl;
			if (current_time == 0.0025) {

				auto norms = error_analysis::compute_error(problem, T, current_time);
				std::cout << "Errors at t=" << current_time
					<< ": L2_abs=" << norms.L2_abs
					<< ", L2_rel=" << norms.L2_rel
					<< ", C_abs=" << norms.C_abs
					<< ", C_rel=" << norms.C_rel
					<< std::endl;
			};

			// Постпроцесс: градиенты и потоки по элементам
			postprocess(problem, T, elem_qs, elem_gs);

			// Перенос значений по элементам на узлы (усреднение)
			to_nodes(problem.mesh, elem_qs, node_qs);
			to_nodes(problem.mesh, elem_gs, node_gs);
			
			// Запись зоны для текущего времени
			// N = number of nodes, E = number of elements
			tec_out << "ZONE T=\"Time=" << current_time << "\" "
				<< "STRANDID=1 SOLUTIONTIME=" << current_time << " "
				<< "DATAPACKING=POINT ZONETYPE=FETRIANGLE "
				<< "N=" << Nd << " E=" << Ne << "\n";

			// 1) данные по узлам: X Y T Gx Gy Qx Qy
			for (size_t i = 0; i < Nd; ++i)
			{
				tec_out
					<< problem.mesh.nodes[i][0] << " "  // X
					<< problem.mesh.nodes[i][1] << " "  // Y
					<< T[i] << " ";                     // T

				// безопасно: если node_gs/node_qs пусты, выведем 0
				if (i < node_gs.size()) tec_out << node_gs[i][0] << " " << node_gs[i][1] << " ";
				else tec_out << 0.0 << " " << 0.0 << " ";

				if (i < node_qs.size()) tec_out << node_qs[i][0] << " " << node_qs[i][1] << "\n";
				else tec_out << 0.0 << " " << 0.0 << "\n";
			}

			// 2) затем connectivity для треугольников (1-based indexing!)
			for (size_t e = 0; e < Ne; ++e)
			{
				const auto& el = problem.mesh.elems[e];
				tec_out << (el[0] + 1) << " " << (el[1] + 1) << " " << (el[2] + 1) << "\n";
			}

			tec_out << std::flush;
		}

		tec_out.close();
		std::cout << "Unsteady results written to results_unsteady.dat\n";
	}


}


void run()
{
	// thermal_solver::run(); // стационарный
	thermal_solver::run_unsteady(); // нестационарный
}

int main()
{
	try
	{
		run();
	}
	catch (std::logic_error e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}

void skip_line(std::istream& in)
{
	in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

bool is_valid_first_name_char(char c)
{
	static std::locale loc;
	return std::isalpha(c, loc) || c == '_';
}

bool is_valid_name_char(char c)
{
	static std::locale loc;
	return std::isalnum(c, loc) || c == '_';// || c=='.';
}

std::string read_name(std::istream& in)
{
	std::string name;
	char c = 0;
	in >> c; //пропускает пробелы
	if (in && is_valid_first_name_char(c))
	{
		do {
			name += c;
			in.get(c);
		} while (in && is_valid_name_char(c));
		if (!in.eof()) in.putback(c);
	}
	else
	{
		if (!in.eof()) in.putback(c);
	}
	return name;
}

std::map<std::string, double*> vars;

void VAR(double& var, const char* name)
{
	vars[name] = &var;
}

void read_vars(std::istream& in)
{
	while (in)
	{
		auto name = read_name(in);

		char sym;

		if (!(in >> sym) || sym != '=') break;

		if (auto it = vars.find(name); it != vars.end())
		{
			double x;
			if (in >> x)
			{
				*(it->second) = x;
			}
			else
			{
				skip_line(in); //не удалось считать значение - пропускаем до конца строки
			}
		}
		else
		{
			skip_line(in); //неизвестное имя переменной - пропускаем до конца строки
		}
	}
}

void write_vars(std::ostream& out)
{

	for (auto& [name, var] : vars)
	{
		out << name << " = " << *var << endl;
	}
}



