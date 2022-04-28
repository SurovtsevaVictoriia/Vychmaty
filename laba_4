#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <locale>
#include<string>

void print(std::vector<double> vec) {


	for (int i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << " ";
	}

	std::cout << "\n" << std::endl;
}

double max_dif(std::vector<double> a, std::vector<double >b) {

	double dif = abs(a[0] - b[0]);
	if (a.size() != b.size()) {
		std::cout << "sizes incompstible \n";
		return -1;
	}
	else {
		for (int i = 1; i < a.size(); i++) {
			if (abs(a[i] - b[i]) > dif) dif = abs(a[i] - b[i]);
		}
	}

	return dif;

}

void get_dots(std::vector<std::vector<double>> vec, std::string name) {//works

	using namespace std;

	string outFile = name + ".csv";
	ofstream out;
	out.open(outFile);

	if (!out.is_open()) {
		cerr << "cannot /*not a comment*/ open outFile"; /*<< outFile << endl;
		return EXIT_FAILURE;*/ int a = 0; //comment
	}

	out << "x;y;\n";

	for (int i = 0; i < vec[0].size(); ++i) {

		out << vec[0][i] << ";" << vec[1][i] << ";" << "\n";

	}

	out.close();
}

double k_left_f(double x) {
	return sin(x * x) + 1;
}
double k_right_f(double x) {
	return sin(x * x) + 1;
}
double q_left_f(double x) {
	return x;
}
double q_right_f(double x) {
	return x * x * x;
}
double f_left_f(double x) {
	return 1;
}
double f_right_f(double x) {
	return x * x - 1;
}


class Solver {
public:
	double u_0;
	double u_L;
	double x0;

	double a0;
	double b0;

	double e; 
	int N;

	Solver(double _u_0, double _u_L, double _x0, double _a0, double _b0, double _e, double _N  ) {
		u_0 = _u_0;
		u_L = _u_L;
		x0 = _x0;
		a0 = _a0;
		b0 = _b0;
		e = _e;
		N = _N;


	}
	class Grid {
	public:	
		int L;
		double h;
		std::vector<std::vector<double>> grid;


		std::vector<double> a;
		std::vector<double> b;
		std::vector<double> c;
		std::vector<double> d;
		std::vector<double> alpha;
		std::vector<double> beta;

		std::vector<double> u;
		Grid(Solver& s, Grid& g) {
			model_get_coeffs(s, g);
			u = get_u(s, g);
		}

		void model_get_coeffs(Solver& s, Grid& g) {

			int l_alpha = s.x0 / g.h;
			int l_beta = s.x0 / g.h + 1;

			a.push_back(0);
			b.push_back(0);
			c.push_back(0);
			d.push_back(0);

			for (int l = 1; l <= l_alpha - 1; l++) { //1 ... l_a - 1
				a.push_back(k_left_f(s.x0));
				b.push_back(-2 * k_left_f(s.x0) - q_left_f(s.x0) * g.h * g.h);
				c.push_back(k_left_f(s.x0));
				d.push_back(-f_left_f(s.x0) * g.h * g.h);
			}
			for (int l = l_alpha; l <= l_beta; l++) {
				a.push_back(0);
				b.push_back(0);
				c.push_back(0);
				d.push_back(0);
			}
			for (int l = l_beta + 1; l <= g.L - 1; l++) { // l_b + 1 ...  L - 1
				a.push_back(k_right_f(s.x0));
				b.push_back(-2 * k_right_f(s.x0) - q_right_f(s.x0) * g.h * g.h);
				c.push_back(k_right_f(s.x0));
				d.push_back(-f_right_f(s.x0) * g.h * g.h);
			}

			a.push_back(0);
			b.push_back(0);
			c.push_back(0);
			d.push_back(0);


			alpha.push_back(0);
			alpha.push_back(-a[1] / b[1]);
			for (int l = 2; l <= l_alpha - 1; l++) {
				alpha.push_back(-a[l] / (b[l] + c[l] * alpha[l - 1]));
			}
			for (int l = l_alpha; l <= g.L - 2; l++) {
				alpha.push_back(0);
			}
			alpha.push_back(-c[g.L - 1] / b[g.L - 1]);
			for (int l = g.L - 2; l >= l_beta + 1; l--) {
				alpha[l] = (-c[l] / (b[l] + a[l] * alpha[l + 1]));
			}
			alpha.push_back(0);


			beta.push_back(0);
			beta.push_back((d[1] - c[1] * s.u_0) / b[1]);
			for (int l = 2; l <= l_alpha - 1; l++) {
				beta.push_back((d[l] - c[l] * beta[l - 1]) / (b[l] + c[l] * alpha[l - 1]));
			}
			for (int l = l_alpha; l <= g.L - 2; l++) {
				beta.push_back(0);
			}
			beta.push_back((d[g.L - 1] - c[g.L - 1] * s.u_L) / b[g.L - 1]);
			for (int l = g.L - 2; l >= l_beta + 1; l--) {
				beta[l] = (d[l] - a[l] * beta[l + 1]) / (b[l] + a[l] * alpha[l + 1]);
			}
			beta.push_back(0);


		}
		void actual_get_coeffs(Solver& s, Grid& g) {

			int l_alpha = s.x0 / g.h;
			int l_beta = s.x0 / g.h + 1;

			a.push_back(0);
			b.push_back(0);
			c.push_back(0);
			d.push_back(0);


			for (int l = 1; l <= l_alpha - 1; l++) { //1 ...l_alpha - 1
				a.push_back(k_left_f(g.h * (l + 0.5)));
				b.push_back(-(k_left_f(g.h * (l + 0.5)) + k_left_f(g.h * (l - 0.5)) + q_left_f(g.h * l) * g.h * g.h));
				c.push_back(k_left_f(g.h * (l - 0.5)));
				d.push_back(-f_left_f(g.h * l) * g.h * g.h);
			}

			for (int i = 0; i < 2; i++) {// l_alpha, l_beta
				a.push_back(0);
				b.push_back(0);
				c.push_back(0);
				d.push_back(0);
			}

			for (int l = l_beta + 1; l <= g.L - 1; l++ ) {//l_beta + 1... L-1

				a.push_back(k_right_f(g.h * (l + 0.5)));
				b.push_back(-(k_right_f(g.h * (l + 0.5)) + k_right_f(g.h * (l - 0.5)) + q_right_f(g.h * l) * g.h * g.h));
				c.push_back(k_right_f(g.h * (l - 0.5)));
				d.push_back(-f_right_f(g.h * l) * g.h * g.h);
			}

			a.push_back(0);
			b.push_back(0);
			c.push_back(0);
			d.push_back(0);

			//now alpha and beta

			alpha.push_back(0);
			alpha.push_back(-a[1] / b[1]);
			for (int l = 2; l <= l_alpha - 1; l++) {
				alpha.push_back(-a[l] / (b[l] + c[l] * alpha[l - 1]));
			}
			for (int l = l_alpha; l <= g.L - 2; l++) {
				alpha.push_back(0);
			}
			alpha.push_back(-c[g.L - 1] / b[g.L - 1]);
			for (int l = g.L - 2; l >= s.x0 / g.h + 2; l--) {
				alpha[l] = (-c[l] / (b[l] + a[l] * alpha[l + 1]));
			}



			beta.push_back(0);
			beta.push_back((d[1] - c[1] * s.u_0) / b[1]);
			for (int l = 2; l <= l_alpha - 1; l++) {
				beta.push_back((d[l] - c[l] * beta[l - 1]) / (b[l] + c[l] * alpha[l - 1]));
			}
			for (int l = l_alpha; l <= g.L - 2; l++) {
				beta.push_back(0);
			}
			beta.push_back((d[g.L - 1] - c[g.L - 1] * s.u_L) / b[g.L - 1]);
			for (int l = g.L - 2; l >= s.x0 / g.h + 2; l--) {
				beta[l] = (d[l] - a[l] * beta[l + 1]) / (b[l] + a[l] * alpha[l + 1]);
			}







		}

		std::vector<double>get_u(Solver& s, Grid& g) {

			int l_alpha = s.x0 / g.h;
			int l_beta = s.x0 / g.h + 1;

			// 0. fill the vector
			std::vector<double> u;
			u.push_back(s.u_0);
			for (int l = 1; l <= g.L - 1; l++) {
				u.push_back(0);
			}
			u.push_back(s.u_L);

			// 1. get middle values
			u[l_alpha] = (k_left_f(s.x0) * beta[l_alpha - 1] + k_right_f(s.x0) * beta[l_beta + 1]) /
				(k_left_f(s.x0) * (1 - alpha[l_alpha - 1]) + k_right_f(s.x0) * (1 - alpha[l_beta + 1]));
			u[l_beta] = u[l_alpha];
			std::cout << "u[l_alpha] = " << u[l_alpha]  << "\n";
			u[l_alpha - 1] = alpha[l_alpha - 1] * u[l_alpha] + beta[l_alpha - 1];
			u[l_beta + 1] = alpha[l_beta + 1] * u[l_beta] + beta[l_beta + 1];

			// 2. everything else
			for (int l = l_alpha - 2; l >=1 ; l--) {
				u[l] = alpha[l] * u[l + 1] + beta[l];
			}
			for (int l = l_beta + 2; l <= g.L - 1; l++) {
				u[l] = alpha[l] * u[l - 1] + beta[l];
			}

			return u;
		}

		Grid(Solver& s, int _L, bool mode) {
			L = _L;
			h = (s.b0 - s.a0) / L;
			grid.push_back(x_axis_grid(s));

			if (mode) {
				actual_get_coeffs(s, *this);
			}
			else { model_get_coeffs(s, *this); }

			grid.push_back(get_u(s, *this));
		}
		Grid() {

		}
		//returns vector of x knots
		std::vector<double> x_axis_grid(Solver& s) {

			std::vector<double> buffer;
			for (int i = 0; i <= L; i++) {
				buffer.push_back(s.a0 + h * i);
			}
			return buffer;
		}

		


	};
	
	std::vector<std::vector<double >> result(bool mode) {

		int L0 = N - 1;
		Grid grid_2h;
		Grid grid_h;

		grid_2h = Grid(*this, L0, mode); //less dots
		grid_h = Grid(*this, L0, mode);//more dots
		std::vector<std::vector<double >> grid_h_to_2h;

		grid_h_to_2h = grid_2h.grid;

		for (int i = 0; i < grid_2h.grid[0].size(); i++) {

			grid_h_to_2h[1][i] += 10 * e;

		}

		std::cout << "max_dif = " << max_dif(grid_2h.grid[1], grid_h_to_2h[1]) << "L0 = " << L0 << "\n";
		while (max_dif(grid_2h.grid[1], grid_h_to_2h[1])  > e) {
			
			grid_h_to_2h.clear();
			grid_2h = grid_h; //less dots
			grid_h = Grid(*this, 2 * L0, mode);//more dots			

			std::vector<double> buffer_x;
			std::vector<double> buffer_y;

			for (int i = 0; i < grid_2h.grid[0].size(); ++i) {

				buffer_x.push_back(grid_h.grid[0][i * 2]);
				buffer_y.push_back(grid_h.grid[1][i * 2]);

			}

			grid_h_to_2h.push_back(buffer_x);
			grid_h_to_2h.push_back(buffer_y);

			L0 *= 2;
			std::cout << "max_dif = " << max_dif(grid_2h.grid[1], grid_h_to_2h[1]) << "L0 = " << L0 << "\n";
		}

		std::cout << "result acquired!\n";

		get_dots(grid_h.grid, "best");
		return grid_h.grid;		
	}

	std::vector<std::vector<double>> equidistant(bool mode) {

		std::vector<std::vector<double >> buffer;
		std::vector<double > buffer_x;
		std::vector<double> buffer_y;


		std::vector<std::vector<double >> res = this->result(mode);
		int n = res[0].size() - 1;

		for (int i = 0; i <= N - 1; i++) {
			buffer_x.push_back(res[0][i * n / ( N - 1)]);
			buffer_y.push_back(res[1][i * n / (N - 1)]);
		}

		buffer.push_back(buffer_x);
		buffer.push_back(buffer_y);

		return buffer;

	}
};

int main() {


	double u_0 = 1;
	double u_L = 2;
	double x0 = 1/ sqrt(5);

	double a0 = 0;
	double b0 = 1;

	double e = 0.0001;
	int N = 11;

	Solver s = Solver(u_0, u_L, x0, a0, b0, e, N);
	std::vector<std::vector<double>> final_model = s.equidistant(false);
	get_dots(final_model, "model");

	std::vector<std::vector<double>> final_actual = s.equidistant(true);
	get_dots(final_actual, "actual");

	std::cout << "model : \n";
	print(final_model[1]);
	std::cout << "actual : \n";
	print(final_actual[1]);

	std::cout << "done";


}
