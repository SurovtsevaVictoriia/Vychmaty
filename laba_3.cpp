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


class Solver {

public:
	
	double a;
	double b;
	double y0;
	double e;

	int k; // порядок сходимости
	int N; // количество точек

	double f(double x, double y) {
		return(-2 - 4 * x * y) / (x * x) - y * y;
	}


	Solver(double _a, double _b, double _y0, double _e, int _k, int _N) {
		a = _a;
		b = _b;
		y0 = _y0;
		e = _e;		
		k = _k;
		N = _N;
	}

	class Grid_function {//contains the grid dunction itself
	public:

		int L;
		double h;
		std::vector<std::vector<double>> grid;

		Grid_function() {};

		Grid_function(Solver & s, int _L) {
			L = _L;
			h = (s.b - s.a) / L;
			grid.push_back(x_axis_grid(s));
			grid.push_back(y_axis_grid(s));
		}		

		//returns vector of x knots
		std::vector<double> x_axis_grid(Solver & s) {

			std::vector<double> buffer;
			for (int i = 0; i <= L; i++) {
				buffer.push_back(s.a + h * i);
			}
			return buffer;
		}

		//returns vector of y knots
		std::vector<double> y_axis_grid( Solver & s) {

			std::vector<double> buffer;
			buffer.push_back(s.y0);
			for (int i = 1; i <= L; i++) {

				double f1 = s.f(grid[0][i - 1], buffer[i - 1]);
				double f2 = s.f(grid[0][i - 1] + h, buffer[i - 1] + h * f1);

				double y = buffer[i - 1] + h / 2 * (f1 + f2);

				buffer.push_back(y);
			}		
			return buffer;
		};

	};




	std::vector<std::vector<double >> result() {		

		double L0 = N; 
		Grid_function grid_2h;
		Grid_function grid_h;

		grid_2h = Grid_function(*this, L0); //less dots
		grid_h = Grid_function(*this, L0);//more dots
		std::vector<std::vector<double >> grid_h_to_2h ;

		grid_h_to_2h = grid_2h.grid;

		for (int i = 0; i < grid_2h.grid[0].size(); i++) {

			grid_h_to_2h[1][i] += 10 * e;
		
		}


		while (max_dif(grid_2h.grid[1], grid_h_to_2h[1]) / (2 * k - 1)> e ) {

			grid_h_to_2h.clear();
			grid_2h = grid_h; //less dots
			grid_h = Grid_function(*this, 2 * L0);//more dots			

			std::vector<double> buffer_x;
			std::vector<double> buffer_y;

			for (int i = 0; i < grid_2h.grid[0].size(); ++i) {

				buffer_x.push_back(grid_h.grid[0][i * 2]);
				buffer_y.push_back(grid_h.grid[1][i * 2]);

			}

			grid_h_to_2h.push_back(buffer_x);
			grid_h_to_2h.push_back(buffer_y);

			L0 *= 2;
		}

		std::cout << "result acquired!\n";

		return grid_h.grid;
	}

	std::vector<std::vector<double>> equidistant() {

		std::vector<std::vector<double >> buffer;
		std::vector<double > buffer_x;
		std::vector<double> buffer_y;

		
		std::vector<std::vector<double >> res = this->result();
		int n = res[0].size() - 1;

		for (int i = 0; i <= N; i++) {
			buffer_x.push_back(res[0][i * n / N]);
			buffer_y.push_back(res[1][i * n / N]);
		}

		buffer.push_back(buffer_x);
		buffer.push_back(buffer_y);

		return buffer;

	}
};




void get_dots(std::vector<std::vector<double>> vec) {//works

	using namespace std;

	string outFile = "grid.csv";
	ofstream out;
	out.open(outFile);

	if (!out.is_open()) {
		cerr << "cannot /*not a comment*/ open outFile"; /*<< outFile << endl;
		return EXIT_FAILURE;*/ int a = 0; //comment
	}

	out << "x;y;\n";

	for (int i = 0; i < vec[0].size(); ++i) {

		out << vec[0][i]<< ";" << vec[1][i] << ";" << "\n";

	}

	out.close();
}

int main() {	

	
	Solver s = Solver(1, 2, -1, 0.0001, 2, 12);
	
	std::vector<std::vector<double>> final = s.equidistant();
	get_dots(final);
	
}
