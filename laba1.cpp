#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <locale>


void print_as_poly(std::vector<double> vec) {

	std::cout << "f(x) = ";
	for (int i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << "x^" << vec.size() - 1 - i << " + ";
	}

	std::cout << 0 << std::endl;
}
void print(std::vector<double> vec) {

	
	for (int i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << " ";
	}

	std::cout << "\n" << std::endl;
}

double power(double x, int n) {
	double value = 1;
	for (int i = 0; i < n; i++) {
		value *= x;
	}
	return value;
}


double get_value(std::vector<double> coef_, double x) {//from array
	double value = 0;
	for (int i = 0; i < coef_.size(); i++) {
		value += power(x, coef_.size() - 1 - i) * coef_[i];
	}
	return value;
}


class Spline {
public:

	int N;
	std::vector<std::vector<double>> values;
	Spline(std::vector<double> x_values_, std::vector<double> y_values_) {
		
		N = x_values_.size() - 1;
		
			values.push_back(x_values_);
			values.push_back(y_values_);
	
	}

	std::vector<double> new_poly;
	std::vector<std::vector<double>> residue;//(x-x0)(x-x1)(x-x2)
	std::vector<double> coefs;
	std::vector<double> final_coefs;
	void calculate_residue() {//works

		residue.push_back(std::vector<double>{ 1, -values[0][0] });		

		for (int i = 1; i < N + 1; i++) {

			std::vector<double> buffer;

			buffer.push_back(1);
			for (int j = 0; j < residue[i - 1].size() - 1;j++) {
				buffer.push_back(residue[i - 1][j] * (-values[0][i]) + residue[i - 1][j + 1]);
			}
			buffer.push_back(residue[i - 1].back() * (-values[0][i]));
			residue.push_back(buffer);			
		}
	}


	void calculate_coefs() {//recurrent shit for b
		for (int i = 0; i < N + 1; i++) {
			coefs.push_back(ff(i, i));
		}
	}

	
	//would be more effective if i would just calculate buffer vectors but idc
	double ff(int i, int j) {
		if (i == 0) { 
			double a = values[1][0];
			std::cout << " i = " << i << " j = " << j << "\t" << a << "\n";
			return values[1][0]; 
		}
		else if (i == 1) {
			double a = (values[1][j] - values[1][j - 1]) / (values[0][j] - values[0][j - 1]);
			std::cout <<" i = 1, j = " << j << "\t" << a << "\n";
			return a;
		}
		else {
			double a = (ff(i - 1, j) - ff(i - 1, j - 1)) /( values[0][j] - values[0][j - i]);
			std::cout << " i = " << i << " j = " << j << "\t" << a << "\n";
			return a ;
		}
		
	}

	void calculate_coefs_final() {//for 

		for (int i = 0; i < N; i++) {//

			double a = 0;
			for (int j = 0; j <= i; j++) {

				a += residue[N - j - 1][i - j] * coefs[N - j];
								
			}

			final_coefs.push_back(a);
		}

		double a = 0;
		for (int j = 0; j < N; j++) {

			a += residue[N - j - 1][N - j] * coefs[N - j];

		}
		a += coefs[0];
		final_coefs.push_back(a);
	}
};


void get_dots(Spline s) {


	using namespace std;
	
	cout << "placing dots";

	string outFile = "out.csv";
	string inFile = "in.csv";

	//ofstream out(outFile);
	ofstream out;

	//out.imbue(std::locale("German_germany"));

	out.open(inFile);
	if (!out.is_open()) {
		cerr << "cannot /*not a comment*/ open outFile"; /*<< outFile << endl;
		return EXIT_FAILURE;*/ int a = 0; //comment
	}

	out << "x;y;\n";
	for (int i = 0; i <= s.N; i++) {
		out << s.values[0][i] << ";" << s.values[1][i] << ";" << "\n";
	}
	out.close();

	out.open(outFile);

	if (!out.is_open()) {
		cerr << "cannot /*not a comment*/ open outFile"; /*<< outFile << endl;
		return EXIT_FAILURE;*/ int a = 0; //comment
	}

	out << "xx;yy;\n";
	double start = s.values[0][0] - 1;
	double end = s.values[0][s.N] + 1;
	double step = 0.005;

	for (auto i = start; i <= (end / step); ++i) {
		out << i * step << ";" << get_value(s.final_coefs, i * step) << ";" << "\n";
		
	}
	out.close();

	
}


void do_the_thing(std::vector<double> &x, std::vector<double>& y,  Spline & s) {
	
	//сначала делаем по всем, потом по трем соседним
	s.calculate_residue();
	s.calculate_coefs();

	s.calculate_coefs_final();

	print_as_poly(s.residue[2]);
	print(s.coefs);
	print_as_poly(s.final_coefs);

	get_dots(s);
}

int main() {

	std::vector<double> x_test = { 1, 2, 4 };
	std::vector<double> y_test = { 1, 2, 3 };

	std::vector<double> x_test2 = { 0, 1, 2, 3, 4 };
	std::vector<double> y_test2 = { 1, 5, 33, 109, 257 };

	std::vector<double> x_test3 = { 1, 2, 3, 4, 5};
	std::vector<double> y_test3 = { 5, 3, 8, 2, 14 };
	

	std::vector<double> x = { 0.17453, 0.5236, 0.87267, 1.22173, 1.5708, 1.91986 };
	std::vector<double> y = { 0.000038, 0.00052, 0.00254, 0.01013, 0.03636, 0.11699 };

	Spline s = Spline(x, y);
	do_the_thing(x, y, s);
	
	std::cout << "s";

}
