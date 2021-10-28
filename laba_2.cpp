#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <locale>
#include<string>


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

int factorial(int n) {

	switch (n)
	{
	case 1:
		return 1;
		break;
	case 2:
		return 2;
		break;
	case 3:
		return 6;
		break;
	case 4:
		return 24;
		break;
	case 5:
		return 120;
		break;
	case 6:
		return 720;
		break;

	}


	int a = 1;
	for (int i = 1; i <= n; i++) {
		a = a * i;
	}
	return a;
}


double get_value(std::vector<double> coef_, double x) {//from array
	double value = 0;
	for (int i = 0; i < coef_.size(); i++) {
		value += power(x, coef_.size() - 1 - i) * coef_[i];
	}
	return value;
}

auto get_derivative(std::vector<double> &coef_) {//derivative from array

	
	int size = coef_.size() - 1;
	std::vector<double> array_(size);

	for (int i = 0; i < size; i++) {
		array_[i] = coef_[i] * double(factorial(coef_.size() - 1 - i)) / double(factorial(coef_.size() - 1 - i - 1));
	}


	return array_;
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

	
	std::vector<std::vector<double>> residue;//(x-x0)(x-x1)(x-x2)
	std::vector<double> coefs;//b
	std::vector<double> final_coefs;//extrapolating polynom
	std::vector<double> derivative;//derivative
	std::vector<std::vector<double>> extra_coefs;//on two dots[[0123][0123]]...

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

	void calculate_extra_coefs() {
		std::cout << "calculating extra coefs\n";
		for (int i = 0; i < N; i++) {//number of intervals
			std::vector<double> buffer;

			for (int j = 0; j < 4; j++) {
				std::cout << "a: i = " << i << " j = " << j << " \t" << a(i, j) << "\n";
				buffer.push_back(a(i, j));
			}

			extra_coefs.push_back(buffer);
			
			print(buffer);
		}
	}

	double a(int i, int j) {
		std::vector<double> x = values[0];
		std::vector<double> y = values[1];
		switch (j) {
			case 0:
				return (get_value(derivative, x[i + 1]) * (x[i + 1] - x[i])
					- 2 * (y[i + 1] - y[i])
					+ get_value(derivative, x[i]) * (x[i + 1] - x[i]))
					/ power(x[i + 1] - x[i], 3);
				break;

			case 1:
				return  (-get_value(derivative, x[i + 1]) * (x[i + 1] - x[i]) * (x[i + 1] + 2 * x[i])
					+ 3 * (y[i + 1] - y[i]) * (x[i + 1] + x[i])
					- get_value(derivative, x[i]) * (x[i + 1] - x[i]) * (x[i] + 2 * x[i+1]))
					/ power(x[i + 1] - x[i], 3);
				break;

			case 2:
				return  (get_value(derivative, x[i + 1]) * x[i] * (2 * x[i + 1] + x[i]) * (x[i + 1] - x[i])
					- 6 * (y[i + 1] - y[i]) * x[i] * x[i + 1]
					+ get_value(derivative, x[i]) * x[i + 1] * (x[i + 1] + 2 * x[i]) * (x[i + 1] - x[i]))
					/ power(x[i + 1] - x[i], 3);
				break;
			case 3:
				return (-get_value(derivative, x[i + 1]) * x[i] * x[i] * x[i + 1] * (x[i + 1] - x[i])
					+ y[i + 1] * x[i] * x[i] * (3 * x[i + 1] - x[i])
					+ y[i] * x[i + 1] * x[i + 1] * (x[i + 1] - 3 * x[i])
					- get_value(derivative, x[i]) * x[i] * x[i + 1] * x[i + 1] * (x[i + 1] - x[i]))
					/ power(x[i + 1] - x[i], 3);
				break;
		}
	}
};


void get_dots(Spline s) {


	using namespace std;
	
	cout << "placing dots\n";

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
	double start = s.values[0][0] ;
	double end = s.values[0][s.N] ;
	double step = 0.05;

	std::cout << "start: " << start << " \t end: " << end << "\n";
	for (int i = (start / step) - 1; i <= (end / step) + 1; ++i) {
		std::cout << i << " ";
		out << i * step << ";" << get_value(s.final_coefs, i * step) << ";" << "\n";
		
	}
	std::cout << "\n";
	out.close();

	
}

void get_more_dots(Spline s) {

	using namespace std;
	ofstream out;

	for (int i = 0; i < s.N; i++) {

		string outFile = "out" + to_string(i) + ".csv";
		out.open(outFile);

		out << "xxx" + to_string(i) + ";yyy" + to_string(i)+ ";\n";
		double step = 0.05;
		double start = s.values[0][i];
		double end = s.values[0][i + 1];
		std::cout << "start: " << start << " \t end: " << end << "\n";


		for (int j = (start / step)  -1 ; j <= (end / step) + 1 ; ++j) {
			std::cout << j << " ";
			out << j * step << ";" << get_value(s.extra_coefs[i], j * step) << ";" << "\n";
			
		}
		std::cout << "\n";
		out.close();
	}

	//out.close();
	


}
void do_the_thing(Spline & s) {//on all
	
	//сначала делаем по всем, потом по трем соседним
	s.calculate_residue();
	s.calculate_coefs();
	s.calculate_coefs_final();
	
	//print_as_poly(s.residue[2]);
	print(s.coefs);
	print_as_poly(s.final_coefs);

	get_dots(s);
}

void do_another_thing(Spline & s) {//on two 

	s.derivative = get_derivative(s.final_coefs);
	print_as_poly(s.derivative);
	s.calculate_extra_coefs();

	get_more_dots(s);


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

	do_the_thing(s);
	do_another_thing(s);
	
	std::cout << "s";

}
