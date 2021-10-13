#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

//---------------------functions----------------------------
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

double power(double x, int n) {
	double value = 1;
	for (int i = 0; i < n; i++) {
		value *= x;
	}
	return value;
}

int sign(double x) {
	if (x > 0) {
		return 1;
	}
	else if (x < 0) {
		return -1;
	}
	else return 0;
}

//--------------------------------------------


auto get_derivative(std::vector<double> coef_, int n) {//derivarive from array

	int size = coef_.size() - n;
	std::vector<double> array_(size);

	for (int i = 0; i < size; i++) {
		array_[i] = coef_[i] * double(factorial(coef_.size() - 1 - i)) / double(factorial(coef_.size() - 1 - i - n));
	}


	return array_;
}

double get_value(std::vector<double> coef_, double x) {//from array
	double value = 0;
	for (int i = 0; i < coef_.size(); i++) {
		value += power(x, coef_.size() - 1 - i) * coef_[i];
	}
	return value;
}

//------------------------------------------------

class Polynomial {

public:

	std::vector<double> coef;///coefs
	std::vector<std::vector<double>> derivatives;

	int N;//power of polynome (6)
	
	double a;//lower bound
	double b;//upper bound

	int negative;//number of negative roots
	int positive;//number of positive roots

	std::vector<std::vector<double>> intervals;
	std::vector<double> roots;

public:

	Polynomial(std::vector<double> coef_ ) {
		
		for (int i = 0; i < coef_.size(); i++) {
			derivatives.push_back(get_derivative(coef_, i));
		}
		

		N = coef_.size() - 1;
		for (int i = 0; i < N + 1; i++) {
			coef.push_back(coef_[i]);
		}

		this->get_a_b();
		//is->get_number_of_roots();
		//std::cout << "---------- polynomial constructed -------------" << std::endl;
	}

	//prints the polynomial
	void print() {

		std::cout << "f(x) = ";
		for (int i = 0; i < N + 1; i++) {
			std::cout << coef[i] << "x^" << N - i << " + ";
		}

		std::cout << 0 << std::endl;
	}


	void get_a_b() {// upper and lower bound od root circle


		double B;//max 0 - n-1
		double A;//max 1 - n


		B = coef[0];
		if (coef.size() > 1) {
			A = coef[1];
		}
		else {
			A = B;
			return;
		}

		for (int i = 0; i < N; i++) {
			if (coef[i] > B) {
				B = coef[i];
			}
			if (coef[i + 1] > A) {
				A = coef[i + 1];
			}

		}


		a = abs(coef[N]) / (abs(coef[N]) + B);
		b = 1 + A / abs(coef[0]);
	}

	int get_number_of_sign_changes(double x) {//return number of sign changes for given x

		int n = 0;
		std::vector<double> values;

		//std::cout << "\n" << x << " values:";
		for (int i = 0; i < N + 1; i++) {

			values.push_back(get_value(derivatives[i], x));
			//std::cout << i << ": " << values[i] << "       ";
		}

		//std::cout << std::endl;

		for (int i = 1; i < N + 1; i++) {
			if (sign(values[i]) * sign(values[i - 1]) < 0) {
				n++;
			}
		}

		return n;
	}
	
};


//passes to polynomial a number of positive and negative roots it has
void get_number_of_roots(Polynomial & p){

	
	p.negative = p.get_number_of_sign_changes(-p.b) - p.get_number_of_sign_changes(-p.a);
	p.positive = -p.get_number_of_sign_changes(p.b) + p.get_number_of_sign_changes(p.a);
}

//should return an array of intervals with one root only?

auto localizer_negative(Polynomial& p) {
	
	double lower = -p.b;
	double current = (-p.b - p.a) / 2;
	double upper = -p.a;

	int n_initial = p.get_number_of_sign_changes(-p.b);

	if (p.negative == 1) {
		p.intervals.push_back(std::vector<double>{lower, upper});
		return;
	}
	

	for (int i = 0; i < p.negative; i++) {
		while (p.get_number_of_sign_changes(current) != n_initial - 1) {
			current = (current + lower) / 2;
		}

		p.intervals.push_back(std::vector<double>{lower, current});
		lower = current;
		current = (upper + lower) / 2;
		n_initial--;
	}
	
	return;
};

auto localizer_positive(Polynomial& p) {///rewrite!!!
	

	double lower = p.a;
	double current = (p.b + p.a) / 2;
	double upper = p.b;
	std::cout << "lower : " << lower << " upper: " << upper << "\n";

	int n_initial = p.get_number_of_sign_changes(p.a);

	std::cout << "n_initial; " << n_initial << "\n";
	

	if (p.positive == 1) {
		p.intervals.push_back(std::vector<double>{lower, upper});
		return;
	}

	for (int i = 0; i < p.positive; i++) {

		while (p.get_number_of_sign_changes(current) != n_initial - 1) {
			std::cout << "in while : current =  " << current << "\n";
			current = (current + lower) / 2;
			

		}

		p.intervals.push_back(std::vector<double>{lower, current});
		lower = current;
		current = (upper + lower) / 2;
		n_initial--;
		std::cout << "after while " << i << "\n";
	}

	

	for (int i = 0; i < p.intervals.size(); i++) {
		std::cout << p.intervals[i][0] << ";" << p.intervals[i][1] << "\n";
	}

	return;
};


//newton method itself
auto solve(Polynomial& p, double e) {

	double x1;//x n
	double x2;//x n + 1

	for (int i = 0; i < p.positive; i++) {

		std::cout<< "sign changes: " << p.get_number_of_sign_changes(p.intervals[i][0]) << " " << p.get_number_of_sign_changes(p.intervals[i][1]) << "\n";

		if (get_value(p.coef, p.intervals[i][0]) * get_value(p.derivatives[2], p.intervals[i][0]) > 0) {
			x1 = p.intervals[i][0] + 2 * e;
			x2 = p.intervals[i][0];

		}
		else if (get_value(p.coef, p.intervals[i][1]) * get_value(p.derivatives[2], p.intervals[i][1]) > 0) {
			x1 = p.intervals[i][1] + 2 * e;
			x2 = p.intervals[i][1];

		}
		else {
			std::cout << "no suitable primary x \n";
			return;
			
		}

		x1 =( p.intervals[i][0] + p.intervals[i][1])/2 + 2 * e;
		x2 = (p.intervals[i][0] + p.intervals[i][1]) / 2;

		
		while (abs(x2 - x1) > e) {

			x1 = x2;			
			x2 = x1 - get_value(p.coef, x1) / get_value(p.derivatives[1], x1);/////		
			std::cout << "x1:  " << x1 <<"  x2:   " << x2 <<"\n";

		}

		p.roots.push_back(x2);
		std::cout << "root located: " << x2 << "	f(x): "<< get_value(p.derivatives[0], x2)<<"\n";
	}
}


int main() {

	
	std::vector<double> array_1 =  { 8.7 , 5.7, -2.8, 6.5, 6, -2.8, -0.4 };
	std::vector<double> array_2 = { 6, 5, 4, 3, 2, 1, 0 };
	std::vector<double> array_3 = { 1, 1, 1, 1, 1, 7, 0};
	std::vector<double> array_4 = { 1, 0, 0, 0, 0, 0, 1};
	std::vector<double> array_test = { 0.2, -1, -1, 3, -1, 5 };
	std::vector<double> array_mini = { 1, -5, 6 };//2, 3


	std::vector<double> array_my = { 0.0600438, -36.2684, 125.554, 48.7877, 5.81322, 0.0952521, -0.00619868};
	

	Polynomial p = Polynomial(array_my);//constructed polynomial with given coefs

	p.print();
	std::cout << "lower bound " << p.a << " upper bound " << p.b << std::endl;
	
	get_number_of_roots(p);

	std::cout << p.negative << " negative roots, " << p.positive << "positve roots " << std::endl;

	
	localizer_positive(p);
	
	solve(p, 0.00005);

	std::cout << "\n\n\n\n\ roots: \n";
	for (int i = 0; i < p.positive; i++) {
		std::cout << "root located: " << p.roots[i] << "	f(x): " << get_value(p.derivatives[0],p.roots[i]) << "\n";
	}


}


