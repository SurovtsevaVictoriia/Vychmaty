#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

class Row;
class Polynomial;

int N0 = 6;

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

	double B;//max 0 - n-1
	double A;//max 1 - n

	double a;//lower bound
	double b;//upper bound

	int negative;//number of negative roots
	int positive;//number of positive roots

	std::vector<std::vector<double>> intervals;

public:

	Polynomial(std::vector<double> coef_ ) {
		
		for (int i = 0; i < coef_.size(); i++) {
			derivatives.push_back(get_derivative(coef_, i));
		}
		

		N = coef_.size() - 1;
		for (int i = 0; i < N + 1; i++) {
			coef.push_back(coef_[i]);
		}

		this->get_A_B();
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
		//--------
	
	}
		
	void get_A_B() {//some max values of coefs

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
	}

	void get_a_b() {// upper and lower bound od root circle
		a = abs(coef[N]) / (abs(coef[N]) + B);
		b = 1 + A / abs(coef[0]);
	}

	int get_number_of_sign_changes(double x) {//return number of sign changes for given x

		int n = 0;
		std::vector<double> values;

		std::cout << "\n" << x << " values:";
		for (int i = 0; i < N + 1; i++) {

			values.push_back(get_value(derivatives[i], x));
			std::cout << i << ": " << values[i] << "       ";
		}

		std::cout << std::endl;

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
			current = (current - lower) / 2;
		}

		p.intervals.push_back(std::vector<double>{lower, current});
		lower = current;
		current = (upper - lower) / 2;
		n_initial--;
	}
	
	return;
};

auto localizer_positive(Polynomial& p) {
	
	double lower = p.a;
	double current = (p.b + p.a) / 2;
	double upper = p.b;

	int n_initial = p.get_number_of_sign_changes(p.a);

	
	if (p.positive == 1) {
		p.intervals.push_back(std::vector<double>{lower, upper});
		return;
	}

	for (int i = 0; i < p.positive; i++) {
		while (p.get_number_of_sign_changes(current) != n_initial - 1) {
			current = (current - lower) / 2;
		}

		p.intervals.push_back(std::vector<double>{lower, current});
		lower = current;
		current = (upper - lower) / 2;
		n_initial--;
	}


	for (int i = 0; i < p.intervals.size(); i++) {
		std::cout << p.intervals[i][0] << ";" << p.intervals[i][1] << "\n";
	}

	return;
};


class Newton {
	
	// 

	//newton method itself
	auto solve(Polynomial p, double a, double b, double e) {
		
		

		double x0;//x n-1
		double x1;//x n
		double x2;//x n + 1

		while (abs(x2 - x1) > e) {
			x2 = x1;/////
			x0 = x1;
			x1 = x2;
		}

		return x2;
	}
};

int main() {

	/*std::cout << "6! = " << factorial(6) << std::endl;
	std::cout << "2^3 = " << power(2, 3) << std::endl;*/


	std::vector<double> array_1 =  { 8.7 , 5.7, -2.8, 6.5, 6, -2.8, -0.4 };
	std::vector<double> array_2 = { 6, 5, 4, 3, 2, 1, 0 };
	std::vector<double> array_3 = { 1, 1, 1, 1, 1, 7, 0};
	std::vector<double> array_4 = { 1, 0, 0, 0, 0, 0, 1};
	std::vector<double> array_test = { 0.2, -1, -1, 3, -1, 5 };


	std::vector<double> array_my = { 0.0600438, -36.2684, 125.554, 48.7877, 5.81322, 0.0952521, -0.00619868};
	

	Polynomial p = Polynomial(array_my);//constructed polynomial with given coefs

	p.print();
	std::cout << "lower bound " << p.a << " upper bound " << p.b << std::endl;
	
	get_number_of_roots(p);

	std::cout << p.negative << " negative roots, " << p.positive << "positve roots " << std::endl;

	
	//localizing roots
	// we know the\at we have "negative" roots on (-a -b) and "positive" roots on (a b)
	// now we have to separate thems

	//localizer_negative(p);


	localizer_positive(p);
	//this shit is done
	//now for precision


	/*std::cout << "---- shit ----" << std::endl;
	for (int i = -100; i < 100; i += 5) {

		if (i == 0) {
			std::cout << "!!!!!!";
		}

		std::cout << r.get_number_of_sign_changes(double(i) / 100) << " ";
	}
	std::cout << std::endl << "---- shit ----" << std::endl;*/

}


