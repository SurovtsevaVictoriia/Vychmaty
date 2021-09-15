#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

int N0 = 6;

int factorial(int n) {
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


class Polynomial {

public:

	std::vector<double> coef;///coefs

	int N;

	double B;//max 0 - n-1
	double A;//max 1 - n

	double a;//lower bound
	double b;//upper bound

	int negative;
	int positive;


public:

	Polynomial(std::vector<double> coef_ ) {
		
		N = coef_.size() - 1;
		for (int i = 0; i < N + 1; i++) {
			coef.push_back(coef_[i]);
		}

		this->get_A_B();
		this->get_a_b();

		std::cout << "---------- polynomial constructed -------------" << std::endl;
	}

	double get_value(double x) {
		double value = 0;
		for (int i = 0; i < N + 1; i++) {
			value += power(x, N - i) * coef[i];
		}
		return value;
	}

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

	void get_number_of_roots() {
		Row r = Row(*this);

		int negative = r.get_number(-p.b) - r.get_number(-p.a);
		int positive = -r.get_number(p.b) + r.get_number(p.a);
	}
};

auto get_derivative(Polynomial& poly, int n) {//returns Polynomial which is 

	std::vector<double> array_(poly.N + 1 - n);

	for (int i = 0; i < poly.N + 1 - n; i++) {

		array_[i] = poly.coef[i] * double(factorial(poly.N - i)) / double(factorial(poly.N - i - n));
	}

	Polynomial deriv = Polynomial(array_);
	return  deriv;

}

class Row {//Budan - Fourier
	// хранит N полиномов и выдвет их значения по изначальному полиному

public:

	Polynomial* poly;

	std::vector< Polynomial > functions;
	
	int N;

	Row(Polynomial& poly_) {

		poly = &poly_;
		N = poly->N  + 1;
		for (int i = 0; i < poly->N +1; i++) {

			functions.push_back(get_derivative(*poly, i));
			std::cout << i << std::endl;
			functions[i].print();
			
		}

	}

	int get_number(double x) {//return number of sign changes for given x

		int n = 0;
		std::vector<double> values;

		std::cout << "\n" << x << " values:";
		for (int i = 0; i < N ; i++) {

			values.push_back(functions[i].get_value(x));			
			std::cout <<i <<  ": " << values[i] << "       ";
		}

		std::cout << std::endl;

		for (int i = 1; i < N; i++) {
			if (sign(values[i] )* sign(values[i - 1] )< 0) {
				n++;
			}
		}

		return n;
	}


	void rough_localize() {

	}
};



int main() {


	std::cout << "6! = " << factorial(6) << std::endl;
	std::cout << "2^3 = " << power(2, 3) << std::endl;


	std::vector<double> array_1 =  { 8.7 , 5.7, -2.8, 6.5, 6, -2.8, -0.4 };
	std::vector<double> array_2 = { 6, 5, 4, 3, 2, 1, 0 };
	std::vector<double> array_3 = { 1, 1, 1, 1, 1, 7, 0};
	std::vector<double> array_test = { 0.2, -1, -1, 3, -1, 5 };
	

	Polynomial p = Polynomial(array_test);//constructed polynomial with given coefs

	p.print();
	std::cout << "lower bound " << p.a << " upper bound " << p.b << std::endl;
	
	Row r = Row(p);//this is a row from Budan-Fourier theorem

	int negative = r.get_number(-p.b) - r.get_number(-p.a);
	int positive = -r.get_number(p.b) + r.get_number(p.a);

	
	
	std::cout << "---- shit ----" << std::endl;
		for (int i = -100; i < 100; i+=5) {

			if (i == 0) {
				std::cout << "!!!!!!";
			}

			std::cout << r.get_number(double(i)/100) << " ";
		}
	std::cout << std::endl << "---- shit ----" << std::endl;
	
	
	
	std::cout << r.get_number(-p.b) << " " << r.get_number(-p.a) << " " << r.get_number(p.a) << " " <<  r.get_number(p.b) << std::endl;
	/*std::cout << r.get_number(-2.5) << " " << r.get_number(-0.625) << " " << r.get_number(1.5) << " " <<  r.get_number(2) << 
		" "<< r.get_number(5.4) << " " << r.get_number(5.5) << std::endl;*/
	std::cout << negative << " negative roots, " << positive << "positve roots " << std::endl;

	//localizing roots
	// we know the\at we have "negative" roots on (-a -b) and "positive" roots on (a b)
	// now we have to separate thems


	



}


