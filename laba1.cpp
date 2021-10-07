#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

int N0 = 6;

//---------------------functions----------------------------
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

//--------------------------------------------



class Polynomial {

public:

	std::vector<double> coef;///coefs

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
		
		N = coef_.size() - 1;
		for (int i = 0; i < N + 1; i++) {
			coef.push_back(coef_[i]);
		}

		this->get_A_B();
		this->get_a_b();
		//is->get_number_of_roots();
		//std::cout << "---------- polynomial constructed -------------" << std::endl;
	}

	//returns the value of polynomial at given x
	double get_value(double x) {
		double value = 0;
		for (int i = 0; i < N + 1; i++) {
			value += power(x, N - i) * coef[i];
		}
		return value;
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

	
};

auto get_derivative(Polynomial& poly, int n) {//returns Polynomial which is a derivativ of n order

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
			//std::cout << i << std::endl;
			//functions[i].print();
			
		}

	}

	int get_number_of_sign_changes(double x) {//return number of sign changes for given x

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

};



//passes to polynomial a number of positive and negative roots it has
void get_number_of_roots(Polynomial & p){

	Row r = Row(p);
	p.negative = r.get_number_of_sign_changes(-p.b) - r.get_number_of_sign_changes(-p.a);
	p.positive = -r.get_number_of_sign_changes(p.b) + r.get_number_of_sign_changes(p.a);
}

//should return an array of intervals with one root only?


auto localizer_negative(Polynomial& p) {
	//nuber : poly_.negative
	//-p.b <> -p.a 
	// делим пополам и смотрим что получилось?

	


	Row r = Row(p);

	double lower = -p.b;
	double current = (-p.b - p.a) / 2;
	double upper = -p.a;

	int n_initial = r.get_number_of_sign_changes(-p.b);

	if (p.negative == 1) {
		p.intervals.push_back(std::vector<double>{lower, upper});
		return;
	}
	

	for (int i = 0; i < p.negative; i++) {
		while (r.get_number_of_sign_changes(current) != n_initial - 1) {
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
	
	


	Row r = Row(p);

	double lower = p.a;
	double current = (p.b + p.a) / 2;
	double upper = p.b;

	int n_initial = r.get_number_of_sign_changes(p.a);

	
	if (p.positive == 1) {
		p.intervals.push_back(std::vector<double>{lower, upper});
		return;
	}

	for (int i = 0; i < p.positive; i++) {
		while (r.get_number_of_sign_changes(current) != n_initial - 1) {
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



int main() {


	/*std::cout << "6! = " << factorial(6) << std::endl;
	std::cout << "2^3 = " << power(2, 3) << std::endl;*/


	std::vector<double> array_1 =  { 8.7 , 5.7, -2.8, 6.5, 6, -2.8, -0.4 };
	std::vector<double> array_2 = { 6, 5, 4, 3, 2, 1, 0 };
	std::vector<double> array_3 = { 1, 1, 1, 1, 1, 7, 0};
	std::vector<double> array_4 = { 1, 0, 0, 0, 0, 0, 1};
	std::vector<double> array_test = { 0.2, -1, -1, 3, -1, 5 };
	

	Polynomial p = Polynomial(array_test);//constructed polynomial with given coefs

	p.print();
	std::cout << "lower bound " << p.a << " upper bound " << p.b << std::endl;
	
	get_number_of_roots(p);

	std::cout << p.negative << " negative roots, " << p.positive << "positve roots " << std::endl;

	
	//localizing roots
	// we know the\at we have "negative" roots on (-a -b) and "positive" roots on (a b)
	// now we have to separate thems

	localizer_negative(p);
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


