#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <random>
#include "math.h"

#define EPS 0.00001
using namespace std;
typedef double(*pointFunction)(double);

default_random_engine generator(time(0));
uniform_real_distribution <double> distribution(0, 1);

vector<vector<double>> T = { { 0 },
{ -0.5773502692, 0.5773502692 },
{ -0.7745966692, 0.0000000000, 0.7745966692 },
{ -0.8611363115, -0.3399810436, 0.3399810436, 0.8611363115 },
{ -0.9061798459, -0.5384693101, 0.0000000000, 0.5384693101, 0.9061798459 },
{ -0.9324695142, -0.6612093864, -0.2386191861, 0.2386191861, 0.6612093864, 0.9324695142 } };

vector<vector<double>> A = { { 2 },
{ 1, 1 },
{ 0.5555555556,0.8888888888,0.5555555556 },
{ 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 },
{ 0.2369268851, 0.4786286705, 0.5688888888,  0.4786286705,0.2369268851 },
{ 0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924 } };

double inputFunction( double x) {
	//return (2 * pow(x, 3) - 7 * x + 4);
	//return (2 * pow(x, 3) / pow(x, 4));
	//return (pow(x, 2) + x);
	//return sin(x);
	//return 1 / x;
	//return x / (pow(x, 4) + 4);
	return pow(x, 2);
}
long double leftRectangleIntegral(double a, double b, int n, pointFunction f) {
	double h = (b - a) / n;
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += f(a + i*h);
	}
	long double rectangle_integral_result = h*sum;
	return rectangle_integral_result;
}

long double trapezeIntegral(double a, double b, int n, pointFunction f) {
	double h = (b - a) /(double)n;
	double sum = 0.5*(f(a) + f(b));
	for (int i = 1; i < n; i++) {
		double x_i = a + i*h;
		sum +=f(x_i);
		//if (!isnormal(f(x_i))) {
		//	cout << "*****************************************" << endl;
		//};
	}
	long double trapeze_integral_result = sum*h;
	return trapeze_integral_result;
}

long double simpsonIntegral(double a, double b, int n, pointFunction f) {
	double h, sum2 = 0, sum4 = 0, sum = 0;
	h = (b - a) / ( 2*n);//Шаг интегрирования.
	for (int i = 1; i < 2*n; i += 2) {
		sum4 += f(a + h*i);       //Значения с нечётными индексами, которые нужно умножить на 4.
		sum2 += f(a + h*(i + 1)); //Значения с чётными индексами, которые нужно умножить на 2.
	}
	sum = f(a) + 4 * sum4 + 2 * sum2 - f(b);//Отнимаем значение f(b) так как ранее прибавили его дважды. 
	long double simpsom_intrgral_result = (h / 3)*sum;
	return simpsom_intrgral_result;

}

long double gaussIntegral(double a, double b, int n, pointFunction f) {
	//return (b - a) / 2 * (f((a + b) / 2 - (b - a) / (2 * sqrt(3))) + f((a + b) / 2 + (b - a) / (2 * sqrt(3)))); для 2 узлов
	double h = (b - a) / 2;
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += A[n-1][i] * f((a + b) / 2 + ((b - a) / 2)*T[n-1][i]);
	}
	long double gauss_integtal_result = h*sum;
	return gauss_integtal_result;
}

long double monteCarloIntegral(double a, double b, int n, pointFunction f) {
	double h = (b - a) / n;
	double sum = 0, tmp;
	for (int i = 1; i < n; i++) {
		tmp = distribution(generator);
		sum += f(a + tmp*(b - a));
	}
	long double monte_carlo_integral_result = h*sum;
	return monte_carlo_integral_result;
}

int main() {
	start:
	//cout << "Input function: 2x^3-7x+4" << endl;
	//cout << "Input function: 1/x" << endl;

	//cout << "Input function: x / (pow(x, 4) + 4)" << endl;
	cout << "Input function: pow(x,2)" << endl;
	//cout << "If a=0, b=2 --> result = 2" << endl;
	//cout << "If a=0.5, b=3 --> result = 3.584" << endl;

	long double result;
	cout << "Expected result: "; 
	cin >> result;
	double a, b, nR, nT, n, nS, nMC;
	cout << "Enter a: ";
	cin >> a;
	cout << "Enter b: ";
	cin	>> b;
	
	//cout << "***********" << inputFunction(a);
	//cout << inf

	//if (!isnormal(inputFunction(a))|| !isnormal(inputFunction(b))) {
	//	cout << "You input wrong limits of integration, try again!"<<endl;
	//	goto start;
		//system("pause");
		//return 0;
	//}








	cout << "Enter number of splits for Rectangle Method:     ";
	cin >> nR;
	cout << "Enter number of splits for Trapeze Method:       ";
	cin >> nT;		
	cout << "Enter number of splits for Simpson Method:       ";
	cin >> nS;
	cout << "Enter number of nodes for Gauss Method:          ";
	cin >> n;
	cout << "Enter number of splits for Monte-Carlo Method:   ";
	cin >> nMC;

	long double rectangle_result = leftRectangleIntegral(a, b, nR, inputFunction);
	
	if (!isnormal(rectangle_result)) {
		cout << "You function has not any results, try again!" << endl;
		goto start;
	}

	cout << "Rectangle integral result:  " << rectangle_result << " Inaccuracy of measurements: " << abs(rectangle_result - result) << endl;
	long double trapeze_result = trapezeIntegral(a, b, nT, inputFunction);
	cout << "Trapeze integral method:    " << trapeze_result << " Inaccuracy of measurements: " << abs(trapeze_result - result) << endl;
	long double simpson_result = simpsonIntegral(a, b, nS, inputFunction);
	cout << "Simpson integral method:    " << simpson_result << " Inaccuracy of measurements: " << abs(simpson_result - result) << endl;
	long double gauss_result = gaussIntegral(a, b, n, inputFunction);
	cout << "Gauss integral method:      " << gauss_result << " Inaccuracy of measurements: " << abs(gauss_result - result) << endl;
	long double monte_carlo_result = monteCarloIntegral(a, b, nMC, inputFunction);
	cout << "MonteCarlo integral method: " << monte_carlo_result << " Inaccuracy of measurements: " << abs(monte_carlo_result - result) << endl;

	system("pause");
	return 0;


}