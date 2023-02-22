#pragma once
#include <iostream>
#include<cmath>
#include<vector>
#include<string>
#include <time.h>

//VECTOR FUNCTIONS
std::vector<double> vector_addition(std::vector<double> vector1, std::vector<double> vector2);
std::vector<double> vector_subtraction(std::vector<double> vector1, std::vector<double> vector2);
double vector_modulus(std::vector<double> vector1);
double dot_product(std::vector<double> vector1, std::vector<double> vector2);
std::vector<double> normalize(std::vector<double> &unnormalized_vector);
std::vector<double> multiplication_by_scalar(double constant, std::vector<double> vector);
double dist_2_points3d(std::vector<double> v1, std::vector<double> v2);
//OTHER
double factorial(int n);
double choose(int n, int  k);
