#include "math_functions.h"


// a list of some of the auxiliary functions we need for polymer generation and polymer class.
//includes:
//vector functions
//factorial and choose

//to optimize these functions in the future it will probably be useful to pass all arguments 
//by reference and/or use pointers.
std::vector<double> vector_addition(std::vector<double> vector1, std::vector<double> vector2) {

    return { vector1[0] + vector2[0],vector1[1] + vector2[1],vector1[2] + vector2[2] };
}
std::vector<double> vector_subtraction(std::vector<double> vector1, std::vector<double> vector2) {
    return { vector1[0] - vector2[0],vector1[1] - vector2[1],vector1[2] - vector2[2] };

}

double vector_modulus(std::vector<double> vector1) {
    return sqrt(pow(abs(vector1[0]), 2) + pow(abs(vector1[1]), 2) + pow(abs(vector1[2]), 2));
}

double dot_product(std::vector<double> vector1, std::vector<double> vector2) {
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}


std::vector<double> normalize(std::vector<double> &unnormalized_vector) {
    std::vector<double> unit_vector;
    for (int i{ 0 }; i < unnormalized_vector.size(); i++) {
        unit_vector.push_back(unnormalized_vector[i] / vector_modulus(unnormalized_vector));
    }
    return unit_vector;

}
std::vector<double> multiplication_by_scalar(double constant, std::vector<double> vector) {
    return { constant * vector[0],constant * vector[1],constant * vector[2] };
}
double dist_2_points3d(std::vector<double> v1, std::vector<double> v2)
{
    return sqrt(pow(v1[0]-v2[0],2)+ pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2));
}
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

//factorial function returns a double to maximise 
double factorial(int n) {
    if (n == 0) return 1;
    else if (n < 0) {
        std::cout << "Negative integer factorial called." << std::endl;
        exit(0);
    }
    else {
        double summand{ 1 };
        for(int i = 1; i <= n; ++i)
            summand *= i;
        return summand;
    }
}
double choose(int n, int  k) {
    if (k == 0) return 1;
    return (n * choose(n - 1, k - 1)) / k;

}
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//
//**************************************************************************************////**************************************************************************************//

