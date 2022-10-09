//
// File: fmath.cpp
// Author: Stanley Goodwin
// Creation Date: 6/21/2022
// Last Modified: 7/4/2022
//
#include <cstdint>


// Fast Inverse Square Root
// https://en.wikipedia.org/wiki/Fast_inverse_square_root (x32)
// https://stackoverflow.com/questions/11644441/fast-inverse-square-root-on-x64 (x64)
// https://www.youtube.com/watch?v=p8u_k2LIZyo (explaination of x32)
// Replicated to have same structure as x32 for the x64 double architecture vs a float number
// Returns the inverse square root of a number
double f_inv_sqrt(double number)
{
	std::int64_t i;
	double x2, y;
	const double threehalfs = 1.5;

	x2 = number * 0.5;
	y  = number;
	i  = * ( std::int64_t * ) &y;
	i  = 0x5fe6eb50c7b537a9 - (i >> 1);
	y  = * ( double * ) &i;
	y  = y * (threehalfs - (x2 * y * y));
	y  = y * (threehalfs - (x2 * y * y));

	return y;
}


// Fast Square Root
// Returns the square root of a number
double f_sqrt(double number)
{
	return number * f_inv_sqrt(number);
}