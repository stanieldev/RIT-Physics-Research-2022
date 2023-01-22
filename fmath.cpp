/*
 * File:	fmath.cpp
 * Author:	Stanley Goodwin
 * Provides a set of fast math approximation functions.
 */
#include "fmath.hpp"

/*
 * Calculates the inverse square root of a double-precision float.
 * @note	Inspired by the Quake 3 inverse square root algorithm.
 */
double f_inv_sqrt(double _number)
{
	long long i;
	double x, y;
	const double threehalfs = 1.5;

	x = _number * 0.5;
	y = _number;
	i = * ( long long * ) &y;           // Converts double memory address to 64-bit int
	i = 0x5fe6eb50c7b537a9 - (i >> 1);  // Scale & shift integer
	y = * ( double * ) &i;              // Converts 64-bit int memory address to double

	y = y * (threehalfs - (x * y * y));
	y = y * (threehalfs - (x * y * y));

	return y;
}

double f_sqrt(double _number)
{
	return _number * f_inv_sqrt(_number);
}