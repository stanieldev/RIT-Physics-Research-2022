/*
 * File:	fmath.cpp
 * Author:	Stanley Goodwin
 * Provides a set of fast math approximation functions.
 */
#include <cstdint>

/*
 * Calculates the inverse square root of a double-precision float.
 * Inspired by the Quake 3 inverse square root algorithm.
 * @brief	Fast inverse square root.
 * @param	number	double	The number to find the inverse square root of.
 * @return	invsqrt	double	The inverse square root of the number.
 */
double f_inv_sqrt(double number)
{
	// Initialization
	std::int64_t i;
	double x2, y;
	const double threehalfs = 1.5;

	// Memory hacking (logarithm magic)
	x2 = number * 0.5;
	y  = number;
	i  = * ( std::int64_t * ) &y;        // Convert double memory address to 64bit int
	i  = 0x5fe6eb50c7b537a9 - (i >> 1);  // Scale & shift integer
	y  = * ( double * ) &i;              // Convert 64bit int memory address to double

	// Newton iteration(s)
	y  = y * (threehalfs - (x2 * y * y));
	y  = y * (threehalfs - (x2 * y * y));

	return y;
}

/*
 * Calculates the square root of a double-precision float.
 * Avoided doing division by multiplcation.
 * @brief	Fast square root.
 * @param	number	double	The number to find the square root of.
 * @return	sqrt	double	The square root of the number.
 */
double f_sqrt(double number)
{
	return number * f_inv_sqrt(number);
}