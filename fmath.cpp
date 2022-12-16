/*
 * File:	fmath.cpp
 * Author:	Stanley Goodwin
 * Provides a set of fast math approximation functions.
 */
#include <cstdint>
#include "fmath.hpp"

/*
 * Calculates the inverse square root of a double-precision float.
 * @note	Inspired by the Quake 3 inverse square root algorithm.
 * @brief	Fast inverse square root.
 * @param	_number	double	The number to find the inverse square root of.
 * @return	invsqrt	double	The inverse square root of the number.
 */
double f_inv_sqrt(double _number)
{
	// Variables
	std::int64_t i;
	double x, y;
	const double threehalfs = 1.5;

	// Memory hacking (logarithm magic)
	x = _number * 0.5;
	y = _number;
	i = * ( std::int64_t * ) &y;        // Converts double memory address to 64-bit int
	i = 0x5fe6eb50c7b537a9 - (i >> 1);  // Scale & shift integer
	y = * ( double * ) &i;              // Converts 64-bit int memory address to double

	// Newton's method
	y = y * (threehalfs - (x * y * y));
	y = y * (threehalfs - (x * y * y));

	return y;  
}

/*
 * Calculates the square root of a double-precision float.
 * @note	Avoids doing division by doing multiplication instead.
 * @brief	Fast square root.
 * @param	_number	double	The number to find the square root of.
 * @return	sqrt	double	The square root of the number.
 */
double f_sqrt(double _number)
{
	return _number * f_inv_sqrt(_number);
}