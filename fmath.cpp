/**
 * @file	fmath.cpp
 * @brief	Provides fast math operations.
 * 
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 * 
 * Creation Date: 6/21/2022
 * Last Modified: 7/5/2022
 */
#include <cstdint>


/**
 * Calculates the inverse square root of a double-precision float.
 * https://en.wikipedia.org/wiki/Fast_inverse_square_root                       (x32)
 * https://www.youtube.com/watch?v=p8u_k2LIZyo                                  (x32)
 * https://stackoverflow.com/questions/11644441/fast-inverse-square-root-on-x64 (x64)
 * 
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

	// Bit-hacking (logarithm magic)
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

/**
 * Calculates the square root of a double-precision float.
 *
 * @brief	Fast square root.
 * @param	number	double	The number to find the square root of.
 * @return	sqrt	double	The square root of the number.
 */
double f_sqrt(double number)
{
	// Avoiding division with multiplication
	return number * f_inv_sqrt(number);
}
