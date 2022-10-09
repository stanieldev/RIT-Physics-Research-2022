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
#include <iostream>
#include <vector>


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


/**
 * Creates a (n-1)x(n-1) sub-matrix of an input matrix.
 *
 * @brief	Submatrix generator.
 * @param	matrix	Any matrix of sides < 7 to find the submatrix of.
 * @param	output	The submatrix output.
 * @param	a	The row to have removed.
 * @param	b	The column to have removed.
 * @param	n	The dimensions of @param matrix.
 */
void submatrix(double matrix[6][6], double output[6][6], int a, int b, int n)
{
	// Indeces
	int i = 0;
	int j = 0;

	// Filling the submatrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			if (row != a && col != b)
			{
				output[i][j++] = matrix[row][col];
				if (j == n - 1) { j = 0; i++; }
			}
		}
	}
}

/**
 * Finds the determinant of any 6x6 or smaller square matrix.
 *
 * @brief	Determinant of a matrix.
 * @param	matrix	Any matrix of sides < 7 to find the determinant of.
 * @param	n	The dimensions of @param matrix.
 */
double determinant(double matrix[6][6], int n)
{
	// Variables
	int sign = 1;
	double det = 0.0;
	double _submatrix[6][6];

	// 2x2 matrix determinant
	if (n == 2) { return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]; }

	// Larger than 2x2 determinant
	for (int i = 0; i < n; i++)
	{
		submatrix(matrix, _submatrix, 0, i, n);  // Generate the submatrix
		det += sign * matrix[0][i] * determinant(_submatrix, n - 1);
		sign = -sign;
	}

	// Return the determinant
	return det;
}


/**
 * Swaps a column in a 6x6 matrix with a chosen vector.
 *
 * @brief	Swap a column.
 * @param	matrix	The 6x6 matrix to have replaced.
 * @param	output	The new column.
 * @param	column	The column to swap.
 * @param	matrix	The new column-swapped matrix.
 */
void swap_column(double matrix[6][6], double output[6], int column, double new_matrix[6][6])
{
	// Filling the sub matrix
	for (int row = 0; row < 6; row++)
	{
		for (int col = 0; col < 6; col++)
		{
			if (col == column)
				new_matrix[row][col] = output[row];
			else
				new_matrix[row][col] = matrix[row][col];
		}
	}
}

/**
 * Returns a list of 6 coefficients from a matrix and output vector.
 * https://en.wikipedia.org/wiki/LAPACK
 * 
 * @brief	6x6 Matrix Cramer's rule.
 * @param	matrix	The 6x6 matrix.
 * @param	output	The output vector.
 * @param	coeff	The coefficients vector.
 */
void cramer(double matrix[6][6], double output[6], double coeff[6])
{
	double det = determinant(matrix, 6);
	printf("\n");
	double rdet = 0;
	double new_matrix[6][6];

	for (int k = 0; k < 6; k++) {
		swap_column(matrix, output, k, new_matrix);
		rdet = determinant(new_matrix, 6);
		coeff[k] = rdet / det;
	}
}
