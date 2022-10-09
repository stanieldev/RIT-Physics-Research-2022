/**
 * File:	fmath.cpp
 * Brief:	Provides fast math operations.
 *
 * Author:	Stanley Goodwin
 * Contact:	sfg99709akwork@gmail.com
 *
 * Creation Date: 6/21/2022
 * Last Modified: 7/31/2022
 */
#include <cstdint>


// Status: COMPLETE & VERIFIED
/*
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

// Status: COMPLETE & VERIFIED
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












#include <iostream>
// Find a faster determinant or gaussian elimination
// STATUS: INCOMPLETE
/**
 * Creates a (n-1)x(n-1) sub-matrix of an input matrix.
 *
 * @brief	Submatrix generator.
 * @param	matrix	Any matrix of sides <= 6 to find the submatrix of.
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
double determinant(double matrix[6][6], int n)  // O(n!)
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
		// Generate the submatrix
		submatrix(matrix, _submatrix, 0, i, n);

		// Add the component to the total determinant
		if (matrix[0][i] == 0) { sign = -sign; continue; }
		else { det += sign * matrix[0][i] * determinant(_submatrix, n - 1); }

		// Swap the sign
		sign = -sign;
	}

	// Return the determinant
	return det;
}  


void print_array(double matrix[6][6])
{
	std::printf("\n[");
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			std::printf("%f, ", matrix[i][j]);
		}
		if (i != 5) { std::printf("\n"); }
		else { std::printf("]\n"); }
	}
}

void three_matrix_multiply(double U[6][6], double matrix[6][6], double V[6][6], int n)
{
	double output_matrix[6][6] = {
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0}
	};

	for (int i = n; i < 6; i++)
	{
		for (int j = n; j < 6; j++)
		{
			output_matrix[i][j] = 0;
			for (int k = n; k < 6; k++)
			{
				output_matrix[i][j] += U[i][k] * matrix[k][j];
			}
		}
	}

	for (int i = n; i < 6; i++)
	{
		for (int j = n; j < 6; j++)
		{
			matrix[i][j] = 0;
			for (int k = n; k < 6; k++)
			{
				matrix[i][j] += output_matrix[i][k] * V[k][j];
			}
		}
	}
}


double invdeterminant_by_matrix(double matrix[6][6], int n)
{
	int offset = 6 - n;
	double a = -matrix[offset][offset];

	for (int j = 0; j < 6; j++)
	{
		if (matrix[j][j] == 0) {
			return determinant(matrix, 6);
		}
	}

	if (n == 1) { return matrix[5][5]; }

	double U[6][6] = {
		{a, 0, 0, 0, 0, 0},
		{0, a, 0, 0, 0, 0},
		{0, 0, a, 0, 0, 0},
		{0, 0, 0, a, 0, 0},
		{0, 0, 0, 0, a, 0},
		{0, 0, 0, 0, 0, a}
	};

	double V[6][6] = {
		{a, 0, 0, 0, 0, 0},
		{0, a, 0, 0, 0, 0},
		{0, 0, a, 0, 0, 0},
		{0, 0, 0, a, 0, 0},
		{0, 0, 0, 0, a, 0},
		{0, 0, 0, 0, 0, a}
	};

	for (int i = 7 - n; i < 6; i++)
	{
		U[i][offset] = matrix[i][offset];
		V[offset][i] = matrix[offset][i];
	}
	//print_array(U);
	//print_array(matrix);
	//print_array(V);

	three_matrix_multiply(U, matrix, V, n - 1);

	return -a * invdeterminant_by_matrix(matrix, n - 1);
}

void print_array(double matrix[6][7])
{
	std::printf("\n[");
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			std::printf("%f, ", matrix[i][j]);
		}
		if (i != 5) { std::printf("\n"); }
		else { std::printf("]\n"); }
	}
}


// MULTIPLICATION BY 0 IN DIAGONAL
double determinant6(double matrix[6][6])  // O(n^3)
{
	double scalar;
	double det = 1.0;

	//print_array(matrix);

	for (int i = 0; i < 6; i++)
	{
		if (matrix[i][i] == 0)
		{
			return determinant(matrix, 6);
		}
	}

	for (int col = 0; col < 6; col++)
	{
		for (int row = col + 1; row < 6; row++)
		{
			if (matrix[row][col] == 0)
			{
				continue;
			}
			
			scalar = matrix[col][col] / matrix[row][col];
			det /= scalar;

			for (int i = col; i < 6; i++)
			{
				matrix[row][i] = matrix[col][i] - scalar * matrix[row][i];
			}
			//print_array(matrix);
		}
	}

	for (int i = 0; i < 6; i++)
	{
		det *= matrix[i][i];
	}

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
 * In the future, use: https://en.wikipedia.org/wiki/LAPACK
 * https://www.sciencedirect.com/topics/mathematics/partial-pivoting#:~:text=The%20partial%20pivoting%20technique%20is,to%20its%20remaining%20row%20entries.
 *
 * @brief	6x6 Matrix Cramer's rule.
 * @param	matrix	The 6x6 matrix.
 * @param	output	The output vector.
 * @param	coeff	The coefficients vector.
 */
void cramer(double matrix[6][6], double output[6], double coeff[6])
{
	// Variables
	double det = determinant6(matrix);
	double r_det = 0;  // Component determinant
	double r_matrix[6][6];

	// Cramer's rule
	for (int k = 0; k < 6; k++)
	{
		// Swap column k of matrix with output array
		swap_column(matrix, output, k, r_matrix);

		// Calculate the determinant of the new matrix
		//r_det = determinant(r_matrix, 6);
		r_det = determinant6(r_matrix);

		// Place ration into the coefficient array
		coeff[k] = r_det / det;
	}
}


// https://www.cfm.brown.edu/people/dobrush/am34/MuPad/LU.html#:~:text=Gaussian%20elimination%20and%20Gauss%2D%2D,applied%20to%20any%20vector%20b.
void gaussian_elimination(double augmented_matrix[6][7], double coeff[6])
{
	double scalar;

	print_array(augmented_matrix);
	for (int diag = 0; diag < 6; diag++)
	{
		for (int row = diag + 1; row < 6; row++)
		{
			if (augmented_matrix[row][diag] == 0)
			{
				continue;
			}

			scalar = augmented_matrix[diag][diag] / augmented_matrix[row][diag];
			for (int col = diag; col < 7; col++)
			{
				augmented_matrix[row][col] = augmented_matrix[diag][col] - scalar * augmented_matrix[row][col];
			}
			print_array(augmented_matrix);
		}
	}

	for (int i = 0; i < 6; i++)
	{
		coeff[i] = augmented_matrix[i][6];
	}
}
