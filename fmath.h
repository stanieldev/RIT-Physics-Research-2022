/**
 * File:	fmath.h
 * Brief:	Provides fast math operations.
 *
 * Author:	Stanley Goodwin
 * Contact:	sfg99709akwork@gmail.com
 *
 * Creation Date: 6/21/2022
 * Last Modified: 7/31/2022
 */
#pragma once
#ifndef FMATH_H
#define FMATH_H


// Fast square roots
double f_inv_sqrt(double number);  // Fast inverse square root
double f_sqrt(double number);      // Fast square root

// Using cramers rule to find inverse matrix values
void cramer(double matrix[6][6], double output[6], double coeff[6]);
void gaussian_elimination(double augmented_matrix[6][7], double coeff[6]);


#endif FMATH_H