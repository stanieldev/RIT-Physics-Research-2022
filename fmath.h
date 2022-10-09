//
// File: fmath.h
// Author: Stanley Goodwin
// Creation Date: 6/21/2022
// Last Modified: 7/4/2022
//
#pragma once
#ifndef FMATH_H
#define FMATH_H


// A faster form of inverse square roots (for normalization)
double f_inv_sqrt(double number);

// A faster form of square rooting using inverse square root
double f_sqrt(double number);


#endif FMATH_H
