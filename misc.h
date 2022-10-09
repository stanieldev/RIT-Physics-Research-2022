/**
 * @file	misc.h
 * @brief	Random print functions for timing things.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/23/2022
 * Last Modified: 7/5/2022
 */
#pragma once
#ifndef MISC_H
#define MISC_H

#include <chrono>
#define hrc std::chrono::high_resolution_clock
#define time std::chrono::high_resolution_clock::time_point


void print_progress(double percentage);  // Progress bar function (print to console)
std::string _duration_string(time start, time stop);    // Duration function (print to console)


#endif MISC_H