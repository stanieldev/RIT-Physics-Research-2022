/**
 * @file	timing.h
 * @brief	Random print functions for timing things.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/23/2022
 * Last Modified: 7/31/2022
 */
#pragma once
#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <string>
#define hrc std::chrono::high_resolution_clock
#define time std::chrono::high_resolution_clock::time_point


std::string _duration_string(time start, time stop);  // Return duration as a simple string


#endif TIMING_H