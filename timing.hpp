/*
 * File:	timing.hpp
 * Author:	Stanley Goodwin
 * Provides a function that makes a string of an interval of time.
 */
#pragma once
#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <string>
#define hrc std::chrono::high_resolution_clock
#define time std::chrono::high_resolution_clock::time_point

std::string duration_string(time start, time stop);

#endif TIMING_H