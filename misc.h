//
// File: misc.h
// Author: Stanley Goodwin
// Creation Date: 6/23/2022
// Last Modified: 7/4/2022
//
#pragma once
#ifndef MISC_H
#define MISC_H

#include <chrono>

#define start_timer std::chrono::high_resolution_clock::now()


// Prints the current progress % to console
void print_progress(double percentage);

// Takes a start time and prints duration to console
void stop_timer(std::chrono::high_resolution_clock::time_point start_time);


#endif MISC_H