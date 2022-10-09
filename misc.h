//
// File: misc.h
// Author: Stanley Goodwin
// Creation Date: 6/23/2022
// Last Modified: 6/23/2022
//
#pragma once
#ifndef MISC_H
#define MISC_H

#include <chrono>


void printProgress(double percentage);
void printTimeElapsed(std::chrono::high_resolution_clock::time_point start_time);


#endif MISC_H