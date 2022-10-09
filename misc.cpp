/**
 * @file	misc.cpp
 * @brief	Random print functions for timing things.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/23/2022
 * Last Modified: 7/5/2022
 */
#include <iostream>
#include <chrono>
#include "mesh.h"

#define hrc std::chrono::high_resolution_clock
#define time hrc::time_point

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 50


/**
 * Prints a progress bar to the console.
 * https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
 * 
 * @brief	Progress bar function (print to console).
 * @param	percentage	double	The number to find the square root of.
 */
void print_progress(double percentage)
{
    int val  = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


/**
 * Prints the duration between a start time and the execution of the function.
 * https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
 *
 * @brief	Duration function (print to console).
 * @param	start_time	time	The time when the clock started ticking.
 */
void print_duration(time start_time)
{
    // Namespace 
    using namespace std::chrono;

    // The time the function is called (effectively the end-time)
    time stop_time = high_resolution_clock::now();

    // Casting and reducing until variables are made
    auto µ  = duration_cast<microseconds>(stop_time - start_time);
    auto ms = duration_cast<milliseconds>(µ);   µ -= duration_cast<microseconds>(ms);
    auto s  = duration_cast<seconds>(ms);      ms -= duration_cast<milliseconds>(s);
    auto m  = duration_cast<minutes>(s);        s -= duration_cast<seconds>(m);
    auto h  = duration_cast<hours>(m);          m -= duration_cast<minutes>(h);

    // Console output
    std::cout << "Complete! (";
         if ( h.count() != 0) { std::cout << h.count() << "h, " << m.count() <<  "m"; }
    else if ( m.count() != 0) { std::cout << m.count() << "m, " << s.count() <<  "s"; }
    else if ( s.count() != 0) { std::cout << s.count() + ms.count() / 1000.0 <<  "s"; }
    else if (ms.count() != 0) { std::cout << ms.count() + µ.count() / 1000.0 << "ms"; }
    else { std::cout << µ.count() << "us"; }
    std::cout << ")\n";
}
