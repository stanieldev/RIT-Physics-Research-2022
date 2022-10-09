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
#include <string>
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


/*
 * Prints the duration between a start time and the execution of the function.
 * https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
 *
 * @brief	Duration function (print to console).
 * @param	start_time	time	The time when the clock started ticking.
 */
std::string _duration_string(time start, time stop)
{
    // Namespace 
    using namespace std::chrono;

    // Casting and reducing until variables are made
    auto µ  = duration_cast<microseconds>(stop - start);
    auto ms = duration_cast<milliseconds>(µ);   µ -= duration_cast<microseconds>(ms);
    auto s  = duration_cast<seconds>(ms);      ms -= duration_cast<milliseconds>(s);
    auto m  = duration_cast<minutes>(s);        s -= duration_cast<seconds>(m);
    auto h  = duration_cast<hours>(m);          m -= duration_cast<minutes>(h);

    // Variable
    std::string output_string = "";

    // Collapse time to relevant order
    if (h.count() != 0) {
        output_string.append(std::to_string(h.count()) + "h, ");
        output_string.append(std::to_string(m.count()) + "m");
    }
    else if (m.count() != 0) { 
        output_string.append(std::to_string(m.count()) + "m, ");
        output_string.append(std::to_string(s.count()) + "s");
    }
    else if ( s.count() != 0) {
        output_string.append(std::to_string(s.count()) + ".");
        if (ms.count() < 10) {
            output_string.append("00" + std::to_string(ms.count()));
        }
        else if (ms.count() < 100) {
            output_string.append("0" + std::to_string(ms.count()));
        }
        else {
            output_string.append(std::to_string(ms.count()));
        }
        output_string.append("s");
    }
    else if (ms.count() != 0) {
        output_string.append(std::to_string(ms.count()) + ".");
        if (µ.count() < 10) {
            output_string.append("00" + std::to_string(µ.count()));
        }
        else if (µ.count() < 100) {
            output_string.append("0" + std::to_string(µ.count()));
        }
        else {
            output_string.append(std::to_string(µ.count()));
        }
        output_string.append("ms");
    }
    else { 
        output_string.append(std::to_string(µ.count()) + "us");
    }

    // Return string representation
    return output_string;
}
