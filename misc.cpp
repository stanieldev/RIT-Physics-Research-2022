//
// File: misc.cpp
// Author: Stanley Goodwin
// Creation Date: 6/23/2022
// Last Modified: 6/23/2022
//
#include <iostream>
#include <chrono>

using namespace std::chrono; // Namespace for simplification


// Progress bar function
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
void printProgress(double percentage) { // https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


// Print-duration-of-function-execution function
void printTimeElapsed(high_resolution_clock::time_point start_time)
{
    // The time the function is called (effectively the end-time)
    auto stop_time = high_resolution_clock::now();

    // Casting and reducing until variables are made
    auto µ = duration_cast<microseconds>(stop_time - start_time);
    auto ms = duration_cast<milliseconds>(µ);
    µ -= duration_cast<microseconds>(ms);
    auto s = duration_cast<seconds>(ms);
    ms -= duration_cast<milliseconds>(s);
    auto m = duration_cast<minutes>(s);
    s -= duration_cast<seconds>(m);
    auto h = duration_cast<hours>(m);
    m -= duration_cast<minutes>(h);

    // Console output
    std::cout << "Complete! (";
    if (h.count() != 0) {
        std::cout << h.count() << "h, " << m.count() << "m";
    }
    else if (m.count() != 0) {
        std::cout << m.count() << "m, " << s.count() << "s";
    }
    else if (s.count() != 0) {
        std::cout << s.count() + ms.count() / 1000.0 << "s";
    }
    else if (ms.count() != 0) {
        std::cout << ms.count() + µ.count() / 1000.0 << "ms";
    }
    else {
        std::cout << µ.count() << "us";
    }
    std::cout << ")\n";
}
