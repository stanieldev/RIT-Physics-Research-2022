/**
 * @file	timing.cpp
 * @brief	Duration function.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/23/2022
 * Last Modified: 7/31/2022
 */

#include "timing.h"


/*
 * Returns the duration between a start time and stop time in a neat fashion.
 *
 * @brief	Duration string generator.
 * @param	start	time	The time when the clock started ticking.
 * @param	stop	time	The time when the clock stopped ticking.
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
