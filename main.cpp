/**
 * @file	main.cpp
 * @brief	The main simulation file/function.
 * https://developer.lsst.io/v/DM-5063/docs/cpp_docs.html#file-description-comment-for-header-files
 *
 * @author	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 9/17/2021
 * Last Modified: 7/5/2022
 */
#include <iostream>
#include "mesh.h"
#include "fmath.h"


/** Runs the droplet simulation.
 *
 * @param argc Arguments count
 * @param argv Arguments vector
 * @return 0 on success, something else on failure
*/
int main(int argc, char* argv[])
{
    // Create mesh
    Mesh mesh;

    /*double new_array[6][6] = {
        { 6.49750, 0.50420, 0.18891, 0.52427, 7.61324, 3.81469 },
        { 4.23276, 7.66484, 5.89937, 5.68359, 0.51517, 2.60055 },
        { 7.85275, 0.93038, 6.00715, 0.70866, 1.00814, 1.10923 },
        { 0.53463, 1.01096, 0.83429, 0.00520, 1.96005, 5.80638 },
        { 2.54015, 0.17862, 1.47831, 4.16291, 5.61333, 1.55178 },
        { 4.05010, 5.78871, 4.34007, 2.69967, 2.96173, 7.37391 }
    };

    double output_array[6] = { 1.07190, 6.42021, 4.07980, 5.94209, 8.59383, 1.76560 };
    double coeff[6];


    cramer(new_array, output_array, coeff);
    for (int i = 0; i < 6; i++) {
        std::cout << coeff[i] << std::endl;
    }*/

    /*
    x1 = 7.484036599381260388
    x2 = -6.1349221986463130493
    x3 = -9.0061845204630784139
    x4 = 11.353647627736594218
    x5 = -9.3211120335752221059
    x6 = 5.8328316121470346201
    */
    


    // Mesh functions
    mesh.initialize();
    mesh.iterate();

    return 0;
}
