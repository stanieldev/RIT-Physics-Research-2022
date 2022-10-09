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

    // Mesh functions
    mesh.InitializeNodes();
    mesh.Iterate();

    return 0;
}
