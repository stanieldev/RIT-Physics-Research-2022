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
 * Last Modified: 7/18/2022
 */
#include "mesh.h"


/* Runs the droplet simulation.
 * @return 0 on success, something else on failure.
 */
int main()
{
    // Create mesh
    Mesh mesh;

    // Mesh functions
    mesh.initialize();
    mesh.iterate(2500);
    printf("%f, %f\n", mesh.volume(), mesh.drop_κ);
    mesh.iterate(25000 - 2500);
    printf("%f, %f\n", mesh.volume(), mesh.drop_κ);
    mesh.iterate(25000);
    printf("%f, %f\n", mesh.volume(), mesh.drop_κ);

    // Return 0 if successful
    return 0;
}
