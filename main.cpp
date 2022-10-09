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
#include "new_droplet.h"
#include "new_substrate.h"


/* Runs the droplet simulation.
 * @return 0 on success, something else on failure.
 */
int main()
{
    // Create mesh
    Droplet droplet(1.9407025E-3, 3.0E-9);
    Mesh mesh;
    
    Substrate surface(
        (77.0 / 180.0) * PI, 0.5 * 0.001 / 1.9407025E-3,
        (45.0 / 180.0) * PI, 0.5 * 0.001 / 1.9407025E-3
    );
    



    // Initialize mesh nodes
    const double θi = (45.0 / 180.0) * PI;  // Spherical cap initial contact angle [To be depreciated]
    mesh.initialize(θi);

    // Iterate nodes
    mesh.iterate(5000, surface);

    // Iterate nodes
    /*mesh.drop_κ *= 0.97;
    mesh.iterate(1000);*/


    // Return 0 if successful
    return 0;
}