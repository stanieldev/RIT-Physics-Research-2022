/**
 * File:	main.cpp
 * Brief:	The execution file.
 *
 * Authors:	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 9/17/2021
 * Last Modified: 7/31/2022
 */
#include "mesh.h"
#include "droplet.h"
#include "substrate.h"


/* Runs the droplet simulation.
 * @return 0 on success, something else on failure.
 */
int main()
{
    // Characteristics of the droplet
    Droplet droplet(
        1.9407025E-3,  // Expected contact radius
        3.0E-9         // Expected final volume
    );

    // Characteristics of the substrate
    Substrate surface(
        (77.0 / 180.0) * PI,         // Receding contact angle of the printed regions
        0.5 * 0.001 / 1.9407025E-3,  // Width of the printed regions
        (5.0 / 180.0) * PI,          // Receding contact angle of the gap regions
        0.5 * 0.001 / 1.9407025E-3   // Width of the gap regions
    );

    // Create mesh using properties above
    int node_print_interval = 100;  // The interval of the # of iterations before an intermediate print
    Mesh mesh(121, droplet, surface, node_print_interval);  // Resolution, Droplet, Substrate

    
    // Initialize mesh nodes
    double θi = (85.0 / 180.0) * PI;  // Spherical cap initial contact angle
    mesh.initialize(θi);

    // Iterate nodes
    for (int i = 1; i <= 4; i++)
    {
        mesh.iterate(1000 / i);
        mesh.rescale(0.99);
    }


    // Return 0 if successful
    return 0;
}