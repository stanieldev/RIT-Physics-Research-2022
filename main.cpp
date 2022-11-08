/*
 * File:	main.cpp
 * Author:	Stanley Goodwin
 * The main file to do the simulation execution.
 */
#include "droplet.hpp"



#include "mesh.h"
#include "substrate.h"


/* Runs the droplet simulation.
 * @return 0 on success, 1 on failure.
 */
int main()
{
    Droplet expected_droplet(
        1.9407025E-3,  // Contact radius
        3.0000000E-9   // Expected volume
    );








    Substrate surface(
        (77.0 / 180.0) * PI,         // Receding contact angle of the printed regions
        0.5 * 0.001 / 1.9407025E-3,  // Width of the printed regions
        (5.0 / 180.0) * PI,          // Receding contact angle of the gap regions
        0.5 * 0.001 / 1.9407025E-3   // Width of the gap regions
    );

    // Create mesh using properties above
    const int node_resolution = 121;      // The resolution of the surface's nodes
    const int node_print_interval = 100;  // The interval of the # of iterations before an intermediate print
    Mesh mesh(node_resolution, expected_droplet, surface, node_print_interval);

    // Initialize mesh nodes
    const double θi = (85.0 / 180.0) * PI;  // Spherical cap initial contact angle
    mesh.initialize(θi);

    // Iterate nodes through changing volume
    for (int i = 1; i <= 4; i++)
    {
        mesh.iterate(1000 / i);  // Forward approximate the surface
        mesh.rescale(0.99);      // Rescale surface to enclose new volume
    }

    return 0;
}