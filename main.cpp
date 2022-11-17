/*
 * File:	main.cpp
 * Author:	Stanley Goodwin
 * The main file to do the simulation execution.
 */
#include "substrate.hpp"
#include "droplet.hpp"
#include "mesh.hpp"


/* Runs the droplet simulation.
 * @return 0 on success, 1 on failure.
 */
int main()
{
    // User inputs
    const int node_resolution = 121;      // The resolution of the surface's nodes

    Droplet droplet(
        1.9407025E-3,  // Contact radius
        3.0000000E-9   // Expected volume
    );

    const double printed_width = 0.5 * 0.001 / droplet.contact_radius;
    const double unprinted_width = 0.5 * 0.001 / droplet.contact_radius;
    Substrate surface(
        77.0, printed_width,  // Printed region  [ Receding contact angle (degrees), Width of Region ]
        5.0,  unprinted_width // Unprinted region  [ Receding contact angle (degrees), Width of Region ]
    );


    // Create mesh using properties above
    Mesh mesh(node_resolution, droplet, surface);


    // Initialize mesh nodes
    const double θi_degrees = 85.0;  // Spherical cap initial contact angle
    mesh.initialize(θi_degrees);

    // Iterate 1000 steps
    mesh.iterate(1000);


    return 0;
}