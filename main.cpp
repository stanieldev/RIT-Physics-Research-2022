/*
 * File:	main.cpp
 * Author:	Stanley Goodwin
 * The main file to do the simulation execution.
 */
#include "substrate.hpp"
#include "droplet.hpp"
#include "mesh.hpp"


/* 
 * Runs the droplet simulation.
 * @return 0 on success, 1 on failure.
 */
int main()
{
    // Create Droplet characteristic
    Droplet droplet(
        1.9407025E-3,  // Contact radius
        3.0000000E-9   // Expected volume
    );

    // Create Substrate characteristic
    Substrate surface(
        77.0, 0.5 * 0.001 / droplet.contact_radius, // Printed region    [ Receding contact angle (degrees), Width of Region ]
         5.0, 0.5 * 0.001 / droplet.contact_radius  // Unprinted region  [ Receding contact angle (degrees), Width of Region ]
    );

    // Create mesh using properties above
    const int node_resolution = 121;  // The resolution of the surface's nodes
    Mesh mesh(node_resolution, droplet, surface);

    // Mesh settings
    mesh.print_nodes_bool = true;
    mesh.output_folder = "";

    // Mesh iterations
    mesh.initialize(85.0);  // In degrees
    mesh.iterate(1000);     // # of iterations to do
    droplet.volume_ratio *= 0.97;  // Reduce volume ratio (simulates evaporation)

    return 0;
}